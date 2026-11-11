"""Every experiment of the evaluation, in order. A step saves results/<step>.json and is skipped
when that file exists; the external tools also cache each program in .cache/.

    uv run python bench/run.py && uv run python bench/tables.py
"""

import json
import math
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(BENCH))

import sample  # noqa: E402
import speed  # noqa: E402
import variants  # noqa: E402
from milo.analysis import BOUNDS, FAIL, set_format  # noqa: E402
from milo.egraph import EGGLOG, build, guarded  # noqa: E402
from milo.extract import best, extract  # noqa: E402
from milo.fpcore import expr, fpcore, parse, read, show  # noqa: E402
from tools import daisy, fptaylor, gappa, herbie, precisa  # noqa: E402

RESULTS = os.path.join(BENCH, "results")
JOBS = os.cpu_count()
EGGLOG_LIMIT, EXTRACT_LIMIT, TOOL_LIMIT = 300, 120, 300
HEADROOM = 3
SIZES = (1, 10, 100)


def cores(directory):
    return sorted(f"{directory}/{f}" for f in os.listdir(os.path.join(BENCH, directory)) if f.endswith(".fpcore"))


def load(name):
    with open(os.path.join(BENCH, name)) as f:
        core = parse(f.read())
    set_format(core.precision)
    return core


def program(text):
    return expr(read(text)[0])


def step(f):
    loaded = []

    def cached():
        path = os.path.join(RESULTS, f.__name__ + ".json")
        if not os.path.exists(path):
            print(f.__name__, flush=True)
            result = f()
            os.makedirs(RESULTS, exist_ok=True)
            with open(path, "w") as out:
                json.dump(result, out, indent=1)
        if not loaded:
            with open(path) as saved:
                loaded.append(json.load(saved))
        return loaded[0]
    return cached


def each(f, *args):
    """f over the zipped args in parallel, keyed by the first; an exception becomes {"error"}."""
    with ProcessPoolExecutor(JOBS) as pool:
        return dict(zip(args[0], pool.map(safe, [f] * len(args[0]), *args)))


def safe(f, *args):
    start = time.time()
    try:
        return f(*args)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"[:300], "seconds": time.time() - start}


def bound(e):
    return {m: math.inf if e is FAIL else float(b(e)) for m, b in BOUNDS.items()}


def metric(bounds):
    """The metric a program is compared on, from the best bound on each: relative error where it is
    below 1, since a relative bound of 1 or more says nothing about the error's size."""
    return "rel" if bounds["rel"] < 1 else "abs"


def pick(bounds):
    """The metric our program is picked by, from the best bound on each: relative error wherever it is
    finite. Even a relative bound of 1 or more makes the error shrink with the real value near its zeros,
    which no absolute bound does."""
    return "rel" if math.isfinite(bounds["rel"]) else "abs"


def readout(m):
    """The m-bound of a cost, with ties broken by the bound in the other metric."""
    other = "abs" if m == "rel" else "rel"
    return lambda e: (BOUNDS[m](e), BOUNDS[other](e))


def ours_metric(o):
    """metric() of our program on a core, which is the input when we found none."""
    return metric(o["bounds"] if o.get("programs") else o.get("original") or o.get("bound") or bound(FAIL))


def failure(e):
    return "timeout" if isinstance(e, subprocess.TimeoutExpired) else f"{type(e).__name__}: {e}"[:300]


def bound_input(core, body):
    """Our bound on a program as it is: no rewriting, no new branches."""
    g = build(core._replace(body=body), 0, False)
    F, _, _ = extract(g)
    return bound(F[g.root][0][0] if F[g.root] else FAIL)


def ours_one(name):
    """Our bound on the input, then, for a rewriting core, our rewrite of it. When the rewrite fails,
    the input's bound survives and stands for the input's bound in our e-graph."""
    core = load(name)
    start = time.time()
    out = {"bound": bound_input(core, core.body), "bound_s": time.time() - start}
    if not name.startswith("rewriting_cores/"):
        return out
    start = time.time()
    try:
        g = build(core, 4, False, EGGLOG_LIMIT)
        built = time.time()
        F, k, cut = extract(g, EXTRACT_LIMIT)
        chosen = {m: best(F[g.root], readout(m)) for m in BOUNDS} if F[g.root] else {}
        return dict(out, egglog_s=built - start, extract_s=time.time() - built, rounds=k, cut=cut,
                    nodes=sum(map(len, g.nodes.values())), classes=len(g.nodes), root_frontier=len(F[g.root]),
                    max_frontier=max(map(len, F.values())),
                    original=bound(variants.analyze(g, guarded(core.body, core.box, core.precision))),
                    programs={m: show(p) for m, (e, p) in chosen.items()},
                    bounds={m: bound(e)[m] for m, (e, p) in chosen.items()} or bound(FAIL))
    except Exception as e:
        return dict(out, egglog_s=time.time() - start, extract_s=0, failed=failure(e), original=out["bound"],
                    programs={}, bounds=bound(FAIL))


@step
def ours():
    """Benchmark (1): our bound on each input, and our rewrite without branches of Benchmark (2)'s."""
    return each(ours_one, benchmark1())


def original_one(name):
    core = load(name)
    return sample.measure(core, core.body, {"original": core.body})["original"]


@step
def originals():
    """Benchmark (1): the measured error of each input."""
    return each(original_one, benchmark1())


def benchmark1():
    """Every core: the ones we only bound, and the ones we also rewrite."""
    return sorted(cores("bounding_cores") + cores("rewriting_cores"), key=os.path.basename)


def benchmark2():
    """The rewriting cores: those whose input loses more than HEADROOM ulps somewhere, which the
    measured error of the inputs confirms."""
    names = cores("rewriting_cores")
    lossy = {n for n, r in originals().items() if (r.get("ulps") or 0) > HEADROOM}
    if lossy != set(names):
        print(f"warning: rewriting_cores/ is not the cores above {HEADROOM} ulps: {sorted(lossy ^ set(names))}")
    return names


@step
def rewrites():
    """Benchmark (2): Herbie's and Daisy's rewrites, without branches."""
    texts = {n: fpcore(load(n)) for n in benchmark2()}
    with tempfile.TemporaryDirectory() as d:
        found = {"herbie": herbie.rewrite(list(texts.values()), d, TOOL_LIMIT, JOBS, branch=False),
                 "daisy": daisy.rewrite(list(texts.values()), d, TOOL_LIMIT, JOBS)}
    return {n: {tool: got[t] for tool, got in found.items()} for n, t in texts.items()}


def variants_one(name):
    core = load(name)
    start = time.time()
    g = build(core, 4, False, EGGLOG_LIMIT)
    built = time.time() - start

    def chosen(p):
        return {"program": show(p), **bound(variants.analyze(g, p))} if p else bound(FAIL)

    out = {}
    for n in SIZES:
        start = time.time()
        found = [(p, bound(variants.analyze(g, p))) for p in variants.smallest(g, n, EXTRACT_LIMIT)]
        out[f"size{n}"] = {m: chosen(min(found, key=lambda f: f[1][m])[0] if found else None) for m in BOUNDS}
        out[f"size{n}"]["seconds"] = built + time.time() - start
    start = time.time()
    out["greedy"] = {}
    for m in BOUNDS:
        F, _, _ = extract(g, EXTRACT_LIMIT, leq=variants.greedy(m))
        out["greedy"][m] = chosen(F[g.root][0][1] if F[g.root] else None)
    out["greedy"]["seconds"] = built + time.time() - start
    for m in BOUNDS:
        start = time.time()
        p = chosen(variants.single(g, m, EXTRACT_LIMIT))
        out[f"{m}-only"] = {"abs": p, "rel": p, "seconds": built + time.time() - start}
    return out


@step
def alternatives():
    """Benchmark (2): extraction by size, greedy readout, and single-metric analyses."""
    return each(variants_one, benchmark2())


def candidates(name):
    """The rewrites of a Benchmark (2) core by our tool, Herbie and Daisy, by who produced them."""
    o, r = ours()[name], rewrites()[name]
    out = {"ours": o.get("programs", {}).get(pick(o["bounds"])), **{tool: r[tool].get("program") for tool in r}}
    return {who: p for who, p in out.items() if p}


def measured_one(name, texts):
    core = load(name)
    return sample.measure(core, core.body, {"original": core.body, **{w: program(t) for w, t in texts.items()}})


@step
def measured():
    """Benchmark (2): the measured error of every compared program."""
    names = benchmark2()
    return each(measured_one, names, [candidates(n) for n in names])


@step
def bounds():
    """Sound bounds from FPTaylor, Gappa, PRECiSA and Daisy: on every input, and on Benchmark (2)'s rewrites."""
    texts, b2 = {}, set(benchmark2())
    for n in benchmark1():
        core = load(n)
        texts[n] = {"original": fpcore(core)}
        if n in b2:
            texts[n].update({w: fpcore(core, program(t)) for w, t in candidates(n).items()})
    every = [t for ts in texts.values() for t in ts.values()]
    with tempfile.TemporaryDirectory() as d:
        found = {"fptaylor": fptaylor.bound(every, d, TOOL_LIMIT, JOBS),
                 "gappa": gappa.bound(every, d, TOOL_LIMIT, JOBS),
                 "precisa": precisa.bound(every, d, TOOL_LIMIT, JOBS),
                 "daisy": daisy.bound([ts["original"] for ts in texts.values()], d, TOOL_LIMIT, JOBS)}
    return {n: {tool: {w: got[t] for w, t in ts.items() if t in got} for tool, got in found.items()}
            for n, ts in texts.items()}


@step
def startup():
    """Each tool's startup, left out of its running times: for the tools launched once per program, the
    least time of a few runs, one at a time, on a trivial program."""
    texts = [f'(FPCore (x) :name "startup {i}" :precision binary64 :pre (<= 1 x 2) (+ x 1))' for i in range(5)]
    with tempfile.TemporaryDirectory() as d:
        found = {"fptaylor": fptaylor.bound(texts, d, TOOL_LIMIT, 1), "gappa": gappa.bound(texts, d, TOOL_LIMIT, 1),
                 "precisa": precisa.bound(texts, d, TOOL_LIMIT, 1), "daisy": daisy.bound(texts, d, TOOL_LIMIT, 1),
                 "daisy-rewrite": daisy.rewrite(texts, d, TOOL_LIMIT, 1)}
        out = {tool: min(got[t]["seconds"] for t in texts) for tool, got in found.items()}
        empty = os.path.join(d, "empty.egg")
        open(empty, "w").close()
        egglog = []
        for _ in texts:
            start = time.time()
            subprocess.run([EGGLOG, empty], capture_output=True)
            egglog.append(time.time() - start)
    return dict(out, egglog=min(egglog))


def speed_one(name, texts):
    core = load(name)
    programs = {"original": core.body, **{w: program(t) for w, t in texts.items()}}
    return speed.measure(core, programs, sample.points(core, 1000))


@step
def speeds():
    """Benchmark (2): run time of the input and of each rewrite, one core at a time."""
    return {n: safe(speed_one, n, candidates(n))
            for n in benchmark2()}


def branching_one(name, reports):
    core = load(name)
    programs = {"original": core.body}
    out = {"expected": core.props.get(":milo-expected", "").strip('"')}
    for who, branch in (("ours-plain", False), ("ours", True)):
        try:
            g = build(core, 4, branch, EGGLOG_LIMIT)
            F, _, _ = extract(g, EXTRACT_LIMIT)
        except Exception as e:
            out[who] = {"status": failure(e)}
            continue
        if not F[g.root]:
            out[who] = {"status": "no program"}
            continue
        chosen = {m: best(F[g.root], readout(m)) for m in BOUNDS}
        programs[who] = chosen[pick({m: bound(e)[m] for m, (e, p) in chosen.items()})][1]
    for who, report in reports.items():
        if report.get("program"):
            programs[who] = program(report["program"])
        else:
            out[who] = {"status": report.get("status", "fail")}
    measured = sample.measure(core, core.body, programs)
    for who, p in programs.items():
        out[who] = {"program": show(p), "measured": measured[who], "bound": safe(bound_input, core, p)}
    return out


@step
def branching():
    """Benchmark (3): our branches against Herbie's and Daisy's regimes, timed one core at a time."""
    names = cores("branching")
    texts = [fpcore(load(n)) for n in names]
    with tempfile.TemporaryDirectory() as d:
        found = {"herbie": herbie.rewrite(texts, d, TOOL_LIMIT, JOBS, branch=True),
                 "daisy": daisy.rewrite(texts, d, TOOL_LIMIT, JOBS, regimes=True)}
    out = each(branching_one, names, [{who: f[t] for who, f in found.items()} for t in texts])
    for name, r in out.items():
        whos = [w for w in r if isinstance(r[w], dict) and "program" in r[w]]
        if whos:
            times = speed.measure(load(name), {w: program(r[w]["program"]) for w in whos},
                                  sample.points(load(name), 1000))
            for w in whos:
                r[w]["ns"] = times[w]
    return out


if __name__ == "__main__":
    for f in (ours, originals, rewrites, alternatives, measured, bounds, startup, speeds, branching):
        f()
