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

import defined  # noqa: E402
import sample  # noqa: E402
import speed  # noqa: E402
import variants  # noqa: E402
from milo.analysis import BOUNDS, FAIL, set_format  # noqa: E402
from milo.egraph import EGGLOG, build, guarded  # noqa: E402
from milo.extract import analyze, best, extract  # noqa: E402
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
    """The metric our program is picked by: relative error wherever its bound is finite, since even a
    bound of 1 or more makes the error shrink with the value near its zeros."""
    return "rel" if math.isfinite(bounds["rel"]) else "abs"


def readout(m):
    """The m-bound of a cost, with ties broken by the bound in the other metric."""
    other = "abs" if m == "rel" else "rel"
    return lambda e: (BOUNDS[m](e), BOUNDS[other](e))


def ours_metric(o):
    """metric() of our program on a core, or of the input when we found none."""
    return metric(o["bounds"] if o["programs"] else o["original"])


def failure(e):
    return "timeout" if isinstance(e, subprocess.TimeoutExpired) else f"{type(e).__name__}: {e}"[:300]


def bound_input(core, body):
    """Our bound on a program as it is: no rewriting, no new branches."""
    g = build(core._replace(body=body), 0, False)
    F, _, _ = extract(g)
    return bound(F[g.root][0][0] if F[g.root] else FAIL)


def ours_one(name):
    """Our rewrite without branches; original is the input's bound in the e-graph, or on its own if
    that fails. cut is whether extraction stopped before its fixed point."""
    core = load(name)
    start = time.time()
    try:
        g = build(core, 4, False, EGGLOG_LIMIT)
        built = time.time()
        F, _, cut = extract(g, EXTRACT_LIMIT)
        chosen = {m: best(F[g.root], readout(m)) for m in BOUNDS} if F[g.root] else {}
        return dict(egglog_s=built - start, extract_s=time.time() - built, cut=cut,
                    original=bound(analyze(g, guarded(core.body, core.box, core.precision))),
                    programs={m: show(p) for m, (e, p) in chosen.items()},
                    bounds={m: bound(e)[m] for m, (e, p) in chosen.items()} or bound(FAIL))
    except Exception as e:
        return dict(egglog_s=time.time() - start, extract_s=0, failed=failure(e),
                    original=bound_input(core, core.body), programs={}, bounds=bound(FAIL))


@step
def ours():
    """Benchmark (2): our rewrite of each input, without branches."""
    return each(ours_one, benchmark2())


def original_one(name):
    core = load(name)
    return sample.measure(core, core.body, {"original": core.body})["original"]


def bound_only_one(name):
    """milo --bound-only: the input's bound in its e-graph, built as usual (4 rounds, branching on)."""
    core = load(name)
    start = time.time()
    g = build(core, 4, True, EGGLOG_LIMIT)
    out = bound(analyze(g, guarded(core.body, core.box, core.precision)))
    return dict(out, seconds=time.time() - start)


@step
def bound_only():
    """Benchmark (1): our bound on each input in the e-graph, as milo --bound-only gives it."""
    return each(bound_only_one, benchmark1())


@step
def originals():
    """Benchmark (1): the measured error of each input."""
    return each(original_one, benchmark1())


def benchmark1():
    """Every core: the ones we only bound, and the ones we also rewrite."""
    return sorted(cores("bounding_programs") + cores("rewriting_programs"), key=os.path.basename)


def benchmark2():
    """The rewriting cores: those whose input loses more than HEADROOM ulps somewhere, which the
    measured error of the inputs confirms."""
    names = cores("rewriting_programs")
    lossy = {n for n, r in originals().items() if (r.get("ulps") or 0) > HEADROOM}
    if lossy != set(names):
        print(f"warning: rewriting_programs/ is not the programs above {HEADROOM} ulps: {sorted(lossy ^ set(names))}")
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
        return {"program": show(p), **bound(analyze(g, p))} if p else bound(FAIL)

    out = {}
    for n in SIZES:
        start = time.time()
        found = [(p, bound(analyze(g, p))) for p in variants.smallest(g, n, EXTRACT_LIMIT)]
        out[f"size{n}"] = {m: chosen(min(found, key=lambda f: f[1][m])[0] if found else None) for m in BOUNDS}
        out[f"size{n}"]["seconds"] = built + time.time() - start
    return out


@step
def alternatives():
    """Benchmark (2): the best of the smallest programs, for each size in SIZES."""
    return each(variants_one, benchmark2())


def single_one(name):
    core = load(name)
    start = time.time()
    g = build(core, 4, False, EGGLOG_LIMIT)
    built = time.time() - start
    out = {}
    for m in BOUNDS:
        start = time.time()
        p, own = variants.single_bound(g, m, EXTRACT_LIMIT)
        out[m] = {"program": show(p) if p else None, "own": own, "seconds": built + time.time() - start}
    return out


@step
def single_bounds():
    """Benchmark (2): the program extracted by each single-metric analysis, with its bound there."""
    return each(single_one, benchmark2())


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


def defined_one(name, texts):
    core = load(name)
    programs = {"original": core.body, **{w: program(t) for w, t in texts.items()}}
    return {w: defined.check(p, core.box) for w, p in programs.items()}


@step
def definedness():
    """Benchmark (2): whether each compared program is defined over the reals on all of its box. A
    rewrite that is undefined somewhere there is not equivalent to its input, which is defined."""
    names = benchmark2()
    return each(defined_one, names, [candidates(n) for n in names])


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
    """Each tool's startup, left out of its running times: its least time on a trivial program."""
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
    out, extra = {}, {}
    for who, branch in (("ours-plain", False), ("ours", True)):
        start = time.time()
        try:
            g = build(core, 4, branch, EGGLOG_LIMIT)
            F, _, _ = extract(g, EXTRACT_LIMIT)
        except Exception as e:
            out[who] = {"status": failure(e), "seconds": time.time() - start}
            continue
        if not F[g.root]:
            out[who] = {"status": "no program", "seconds": time.time() - start}
            continue
        chosen = {m: best(F[g.root], readout(m)) for m in BOUNDS}
        bounds = {m: bound(e)[m] for m, (e, p) in chosen.items()}
        programs[who] = chosen[pick(bounds)][1]
        extra[who] = {"seconds": time.time() - start, "bounds": bounds}
    for who, report in reports.items():
        if report.get("program"):
            programs[who] = program(report["program"])
            extra[who] = {"seconds": report.get("seconds")}
        else:
            out[who] = {"status": report.get("status", "fail"), "seconds": report.get("seconds")}
    measured = sample.measure(core, core.body, programs)
    for who, p in programs.items():
        out[who] = {"program": show(p), "measured": measured[who], "bound": safe(bound_input, core, p),
                    "defined": defined.check(p, core.box), **extra.get(who, {})}
    return out


def first_condition(e):
    """The condition of the first if-statement in e, or None if it has none."""
    match e:
        case ("If", c, _, _):
            return c
        case ("Var", _) | ("Num", _):
            return None
    return next((c for a in e[1:] if (c := first_condition(a))), None)


def resolve(e, c, taken):
    """e with every if-statement on condition c replaced by the branch it takes."""
    match e:
        case ("If", d, a, b):
            if d == c:
                return resolve(a if taken else b, c, taken)
            return ("If", d, resolve(a, c, taken), resolve(b, c, taken))
        case ("Var", _) | ("Num", _):
            return e
    return (e[0], *(resolve(a, c, taken) for a in e[1:]))


def side(box, c, taken, fmt):
    """The box cut to the floats that take (or skip) the condition x op k, comparing x with k rounded
    to the format as the program does; None if no float does."""
    op, x, k = c
    kf = sample.rounded(k, fmt)
    above, strict = (op[0] == ">") == taken, (op in (">", "<")) == taken
    lo, hi = box[x]
    if above:
        lo = max(lo, sample.from_ordinal(sample.ordinal(kf, fmt) + 1, fmt) if strict else kf)
    else:
        hi = min(hi, sample.from_ordinal(sample.ordinal(kf, fmt) - 1, fmt) if strict else kf)
    return dict(box, **{x: (lo, hi)}) if lo <= hi else None


def pieces(core, body):
    """(box, straight-line program) for each way through body's if-statements, where box holds the floats
    that go that way. The boxes cover core's box, so the largest bound over the pieces bounds body."""
    c = first_condition(body)
    if c is None:
        return [(core.box, body)]
    out = []
    for taken in (True, False):
        box = side(core.box, c, taken, core.precision)
        if box is not None:
            out += pieces(core._replace(box=box), resolve(body, c, taken))
    return out


@step
def branching():
    """Benchmark (3): our branches against Herbie's and Daisy's regimes, timed one core at a time."""
    names = cores("branching_programs")
    texts = [fpcore(load(n)) for n in names]
    with tempfile.TemporaryDirectory() as d:
        found = {"herbie": herbie.rewrite(texts, d, TOOL_LIMIT, JOBS, branch=True),
                 "daisy": daisy.rewrite(texts, d, TOOL_LIMIT, JOBS, regimes=True)}
    out = each(branching_one, names, [{who: f[t] for who, f in found.items()} for t in texts])
    for name, r in out.items():
        whos = [w for w in r if isinstance(r[w], dict) and "program" in r[w]]
        if whos:
            times = safe(speed.measure, load(name), {w: program(r[w]["program"]) for w in whos},
                         sample.points(load(name), 1000))
            for w in whos:
                if w in times:
                    r[w]["ns"] = times[w]
                else:
                    r[w]["speed_error"] = times.get("error")
    return out


@step
def branching_bounds():
    """Sound bounds on Benchmark (3)'s programs from FPTaylor, Gappa and PRECiSA, which take no
    if-statements: each program's pieces(), bounded one at a time, as a list per tool."""
    texts = {}
    for name, r in branching().items():
        core = load(name)
        for who, p in r.items():
            if "program" in p:
                texts[name, who] = [fpcore(core._replace(box=box, name=f"{core.name} {i}"), body)
                                    for i, (box, body) in enumerate(pieces(core, program(p["program"])))]
    every = sorted({t for ts in texts.values() for t in ts})
    with tempfile.TemporaryDirectory() as d:
        found = {"fptaylor": fptaylor.bound(every, d, TOOL_LIMIT, JOBS), "gappa": gappa.bound(every, d, TOOL_LIMIT, JOBS),
                 "precisa": precisa.bound(every, d, TOOL_LIMIT, JOBS)}
    out = {}
    for (name, who), ts in texts.items():
        out.setdefault(name, {})[who] = {tool: [got[t] for t in ts] for tool, got in found.items()}
    return out


if __name__ == "__main__":
    for f in (originals, ours, bound_only, rewrites, alternatives, measured, definedness, bounds, startup, speeds,
              branching, branching_bounds):
        f()
