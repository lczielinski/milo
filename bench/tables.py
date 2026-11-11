"""The tables of the evaluation, from results/*.json into results/tables/<name>.tex. A tool that fails
counts as returning its input, and each table says how many did. A relative bound of 1 or more says
nothing, so it counts as no bound. Each table compares a row with a baseline, where lower is better."""

import math
import os
import statistics

import run
from run import SIZES

INF = math.inf
LN10 = math.log(10)
NAMES = {"fptaylor": "FPTaylor", "gappa": "Gappa", "precisa": "PRECiSA", "daisy": "Daisy", "herbie": "Herbie",
         "ours": "Ours", "ours-plain": "No branch", "greedy": "Ours, greedy", "abs-only": "Absolute error only",
         "rel-only": "Relative error only"}
SOUND = ("fptaylor", "gappa", "precisa")
BSW = ("b/s/w counts the programs where {who} is better than {base} by at least one bit (a factor of 2), within one "
       "bit of it, or worse by at least one bit")


def num(x):
    if x is None:
        return "--"
    return str(int(x)) if float(x).is_integer() and abs(x) < 1e6 else f"{x:.3g}"


def lognum(x):
    """A ratio given by its log, written out even past the float range."""
    if x is None:
        return "--"
    if math.isinf(x):
        return "inf" if x > 0 else "0"
    e = math.floor(x / LN10)
    return num(math.exp(x)) if abs(e) < 300 else f"{math.exp(x - e * LN10):.3g}e{e:+d}"


def median(xs):
    return statistics.median(xs) if xs else None


def geomean(xs):
    return math.exp(statistics.fmean(map(math.log, xs))) if xs else None


def said(m, x):
    """A bound, or INF for a missing one and for a relative bound of 1 or more."""
    return INF if x is None or (m == "rel" and x >= 1) else x


def change(a, b):
    """0 if a is at least 2x below b, 2 if at least 2x above, and 1 between, where INF is above and 0
    below any other value."""
    if a == b:
        return 1
    if min(a, b) == 0 or max(a, b) == INF:
        return 0 if a < b else 2
    return 0 if a <= b / 2 else 2 if a >= 2 * b else 1


def quotient(a, b):
    """log(a / b), where 0 is below and INF above any other value, and equal values give 0."""
    if a == b:
        return 0.0
    if min(a, b) == 0 or max(a, b) == INF:
        return -INF if a < b else INF
    return math.log(a) - math.log(b)


def middle(logs):
    """The median of log ratios that may be infinite: for an even count, the mean of the two middle
    ones when both are finite, and the lower one otherwise."""
    xs = sorted(logs)
    if not xs:
        return None
    lo, hi = xs[(len(xs) - 1) // 2], xs[len(xs) // 2]
    return (lo + hi) / 2 if math.isfinite(lo) and math.isfinite(hi) else lo


def compare(pairs):
    """For pairs (value, baseline), where lower is better: how many values are better, the same and
    worse by change(), and the median of value / baseline."""
    pairs = list(pairs)
    found = [change(a, b) for a, b in pairs]
    return ["/".join(str(found.count(k)) for k in range(3)), lognum(middle([quotient(a, b) for a, b in pairs]))]


def table(caption, label, columns, head, rows):
    """A row is a list of cells, or a string such as \\midrule."""
    body = "\n".join(r if isinstance(r, str) else " & ".join(map(str, r)) + " \\\\" for r in rows)
    return (f"\\begin{{table}}[H] \\centering \\small\n\\caption{{{caption}}}\\label{{{label}}}\n"
            f"\\begin{{tabular}}{{{columns}}}\n\\toprule\n{head} \\\\\n\\midrule\n{body}\n"
            "\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def metric(name):
    return run.ours_metric(run.ours()[name])


def sound(name, who):
    """The smallest of FPTaylor's, Gappa's and PRECiSA's bounds on who's program, on the core's metric."""
    m = metric(name)
    return min(said(m, run.bounds()[name][t].get(who, {}).get(m)) for t in SOUND)


def accuracy():
    measured, speeds = run.measured(), run.speeds()
    names = [n for n in run.benchmark2() if "error" not in measured[n]]
    tools = ("herbie", "daisy", "ours")

    def got(t, n):
        """t's program, or the input if t returned none or one that disagrees with it over the reals."""
        return t if measured[n].get(t, {}).get("equivalent") else "original"

    def speed(t):
        return geomean([speeds[n][got(t, n)] / speeds[n]["original"] for n in names
                        if speeds[n].get(got(t, n)) and speeds[n].get("original")])

    columns = (lambda t, n: (sound(n, got(t, n)), sound(n, "original")),
               lambda t, n: (measured[n][got(t, n)]["ulps"], measured[n]["original"]["ulps"]),
               lambda t, n: (measured[n][got(t, n)]["mean"], measured[n]["original"]["mean"]))
    rows = [[NAMES[t], *(cell for f in columns for cell in compare(f(t, n) for n in names)), num(speed(t)),
             sum(got(t, n) == "original" for n in names)] for t in tools]
    head = ("& \\multicolumn{2}{c}{Sound Bound} & \\multicolumn{2}{c}{Max Sampled Error} & "
            "\\multicolumn{2}{c}{Avg Sampled Error} & & \\\\\n"
            "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n"
            "Tool" + " & b/s/w & Median" * 3 + " & Speed & Failed")
    left = len(run.benchmark2()) - len(names)
    disagree = sum(t in measured[n] and not measured[n][t]["equivalent"] for t in tools for n in names)
    return table(f"Accuracy of the rewrites of the {len(names)} programs of Benchmark (2), against the input. "
                 + BSW.format(who="the rewrite's bound or error", base="the input's")
                 + ", where a missing bound is infinite. Medians are of the rewrite's bound or error divided by the "
                 "input's (lower is better). A tool that returns no program"
                 + ", or one that disagrees with its input over the reals at one of 20 sampled inputs," * bool(disagree)
                 + " counts as returning its input (Failed). Speed is the geomean of run time divided by the input's "
                 "(lower is better)." + f" Programs whose measurement failed are left out ({left})." * bool(left),
                 "tab:accuracy", "lcrcrcrrr", head, rows)


def tightness():
    ours, bounds = run.ours(), run.bounds()
    mine = {m: {n: said(m, ours[n].get("bound", {}).get(m)) for n in ours} for m in ("abs", "rel")}

    def cells(t, m):
        b = mine[m]
        if t == "ours":
            return [sum(x < INF for x in b.values()), "--", "--"]
        a = {n: said(m, bounds[n][t].get("original", {}).get(m)) for n in ours}
        return [sum(x < INF for x in a.values()), *compare((a[n], b[n]) for n in ours if min(a[n], b[n]) < INF)]

    rows = [[NAMES[t], *cells(t, "abs"), *cells(t, "rel")] for t in ("ours", *SOUND, "daisy")]
    head = ("& \\multicolumn{3}{c}{Absolute Error} & \\multicolumn{3}{c}{Relative Error} \\\\\n"
            "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
            "Tool" + " & Bounded & b/s/w & Median" * 2)
    return table(f"Sound bounds on the {len(ours)} programs of Benchmark (1). Bounded counts the programs each tool "
                 "bounds. The other columns cover the programs either tool bounds: "
                 + BSW.format(who="the tool's bound", base="ours")
                 + ", where a missing bound is infinite, and medians are of the tool's bound divided by ours (lower "
                 "is better).", "tab:tightness", "lrcrrcr", head, rows)


def raw(name, who, m):
    """The m-bound of who's program on a Benchmark (2) core, with our analysis in our e-graph."""
    o = run.ours()[name]
    if who == "ours":
        x = o["bounds"][m] if o.get("programs") else None
    elif who == "original":
        x = (o.get("original") or o.get("bound") or {}).get(m)
    else:
        x = run.alternatives()[name].get(who, {}).get(m, {}).get(m)
    return said(m, x)


def picked(name, who, m):
    """who, or the input if who has no program with a finite m-bound."""
    return who if raw(name, who, m) < INF else "original"


def bound(name, who, m):
    return raw(name, picked(name, who, m), m)


def seconds(name, who):
    o, a = run.ours()[name], run.alternatives()[name]
    if who == "ours":
        return o["egglog_s"] + o["extract_s"] if "egglog_s" in o else o.get("seconds", 0)
    return a[who]["seconds"] if who in a else a.get("seconds", 0)


def extraction():
    names = run.benchmark2()
    m = {n: metric(n) for n in names}

    def time(w):
        return num(round(sum(seconds(n, w) for n in names)))

    def failed(w):
        return sum(picked(n, w, m[n]) == "original" for n in names)

    optimal = sum("rounds" in run.ours()[n] and not run.ours()[n]["cut"] for n in names)
    rows = [[NAMES.get(w, f"Size, top-{w[4:]}"), *compare((bound(n, w, m[n]), bound(n, "ours", m[n])) for n in names),
             "--", failed(w), time(w)] for w in [f"size{k}" for k in SIZES] + ["greedy"]]
    rows.append(["Ours", "--", "--", optimal, failed("ours"), time("ours")])
    return table(f"Other ways of choosing a program from the same e-graph as our extraction, on the {len(names)} "
                 "programs of Benchmark (2): extracting by size and then bounding, and greedy readout. Bounds are in "
                 "the error we report for each program. " + BSW.format(who="the method's bound", base="ours")
                 + ", and medians are of the method's bound divided by ours (lower is better). A method with no "
                 "program whose bound is finite counts as returning the input (Failed). Proved optimal counts the "
                 "programs where our extraction reaches a fixed point, so that no program in the e-graph has a smaller "
                 "bound. Time is the total over all programs, including building the e-graph (lower is better).",
                 "tab:extraction", "lcrrrr", "Method & b/s/w & Median & Proved optimal & Failed & Time (s)", rows)


def runtime():
    ours, bounds, rewrites, names, start = run.ours(), run.bounds(), run.rewrites(), run.benchmark2(), run.startup()

    def cells(reports, startup):
        xs = [max(0.0, r["seconds"] - startup) for r in reports if r.get("seconds")]
        return [*map(num, (median(xs), max(xs, default=None))), sum(r.get("status") != "ok" for r in reports)]

    bounding = {t: cells([bounds[n][t]["original"] for n in bounds], start[t]) for t in (*SOUND, "daisy")}
    bounding["ours"] = cells([{"seconds": o.get("bound_s"), "status": "ok" if said("abs", o["bound"]["abs"]) < INF
                               else "fail"} for o in ours.values()], start["egglog"])
    rewriting = {"herbie": cells([rewrites[n]["herbie"] for n in names], 0),
                 "daisy": cells([rewrites[n]["daisy"] for n in names], start["daisy-rewrite"]),
                 "ours": cells([{"seconds": seconds(n, "ours"), "status": "ok" if ours[n].get("programs") else "fail"}
                                for n in names], start["egglog"])}
    rows = [[NAMES[t], *bounding.get(t, ["--"] * 3), *rewriting.get(t, ["--"] * 3)]
            for t in (*SOUND, "herbie", "daisy", "ours")]
    head = ("& \\multicolumn{3}{c}{Bounding} & \\multicolumn{3}{c}{Rewriting} \\\\\n"
            "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
            "Tool" + " & Median & Max & Failed" * 2)
    return table(f"Running times in seconds (lower is better): bounding the {len(bounds)} programs of Benchmark (1), "
                 f"and rewriting the {len(names)} of Benchmark (2). Failed counts the runs that give no bound or no "
                 "program, including those stopped by a time limit; our extraction stopped by its limit still "
                 "returns the minimal-cost set of the last round it finished. Times leave out each tool's startup, "
                 "its least time on a trivial program.", "tab:runtime", "lrrrrrr", head, rows)


def conditions(text):
    """The conditions of the ifs in a program's text."""
    out, i = [], text.find("(if (")
    while i >= 0:
        out.append(text[i + 4:text.index(")", i + 4) + 1])
        i = text.find("(if (", i + 1)
    return out


def branching():
    results = run.branching()
    rows, disagree = [], 0
    for who in ("ours-plain", "ours", "herbie", "daisy"):
        branched = failed = 0
        errors, speeds = [], []
        for r in results.values():
            p = r.get(who) or {}
            disagree += "program" in p and not p["measured"]["equivalent"]
            if "program" not in p or not p["measured"]["equivalent"]:
                failed += 1
                p = None
            branched += bool(p and conditions(p["program"]))
            errors.append(((p or r["original"])["measured"]["ulps"], r["original"]["measured"]["ulps"]))
            base = (r.get("ours-plain") or {}).get("ns")
            if p and base and "ns" in p:
                speeds.append(p["ns"] / base)
        rows.append([NAMES[who], branched, compare(errors)[0], num(median([a for a, _ in errors])),
                     num(median(speeds)), failed])
    return table(f"Branching on the {len(results)} programs of Benchmark (3). Branched counts the programs where the "
                 "tool's program has an if-statement. Max sampled error is in ulps (lower is better): "
                 + BSW.format(who="it", base="the input's") + ", and its median is over all programs. Speed is the "
                 "median of run time divided by our unbranched program's (lower is better). A tool that returns no "
                 "program" + ", or one that disagrees with its input over the reals," * bool(disagree)
                 + " counts as returning its input (Failed).", "tab:branching", "lrcrrr",
                 "Tool & Branched & b/s/w & Median max error & Speed & Failed", rows)


def ablation():
    names, variants = run.benchmark2(), ("abs-only", "rel-only", "ours")
    tracks = {"abs": ("abs-only", "ours"), "rel": ("rel-only", "ours")}

    def cells(w, m):
        if w not in tracks[m]:
            return ["--"] * 3
        return [*compare((bound(n, w, m), raw(n, "original", m)) for n in names),
                sum(picked(n, w, m) == "original" for n in names)]

    rows = [[NAMES[w], *cells(w, "abs"), *cells(w, "rel"), num(round(sum(seconds(n, w) for n in names)))]
            for w in variants]
    head = ("& \\multicolumn{3}{c}{Absolute Bound} & \\multicolumn{3}{c}{Relative Bound} & \\\\\n"
            "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
            "Variant" + " & b/s/w & Median & Failed" * 2 + " & Time (s)")
    return table(f"Ablations on the {len(names)} programs of Benchmark (2), without branching, against the input. "
                 + BSW.format(who="the variant's bound", base="the input's")
                 + ", where a missing bound is infinite, and medians are of the variant's bound divided by the "
                 "input's (lower is better). A variant with no program whose bound is finite counts as returning the "
                 "input (Failed). Each error compares the variants that track it. Time is the total over all programs "
                 "(lower is better).", "tab:ablation", "lcrrcrrr", head, rows)


if __name__ == "__main__":
    broken = [n for n, r in run.originals().items() if "error" in r]
    if broken:
        print(f"warning: {len(broken)} inputs could not be measured, so Benchmark (2) leaves them out: {broken}")
    os.makedirs(os.path.join(run.RESULTS, "tables"), exist_ok=True)
    for f in (accuracy, tightness, extraction, runtime, branching, ablation):
        with open(os.path.join(run.RESULTS, "tables", f.__name__ + ".tex"), "w") as out:
            out.write(f())
