"""The numbers the evaluation's text gives, and the data of its figures, from results/*.json:
`uv run python bench/tables.py [accuracy ...]`. A tool that fails counts as returning its input, and a
relative bound of 1 or more counts as no bound."""

import math
import statistics
import sys

import run
from run import SIZES

INF = math.inf
NAMES = {"fptaylor": "FPTaylor", "gappa": "Gappa", "precisa": "PRECiSA", "daisy": "Daisy", "herbie": "Herbie",
         "ours": "Milo", "bound-only": "Milo", "ours-plain": "Milo, no branch",
         "abs-only": "Absolute error only", "rel-only": "Relative error only"}
SOUND = ("fptaylor", "gappa", "precisa")
BELOW, LOW, HIGH = 1e-3, 1e-5, 1e20
WINDOW = {"abs": (1e-20, 1e20), "rel": (1e-19, 1)}


def num(x):
    if x is None:
        return "--"
    return str(int(x)) if float(x).is_integer() and abs(x) < 1e6 else f"{x:.3g}"


def percent(ratio):
    """A ratio as the percent change it is, such as 8% for 1.08."""
    if ratio is None:
        return "--"
    return f"{round((ratio - 1) * 100)}%"


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


def gain_loss(pairs):
    """For pairs (value, baseline), where lower is better: how many values are better than their baseline
    by change(), and how many worse."""
    found = [change(a, b) for a, b in pairs]
    return [found.count(0), found.count(2)]


def error_gain_loss(pairs):
    """For pairs (value, baseline) of average errors in bits: how many values are at least one bit
    better, and how many at least one bit worse."""
    pairs = list(pairs)
    return [sum(a <= b - 1 for a, b in pairs), sum(a >= b + 1 for a, b in pairs)]


def table(head, rows):
    """Rows of cells as aligned columns under a header row."""
    rows = [head, *[[str(c) for c in r] for r in rows]]
    widths = [max(map(len, column)) for column in zip(*rows)]
    return "\n".join("  ".join(c.ljust(w) for c, w in zip(r, widths)).rstrip() for r in rows) + "\n"


def metric(name):
    """run.metric() of our program on a Rewriting program, or of the input when we found none."""
    o = run.ours()[name]
    return run.metric(o["bounds"] if o["programs"] else o["original"])


def sound(name, who):
    """The smallest of FPTaylor's, Gappa's and PRECiSA's bounds on who's program, in the program's metric()."""
    m = metric(name)
    return min(said(m, run.bounds()[name][t].get(who, {}).get(m)) for t in SOUND)


def got(t, n):
    """t's program on a Rewriting program, or the input if t returned none or one that is not equivalent
    to it over the reals: one that disagrees at the sampled checks or is undefined somewhere in the box."""
    same = run.measured()[n].get(t, {}).get("equivalent") and run.definedness()[n][t]["status"] != "undefined"
    return t if same else "original"


def accuracy():
    """Rewriting against the input, the numbers Q1's text gives (the paper plots versus() instead)."""
    measured, speeds = run.measured(), run.speeds()
    names = run.rewriting()
    tools = ("herbie", "daisy", "ours")

    def speed(t):
        return geomean([speeds[n][got(t, n)] / speeds[n]["original"] for n in names
                        if speeds[n].get(got(t, n)) and speeds[n].get("original")])

    def bounds(t):
        return [(sound(n, got(t, n)), sound(n, "original")) for n in names]

    def errors(t):
        return [(measured[n][got(t, n)]["bits"], measured[n]["original"]["bits"]) for n in names]

    rows = [[NAMES[t], *gain_loss(bounds(t)), *error_gain_loss(errors(t)), percent(speed(t)),
             median_time(rewriting_seconds(t, n) for n in names), sum(got(t, n) == "original" for n in names)]
            for t in tools]
    return table(["Tool", "Bound gain", "Bound loss", "Error gain", "Error loss", "Slowdown", "Time (s)", "Failed"],
                 rows)


def improvement(base, value):
    """How much a bound improves on base (base over value), as a pgfplots coordinate. Log axes have no 0
    or infinity, so an improvement below BELOW, or a lost bound, is drawn at LOW, and an infinite one (only
    the rewrite has a bound) at HIGH."""
    q = quotient(base, value)
    return f"{HIGH if q == INF else LOW if q < math.log(BELOW) else math.exp(q):.3g}"


def windowed(m, v):
    """An m-bound as a pgfplots coordinate: one below the error's WINDOW, or 0, is drawn a thousandfold
    below it, and one above it, or none, a thousandfold above it."""
    lo, hi = WINDOW[m]
    return f"{lo / 1e3 if v < lo else hi * 1e3 if v > hi else v:.3g}"


def versus():
    """Our rewrite's improvement on each Rewriting input against Herbie's and Daisy's, as pgfplots rows:
    in sound bound by improvement(), and in average sampled error as the bits it drops."""
    measured, out = run.measured(), []
    for t in ("herbie", "daisy"):
        bounds, errors = [f"% {t} bound", "x y"], [f"% {t} error", "x y"]
        for n in run.rewriting():
            base, bits = sound(n, "original"), measured[n]["original"]["bits"]
            bounds.append(" ".join(improvement(base, sound(n, got(w, n))) for w in ("ours", t)))
            errors.append(" ".join(f"{bits - measured[n][got(w, n)]['bits']:.3g}" for w in ("ours", t)))
        out += bounds + errors
    return "\n".join(out) + "\n"


def analyzers():
    """Our bound on each input (milo --bound-only) against each other tool's, per error, as pgfplots rows."""
    only, bounds, out = run.bound_only(), run.bounds(), []
    for m in ("abs", "rel"):
        for t in (*SOUND, "daisy"):
            out += [f"% {t} {m}", "x y"]
            for n in only:
                xs = (only[n].get(m), bounds[n][t].get("original", {}).get(m))
                out.append(" ".join(windowed(m, said(m, v)) for v in xs))
    return "\n".join(out) + "\n"


def median_time(xs):
    """The median of running times, leaving out runs with no time."""
    return num(median([x for x in xs if x is not None]))


def rewriting_seconds(t, n):
    """t's running time on a Rewriting program, less its startup (Herbie's own times have none)."""
    start = run.startup()
    if t == "ours":
        return max(0.0, seconds(n, "ours") - start["egglog"])
    r = run.rewrites()[n][t]
    return max(0.0, r["seconds"] - (start["daisy-rewrite"] if t == "daisy" else 0)) if r.get("seconds") else None


def tightness():
    """How many inputs each tool bounds, in each error, and its median time: the numbers Q2's text gives
    (the paper plots analyzers() instead). Our bounds are milo --bound-only's."""
    only, bounds, start = run.bound_only(), run.bounds(), run.startup()

    def bound(t, m, n):
        return said(m, (only[n] if t == "bound-only" else bounds[n][t].get("original", {})).get(m))

    def seconds(t, n):
        x = (only[n] if t == "bound-only" else bounds[n][t].get("original", {})).get("seconds")
        return None if x is None else max(0.0, x - start["egglog" if t == "bound-only" else t])

    rows = [[NAMES[t], *(sum(bound(t, m, n) < INF for n in only) for m in ("abs", "rel")),
             median_time(seconds(t, n) for n in only)] for t in ("bound-only", *SOUND, "daisy")]
    return table(["Tool", "Absolute", "Relative", "Time (s)"], rows)


def raw(name, who, m):
    """The m-bound of who's program on a Rewriting program, with our analysis in our e-graph."""
    o = run.ours()[name]
    if who == "ours":
        x = o["bounds"][m] if o.get("programs") else None
    elif who == "original":
        x = o["original"][m]
    else:
        x = run.alternatives()[name].get(who, {}).get(m, {}).get(m)
    return said(m, x)


def picked(name, who, m):
    """who, or the input if who has no program with a finite m-bound."""
    return who if raw(name, who, m) < INF else "original"


def bound(name, who, m):
    return raw(name, picked(name, who, m), m)


def seconds(name, who):
    """who's running time on a Rewriting program, building the e-graph included."""
    if who == "ours":
        return run.ours()[name]["egglog_s"] + run.ours()[name]["extract_s"]
    if who.endswith("-only"):
        return run.single_bounds()[name][who[:3]]["seconds"]
    return run.alternatives()[name][who]["seconds"]


def own(n, m):
    """The m-bound of the m-only analysis's program on a Rewriting program, by that analysis."""
    return said(m, run.single_bounds()[n][m]["own"])


def extraction():
    """Gain and loss of each method's bound against ours, in each error: the numbers Q3's text gives.
    Extraction by size is bounded by our analysis and counts as returning the input when no bound is
    finite; the single-metric analyses are bounded by their own."""
    names = run.rewriting()
    sizes = [f"size{k}" for k in SIZES]
    tracks = {"abs": sizes + ["abs-only"], "rel": sizes + ["rel-only"]}

    def cells(w, m):
        if w not in tracks[m]:
            return ["--"] * 2
        if w.endswith("-only"):
            return gain_loss((own(n, m), bound(n, "ours", m)) for n in names)
        return gain_loss((bound(n, w, m), bound(n, "ours", m)) for n in names)

    optimal = sum(not run.ours()[n].get("cut", True) for n in names)
    rows = [[NAMES.get(w, f"Size, top-{w[4:]}"), *cells(w, "abs"), *cells(w, "rel"),
             optimal if w == "ours" else "--", median_time(seconds(n, w) for n in names)]
            for w in sizes + ["abs-only", "rel-only", "ours"]]
    return table(["Method", "Abs gain", "Abs loss", "Rel gain", "Rel loss", "Optimal", "Time (s)"], rows)


def branch_got(t, n):
    """t's program on a Branching program, or the input if t returned none or one not equivalent to it."""
    p = run.branching()[n].get(t) or {}
    same = "program" in p and p["measured"]["equivalent"] and p.get("defined", {}).get("status") != "undefined"
    return t if same else "original"


def branch_sound(n, who):
    """The largest over the pieces of who's program on a Branching program of the smallest of FPTaylor's,
    Gappa's and PRECiSA's bounds, in relative error where our branched program's is below 1."""
    o = run.branching()[n].get("ours") or {}
    m = run.metric(o.get("bounds") or {"rel": run.branching()[n]["original"]["bound"].get("rel", INF)})
    reports = run.branching_bounds()[n].get(who)
    if not reports:
        return INF
    return max(min(said(m, r.get(m)) for r in piece) for piece in zip(*(reports[t] for t in SOUND)))


def branching():
    """Branching as accuracy(), with how many programs branch: the numbers Q4's text gives."""
    results, start = run.branching(), run.startup()
    names = list(results)
    got, sound = branch_got, branch_sound

    def seconds(t, n):
        x = (results[n].get(t) or {}).get("seconds")
        s0 = start["egglog"] if t in ("ours", "ours-plain") else start["daisy-rewrite"] if t == "daisy" else 0
        return None if x is None else max(0.0, x - s0)

    def speed(t):
        return geomean([results[n][got(t, n)]["ns"] / results[n]["original"]["ns"] for n in names
                        if results[n][got(t, n)].get("ns") and results[n]["original"].get("ns")])

    rows = [[NAMES[t], sum("(if " in results[n][got(t, n)]["program"] for n in names),
             *gain_loss((sound(n, got(t, n)), sound(n, "original")) for n in names),
             *error_gain_loss((results[n][got(t, n)]["measured"]["bits"], results[n]["original"]["measured"]["bits"])
                              for n in names),
             percent(speed(t)), median_time(seconds(t, n) for n in names), sum(got(t, n) == "original" for n in names)]
            for t in ("herbie", "daisy", "ours-plain", "ours")]
    return table(["Tool", "Branched", "Bound gain", "Bound loss", "Error gain", "Error loss", "Slowdown", "Time (s)",
                  "Failed"], rows)


if __name__ == "__main__":
    every = (accuracy, tightness, extraction, branching, versus, analyzers)
    for f in [f for f in every if f.__name__ in sys.argv[1:]] or every:
        print(f"% {f.__name__}\n{f()}")
