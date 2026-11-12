import argparse

from .analysis import BOUNDS, FAIL, set_format
from .egraph import build, guarded
from .extract import analyze, best, extract
from .fpcore import parse, show


def main(argv=None):
    parser = argparse.ArgumentParser(prog="milo", description="Extract the programs with the "
                                     "smallest absolute and relative error bounds.")
    parser.add_argument("file", help="an FPCore file")
    parser.add_argument("--iters", type=int, default=4, help="rewriting iterations; 0 bounds the input")
    parser.add_argument("--no-branch", action="store_true", help="do not add if-statements")
    parser.add_argument("--time-limit", type=float, default=120,
                        help="seconds for extraction, after which F_k of the last finished round is reported")
    parser.add_argument("--bound-only", action="store_true",
                        help="build the e-graph as usual, then bound only the input program")
    args = parser.parse_args(argv)
    with open(args.file) as f:
        core = parse(f.read())
    set_format(core.precision)
    g = build(core, args.iters, not args.no_branch)
    I = g.interval[g.root]
    if args.bound_only:
        e = analyze(g, guarded(core.body, core.box, core.precision))
        print(f"{core.name}: {sum(map(len, g.nodes.values()))} e-nodes, {len(g.nodes)} e-classes, range in {I}")
        for name, bound in BOUNDS.items():
            print(f"{name} none: the input may fail" if e is FAIL else f"{name} {float(bound(e)):.3e}")
        return
    F, k, cut = extract(g, args.time_limit)
    print(f"{core.name}: {sum(map(len, g.nodes.values()))} e-nodes, {len(g.nodes)} e-classes, "
          f"range in {I}, minimal-cost set of {len(F[g.root])}")
    if cut:
        print(f"no fixed point: the time limit stopped extraction, so these are optimal only up to height {k}")
    else:
        print(f"fixed point at height {k}: these are optimal among all programs in the e-graph")
    for name, bound in BOUNDS.items():
        if F[g.root]:
            e, p = best(F[g.root], bound)
            print(f"{name} {float(bound(e)):.3e}  {show(p)}")
        else:
            print(f"{name} none: " + (f"no program up to height {k} has a bound" if cut else "every program may fail"))
