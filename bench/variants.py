"""The alternatives to our extraction of the minimal-cost set that the evaluation compares it with."""

import heapq
import time

from milo import analysis as A
from milo.analysis import BOUNDS, FAIL, ONE, TOP, ZERO, err, mix, power_of_two, sterbenz
from milo.extract import extract, participants, program
from milo.interval import Iv


def analyze(g, p):
    """err of the program p as a term of g, with its conditions written as guarded writes them."""
    index = {(n.op, n.payload, n.children): c for c, ns in g.nodes.items() for n in ns}

    def node(op, payload, kids):
        c = index[(op, payload, tuple(k for k, _ in kids))]
        Xs = [g.interval[k] for k, _ in kids]
        return c, err(op, g.interval[c], Xs, [e for _, e in kids])

    def walk(p, guard):
        match p:
            case ("Var", x) if x in guard:
                return node("Ctx", guard[x], [walk(p, {})])
            case ("Var" | "Num", v):
                return node(p[0] if (p[0], v, ()) in index else "Lit", v, [])
            case ("If", (_, x, k), a, b):
                return node("If", (x, k), [walk(a, guard | {x: (x, k, True)}), walk(b, guard | {x: (x, k, False)})])
            case (op, *args):
                return node(op, None, [walk(a, guard) for a in args])

    return walk(p, {})[1]


def smallest(g, n, time_limit=None):
    """The n smallest programs of the root class, counting operations and leaves, found in
    rounds like extraction's."""
    deadline = time_limit and time.monotonic() + time_limit
    nodes, members, readers = participants(g)
    best, shares, todo = {c: [] for c in g.nodes}, {}, range(len(nodes))
    while todo and not (deadline and time.monotonic() > deadline):
        for i in todo:
            shares[i] = sums(nodes[i][1], [best[x] for x in nodes[i][1].children], n)
        changed = {}
        for c in {nodes[i][0] for i in todo}:
            found = heapq.nsmallest(n, (e for i in members[c] for e in shares.get(i, ())))
            if found != best[c]:
                changed[c] = found
        best.update(changed)
        todo = sorted({i for c in changed for i in readers[c]})
    return [p for _, p in best[g.root]]


def sums(m, lists, n):
    """The n smallest terms with root m, from its children's sorted (size, program) lists."""
    if not lists:
        return [(1, program(m, []))]
    if len(lists) == 1:
        return [(1 + s, program(m, [p])) for s, p in lists[0][:n]]
    a, b = lists
    if not a or not b:
        return []
    out, seen, heap = [], {(0, 0)}, [(a[0][0] + b[0][0], 0, 0)]
    while heap and len(out) < n:
        s, i, j = heapq.heappop(heap)
        out.append((1 + s, program(m, [a[i][1], b[j][1]])))
        for i2, j2 in ((i + 1, j), (i, j + 1)):
            if i2 < len(a) and j2 < len(b) and (i2, j2) not in seen:
                seen.add((i2, j2))
                heapq.heappush(heap, (a[i2][0] + b[j2][0], i2, j2))
    return out


def greedy(metric):
    """The order that keeps a single cost per class: the one with the smallest bound."""
    return lambda e1, e2: BOUNDS[metric](e1) <= BOUNDS[metric](e2)


def abs_only(op, I, Xs, Ds):
    """x~ - x in D, tracked alone: no S, no reduction."""
    if any(D is FAIL for D in Ds):
        return FAIL

    def rounded(D, J=TOP, exact=False, scale=None, eta=A.ETA):
        Ih = (I + D) & J
        if Ih.mag > A.OMEGA:
            return FAIL
        if exact or scale is not None and (scale >= 0 or Ih.mig >= A.TINY):
            return D
        h = A.gamma_a(Ih, eta)
        return D + Iv(-h, h)

    def k(i):
        return power_of_two(Xs[i], A.Pair(ONE, Ds[i]))

    match op:
        case "Var" | "Num":
            return ZERO
        case "Lit":
            return rounded(ZERO)
        case "Neg":
            return -Ds[0]
        case "Add":
            Xt, Yt = Xs[0] + Ds[0], Xs[1] + Ds[1]
            return rounded(Ds[0] + Ds[1], Xt + Yt, exact=sterbenz(Xt, Yt), eta=0)
        case "Mul":
            ks = [x for x in (k(0), k(1)) if x is not None]
            D = Xs[0] * Ds[1] + Xs[1] * Ds[0] + Ds[0] * Ds[1]
            return rounded(D, (Xs[0] + Ds[0]) * (Xs[1] + Ds[1]), scale=max(ks, default=None))
        case "Div":
            Yt = Xs[1] + Ds[1]
            if 0 in Yt:
                return FAIL
            scale = None if k(1) is None else -k(1)
            return rounded((Ds[0] - I * Ds[1]) / Yt, (Xs[0] + Ds[0]) / Yt, scale=scale)
        case "Sqrt":
            Xt = Xs[0] + Ds[0]
            if Xt.lo < 0:
                return FAIL
            return rounded(Ds[0] / (Xt.sqrt() + Xs[0].sqrt()), Xt.sqrt())


def rel_only(op, I, Xs, Ss):
    """x~/x in S, tracked alone: no D, no reduction, and nothing where a value can be 0."""
    if any(S is FAIL for S in Ss):
        return FAIL

    def rounded(S, J=TOP, exact=False, scale=None, eta=A.ETA):
        Ih = (TOP if 0 in I else I * S) & J
        if Ih.mag > A.OMEGA:
            return FAIL
        if exact or scale is not None and (scale >= 0 or Ih.mig >= A.TINY):
            return S
        g = A.gamma_r(Ih, eta)
        return S * (ONE + Iv(-g, g))

    def k(i):
        return power_of_two(Xs[i], A.Pair(Ss[i], ZERO))

    def computed(i):
        return TOP if 0 in Xs[i] else Xs[i] * Ss[i]

    match op:
        case "Var" | "Num":
            return ONE
        case "Lit":
            return rounded(ONE)
        case "Neg":
            return Ss[0]
        case "Add":
            S = TOP
            if 0 not in Xs[0] and 0 not in Xs[1]:
                X, Y = Xs
                lam = (X / I) & (ONE - Y / I) & (ONE / (ONE + Y / X))
                S = mix(lam.lo, *Ss) | mix(lam.hi, *Ss)
            Xt, Yt = computed(0), computed(1)
            return rounded(S, Xt + Yt, exact=sterbenz(Xt, Yt), eta=0)
        case "Mul":
            ks = [x for x in (k(0), k(1)) if x is not None]
            return rounded(Ss[0] * Ss[1], computed(0) * computed(1), scale=max(ks, default=None))
        case "Div":
            if 0 in computed(1):
                return FAIL
            scale = None if k(1) is None else -k(1)
            return rounded(Ss[0] / Ss[1], computed(0) / computed(1), scale=scale)
        case "Sqrt":
            if computed(0).lo < 0:
                return FAIL
            return rounded(Ss[0].sqrt(), computed(0).sqrt())


def single(g, metric, time_limit=None):
    """The program extracted by tracking one metric alone."""
    cost, bound = (abs_only, lambda D: D.mag) if metric == "abs" else (rel_only, lambda S: (S - ONE).mag)
    F, _, _ = extract(g, time_limit, cost=cost, leq=lambda a, b: a <= b)
    return min(F[g.root], key=lambda entry: bound(entry[0]))[1] if F[g.root] else None
