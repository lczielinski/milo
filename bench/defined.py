"""Whether a program is defined over the reals on all of its box: no division by zero and no square
root of a negative number. Outward-rounded interval arithmetic proves boxes defined, bisecting those it
cannot. An undefined program's witness is a corner, or a continuous divisor with opposite signs at two
corners."""

import itertools
import math
import random
from fractions import Fraction

import gmpy2
from gmpy2 import mpfr, mpz

import sample
from milo.fpcore import show

DOWN = gmpy2.context(precision=256, round=gmpy2.RoundDown)
UP = gmpy2.context(precision=256, round=gmpy2.RoundUp)
INF = mpfr("inf")
WHOLE = (-INF, INF)
BOXES, CORNERS = 5000, 64


class Undefined(Exception):
    pass


def exact(v):
    return Fraction(*v.as_integer_ratio())


def constant(q):
    n, d = mpz(q.numerator), mpz(q.denominator)
    return DOWN.div(n, d), UP.div(n, d)


def finite(*intervals):
    return all(gmpy2.is_finite(lo) and gmpy2.is_finite(hi) for lo, hi in intervals)


def combine(down, up, a, b):
    if not finite(a, b):
        return WHOLE
    return min(down(x, y) for x in a for y in b), max(up(x, y) for x in a for y in b)


def branches(e):
    return e[0] == "If" or e[0] not in ("Var", "Num") and any(map(branches, e[1:]))


def restrict(box, op, x, k, taken):
    """The box cut down to its points that take (or skip) the condition x op k, or None if none do."""
    above, strict = (op[0] == ">") == taken, (op in (">", "<")) == taken
    lo, hi = box[x]
    klo, khi = constant(k)
    if above:
        if exact(hi) < k or strict and exact(hi) == k:
            return None
        return dict(box, **{x: (max(lo, klo), hi)})
    if exact(lo) > k or strict and exact(lo) == k:
        return None
    return dict(box, **{x: (lo, min(hi, khi))})


def walk(e, box, issues, seen=None):
    """An interval holding e's values on the box. Adds each division and square root that may be
    undefined to issues, as (kind, node, continuous), and raises Undefined for one that is undefined on
    all of the box. seen, if given, records each divisor's interval."""
    match e:
        case ("Var", x):
            return box[x]
        case ("Num", q):
            return constant(q)
        case ("Neg", a):
            lo, hi = walk(a, box, issues, seen)
            return -hi, -lo
        case ("Sqrt", a):
            lo, hi = walk(a, box, issues, seen)
            if hi < 0:
                raise Undefined(f"square root of a negative number at {point(box)}")
            if lo < 0:
                issues.append(("sqrt", e, False))
            return (DOWN.sqrt(lo) if lo > 0 else mpfr(0)), UP.sqrt(hi)
        case ("If", (op, x, k), a, b):
            parts = [walk(c, cut, issues, seen)
                     for c, t in ((a, True), (b, False)) if (cut := restrict(box, op, x, k, t))]
            return min(p[0] for p in parts), max(p[1] for p in parts)
        case ("Div", a, b):
            num, inner = walk(a, box, issues, seen), []
            den = walk(b, box, inner, seen)
            issues += inner
            if seen is not None:
                seen[id(e)] = den
            if den[0] == 0 == den[1]:
                raise Undefined(f"divides {show(b)} = 0 at {point(box)}")
            if den[0] <= 0 <= den[1]:
                issues.append(("div", e, not inner and not branches(b)))
                return WHOLE
            return combine(DOWN.div, UP.div, num, den)
        case ("Add", a, b):
            (alo, ahi), (blo, bhi) = walk(a, box, issues, seen), walk(b, box, issues, seen)
            return (DOWN.add(alo, blo), UP.add(ahi, bhi)) if finite((alo, ahi), (blo, bhi)) else WHOLE
        case ("Mul", a, b):
            return combine(DOWN.mul, UP.mul, walk(a, box, issues, seen), walk(b, box, issues, seen))
    raise ValueError(f"unknown expression {e}")


def point(box):
    return ", ".join(f"{x}={float(lo):.17g}" if lo == hi else f"{x} in [{float(lo):.17g}, {float(hi):.17g}]"
                     for x, (lo, hi) in box.items())


def corners(box):
    every = itertools.product(*[sorted({lo, hi}) for lo, hi in box.values()])
    found = [dict(zip(box, c)) for c in itertools.islice(every, 1 << 16)]
    if len(found) > CORNERS:
        found = random.Random(0).sample(found, CORNERS)
    return found


def witness(e, box, issues):
    """Raises Undefined if a corner of the box is undefined, or if a continuous divisor of an issue
    has opposite signs at two corners that reach it (a corner that takes another branch never does)."""
    values = []
    for c in corners(box):
        seen = {}
        walk(e, {x: (v, v) for x, v in c.items()}, [], seen)
        values.append(seen)
    for kind, node, continuous in issues:
        if kind == "div" and continuous:
            signs = {(den[1] < 0) - (den[0] > 0) for seen in values if (den := seen.get(id(node)))}
            if {-1, 1} <= signs:
                raise Undefined(f"divisor {show(node[2])} changes sign on {point(box)}")


def size(lo, hi):
    """How many floats an interval spans, or, below 2, what fraction of an ulp it spans."""
    if lo == hi:
        return 0
    if float(lo) == lo and float(hi) == hi:
        span = sample.ordinal(float(hi), "binary64") - sample.ordinal(float(lo), "binary64")
        if span >= 2:
            return span
    return float(hi - lo) / math.ulp(float(lo))


def bisect(box):
    """The box's two halves along its widest variable, split at the middle float or, between two
    adjacent floats, at the middle real; None if no variable can be split."""
    sizes = {x: size(lo, hi) for x, (lo, hi) in box.items()}
    x = max(sizes, key=sizes.get)
    lo, hi = box[x]
    if sizes[x] >= 2:
        mid = mpfr(sample.from_ordinal(sample.ordinal(float(lo), "binary64") + sizes[x] // 2, "binary64"))
    else:
        mid = DOWN.div(DOWN.add(lo, hi), 2)
    if not lo < mid < hi:
        return None
    return [dict(box, **{x: (lo, mid)}), dict(box, **{x: (mid, hi)})]


def check(e, box):
    """'defined', 'undefined' or 'unknown' (out of boxes, or a point that is neither), with a witness
    for 'undefined'."""
    work, status = [{x: (mpfr(lo), mpfr(hi)) for x, (lo, hi) in box.items()}], "defined"
    for _ in range(BOXES):
        if not work:
            return {"status": status}
        b = work.pop()
        issues = []
        try:
            walk(e, b, issues)
            if issues:
                witness(e, b, issues)
        except Undefined as u:
            return {"status": "undefined", "witness": str(u)}
        if issues:
            halves = bisect(b)
            if halves:
                work += halves
            else:
                status = "unknown"
    return {"status": "unknown" if work else status}
