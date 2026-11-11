"""Measured error: programs against an exact reference, on the same points of their box.

The points are uniform over the floats in the box, plus the points a hill climb finds
from each program's worst ones, so every program is also scored where its rivals do
worst. Each program is checked against the original over the reals at the first points.
"""

import math
import random
import statistics
import struct
from fractions import Fraction

import gmpy2
from gmpy2 import mpfr

POINTS = 10_000
CLIMBS = 8
STEPS = 4000
CHECKS = 20
FORMATS = {"binary64": ("<d", "<Q", 1 << 63), "binary32": ("<f", "<I", 1 << 31)}


def single(x):
    try:
        return struct.unpack("f", struct.pack("f", x))[0]
    except OverflowError:
        return math.copysign(math.inf, x)


def rounded(q, fmt):
    """A Fraction or mpfr rounded once to the format."""
    ctx = gmpy2.ieee(32 if fmt == "binary32" else 64)
    return float(ctx.div(q.numerator, q.denominator) if isinstance(q, Fraction) else ctx.plus(q))


def ordinal(x, fmt):
    f, i, sign = FORMATS[fmt]
    n = struct.unpack(i, struct.pack(f, x))[0]
    return -(n & ~sign) if n & sign else n


def from_ordinal(n, fmt):
    f, i, sign = FORMATS[fmt]
    return struct.unpack(f, struct.pack(i, (-n | sign) if n < 0 else n))[0]


def taken(op, x, k):
    return x > k or (op == ">=" and x == k) if op[0] == ">" else x < k or (op == "<=" and x == k)


def run(e, env, fmt):
    """The program in floating point, with IEEE infinities and NaNs."""
    rnd = single if fmt == "binary32" else float
    match e:
        case ("Var", x):
            return env[x]
        case ("Num", q):
            return rounded(q, fmt)
        case ("If", (op, x, k), a, b):
            return run(a if taken(op, env[x], rounded(k, fmt)) else b, env, fmt)
        case ("Neg", a):
            return -run(a, env, fmt)
        case ("Sqrt", a):
            a = run(a, env, fmt)
            return rnd(math.sqrt(a)) if a >= 0 else math.nan
    a, b = run(e[1], env, fmt), run(e[2], env, fmt)
    match e[0]:
        case "Add":
            return rnd(a + b)
        case "Mul":
            return rnd(a * b)
        case "Div" if b == 0:
            return math.nan if a == 0 or a != a else math.copysign(math.inf, a) * math.copysign(1, b)
        case "Div":
            return rnd(a / b)


def real(e, env, peak):
    """The program over the reals at the current mpfr precision; peak[0] becomes the largest
    magnitude of any value along the way."""
    match e:
        case ("Var", x):
            v = mpfr(env[x])
        case ("Num", q):
            v = mpfr(q.numerator) / q.denominator
        case ("If", (op, x, k), a, b):
            v = real(a if taken(op, env[x], k) else b, env, peak)
        case ("Neg", a):
            v = -real(a, env, peak)
        case ("Sqrt", a):
            v = gmpy2.sqrt(real(a, env, peak))
        case (op, a, b):
            a, b = real(a, env, peak), real(b, env, peak)
            v = a + b if op == "Add" else a * b if op == "Mul" else a / b
    peak[0] = max(peak[0], abs(v))
    return v


def exact(e, env):
    """The value over the reals; None if it is undefined (some value along the way is not finite,
    since mpfr takes x/0 to an infinity and then 1/inf to 0) or never settles. The precision
    doubles until two results agree and stand clear of the rounding noise of the largest value
    along the way; a result still inside that noise at 16384 bits is 0."""
    before, clear = None, False
    with gmpy2.context(precision=128) as ctx:
        while ctx.precision <= 16384:
            peak = [mpfr(0)]
            now = real(e, env, peak)
            if not (gmpy2.is_finite(now) and gmpy2.is_finite(peak[0])):
                return None
            clear = abs(now) > peak[0] * gmpy2.exp2(64 - ctx.precision)
            if clear and before is not None and abs(now - before) <= abs(now) * 2.0 ** -64:
                return now
            before = now
            ctx.precision *= 2
    return None if clear else mpfr(0)


def score(e, env, ref, fmt):
    """The error in ulps at one point: how many floats lie between the program's result and the
    exact value rounded. An infinite or NaN result is the worst possible."""
    got = run(e, env, fmt)
    if not math.isfinite(got):
        return 1 << (32 if fmt == "binary32" else 64)
    return abs(ordinal(got, fmt) - ordinal(rounded(ref, fmt), fmt))


def climb(f, start, box, fmt, budget):
    """Pattern search for worse points, in steps of float ordinals that halve down to one."""
    best, now, level, found = f(start), start, 1, []
    while budget > 0 and any((b - a) >> level for a, b in box.values()):
        moved = False
        for x, (a, b) in box.items():
            for step in (-1, 1):
                o = ordinal(now[x], fmt) + step * max(1, (b - a) >> level)
                if budget > 0 and a <= o <= b:
                    budget -= 1
                    env = dict(now, **{x: from_ordinal(o, fmt)})
                    got = f(env)
                    if got is not None and (best is None or got > best):
                        best, now, moved = got, env, True
                        found.append(env)
        level += not moved
    return found


def spans(core):
    """Each variable's floats in the box, as a range of ordinals."""
    out, fmt = {}, core.precision
    for x, (lo, hi) in core.box.items():
        a, b = ordinal(lo, fmt), ordinal(hi, fmt)
        out[x] = (a + (from_ordinal(a, fmt) < lo), b - (from_ordinal(b, fmt) > hi))
    return out


def points(core, n=POINTS):
    rng = random.Random(core.name)
    return [{x: from_ordinal(rng.randint(a, b), core.precision) for x, (a, b) in spans(core).items()}
            for _ in range(n)]


def measure(core, original, programs):
    """For each named program: its worst and mean error in ulps, and whether it agreed with the
    original over the reals."""
    fmt, box, refs = core.precision, spans(core), {}

    def ref(env):
        key = tuple(env.values())
        if key not in refs:
            refs[key] = exact(original, env)
        return refs[key]

    uniform = [env for env in points(core) if ref(env) is not None]
    found = []
    for e in programs.values():
        def ulps(env):
            return None if ref(env) is None else score(e, env, ref(env), fmt)
        for env in sorted(uniform, key=ulps)[-CLIMBS:]:
            found += climb(ulps, env, box, fmt, STEPS // CLIMBS)
    pool = uniform + [env for env in found if ref(env) is not None]
    out = {}
    for name, e in programs.items():
        scores = [score(e, env, ref(env), fmt) for env in pool]
        agrees = all((v := exact(e, env)) is not None and abs(v - ref(env)) <= 1e-12 * max(abs(v), abs(ref(env)))
                     for env in uniform[:CHECKS])
        out[name] = {"ulps": max(scores, default=None),
                     "mean": statistics.fmean(scores[:len(uniform)]) if uniform else None, "equivalent": agrees}
    return out
