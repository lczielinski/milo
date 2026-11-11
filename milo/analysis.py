"""The error analysis. For a term with real value x and computed value x~,
err returns FAIL or a pair (S, D) with x~/x in S (when x != 0) and x~ - x in D.

Each operation follows the same steps: compute (S^, D^), the pair for the exact result of the
operation on the computed operands; reduce it; find I^, where that pre-rounding value lies; FAIL if
I^ may overflow; perform one rounding unless exact."""

from collections import namedtuple

import gmpy2

from .interval import ONE, TOP, UP, ZERO, Iv

Pair = namedtuple("Pair", "S D")
FAIL = None
EXACT = Pair(ONE, ZERO)
FORMATS = {"binary64": (53, -1022), "binary32": (24, -126)}


def set_format(name):
    """Set the rounding model for precision p and minimum exponent emin: U = 2^-p,
    ETA = eta = 2^(emin - p), TINY = 2^emin, the smallest normal, OMEGA the largest float,
    and FLOATS, a context that rounds to the format."""
    global U, ETA, TINY, OMEGA, FLOATS
    p, emin = FORMATS[name]
    FLOATS = gmpy2.ieee(64 if name == "binary64" else 32)
    U, ETA, TINY = gmpy2.exp2(-p), gmpy2.exp2(emin - p), gmpy2.exp2(emin)
    OMEGA = (2 - 2 * U) * gmpy2.exp2(1 - emin)


set_format("binary64")


def precedes(e1, e2):
    """The order on costs (S, D) <= (S', D')."""
    return e1.S <= e2.S and e1.D <= e2.D


def computed(X, e):
    """X~ = (X + D) & XS, where the computed value lies, for a range X and a pair e.
    XS needs x != 0, so it is dropped when 0 is in X."""
    return X + e.D if 0 in X else (X + e.D) & (X * e.S)


def reduce(I, S, D):
    """The reduction: (S & (1 + D/I), D & I(S - 1)). The second intersection needs
    x != 0, so it is dropped when 0 is in I."""
    return Pair(S & (ONE + D / I), D if 0 in I else D & (I * (S - ONE)))


# Rounding. Every real r with |r| <= Omega has |fl(r) - r| <= e(r), where e(r) = 0 for a float and
# max(u * ufp(r), eta) otherwise, or u * ufp(r) for a sum. gamma_a and gamma_r take its worst case over I^.

def ufp(r):
    """2^floor(log2 |r|), the value of r's leading bit."""
    return gmpy2.exp2(gmpy2.frexp(r)[0] - 1)


def ufp_below(M):
    """The largest ufp(r) of a non-float r with |r| <= M."""
    return ufp(M) / 2 if ufp(M) == M else ufp(M)


def gamma_a(Ih, eta):
    """sup of e(r) over I^."""
    return max(U * ufp_below(Ih.mag), eta) if Ih.mag else 0


def gamma_r(Ih, eta):
    """min(1, sup of e(r)/|r| over the nonzero r in I^)."""
    m = Ih.mig
    if m == 0:
        return 1 if eta else U
    g = UP.div(U * ufp(m), m) if ufp(m) == ufp_below(Ih.mag) else U
    return min(1, max(g, UP.div(eta, m)))


def rounded(I, S, D, J=TOP, exact=False, scale=None, is_sum=False):
    """The pair of an operation, from its pre-rounding pair (S^, D^) and J = f(X~, Y~), the
    operation on the computed operands. exact is Sterbenz's lemma, and scale is k when the
    operation multiplies by 2^k. A sum has no eta, since a sum of floats in the subnormals is exact."""
    e = reduce(I, S, D)
    Ih = computed(I, e) & J
    if Ih.mag > OMEGA:
        return FAIL
    if scale is not None:
        exact = exact or scale >= 0 or Ih.mig >= TINY  # scaling down may round in the subnormals
    if exact or Ih.lo == Ih.hi == FLOATS.plus(Ih.lo):  # a float rounds to itself
        return e
    eta = 0 if is_sum else ETA
    g, h = gamma_r(Ih, eta), gamma_a(Ih, eta)
    return reduce(I, e.S * (ONE + Iv(-g, g)), e.D + Iv(-h, h))


def power_of_two(X, e):
    """k when the operand is always 2^k or -2^k and computed exactly."""
    if e == EXACT and X.lo == X.hi:
        k, m = gmpy2.frexp(X.lo)
        if abs(m) == 0.5:
            return k - 1


def sterbenz(Xt, Yt):
    """Whether -y/x is in [1/2, 2] for every x in X~ and y in Y~, so that x + y is exact."""
    a, b = (Xt, -Yt) if Xt.lo > 0 else (-Xt, Yt)
    return a.lo > 0 and b.lo > 0 and b.hi <= 2 * a.lo and a.hi <= 2 * b.lo


def mix(lam, S1, S2):
    """lam S1 + (1 - lam) S2."""
    L = Iv(lam, lam)
    return L * S1 + (ONE - L) * S2


def add(I, X, Y, e1, e2):
    # (x~ + y~) - (x + y) = d1 + d2. Where x, y != 0, (x~ + y~)/(x + y) = lam s1 + (1 - lam) s2
    # with lam = x/(x + y) = 1 - y/(x + y) = 1/(1 + y/x), which is extreme at lam's endpoints.
    S = TOP
    if 0 not in X and 0 not in Y:
        lam = (X / I) & (ONE - Y / I) & (ONE / (ONE + Y / X))
        S = mix(lam.lo, e1.S, e2.S) | mix(lam.hi, e1.S, e2.S)
    Xt, Yt = computed(X, e1), computed(Y, e2)
    return rounded(I, S, e1.D + e2.D, Xt + Yt, exact=sterbenz(Xt, Yt), is_sum=True)


def mul(I, X, Y, e1, e2):
    # x~y~/(xy) = s1 s2 and x~y~ - xy = x d2 + y d1 + d1 d2
    ks = [k for k in (power_of_two(X, e1), power_of_two(Y, e2)) if k is not None]
    D = X * e2.D + Y * e1.D + e1.D * e2.D
    return rounded(I, e1.S * e2.S, D, computed(X, e1) * computed(Y, e2), scale=max(ks, default=None))


def div(I, X, Y, e1, e2):
    # (x~/y~)/(x/y) = s1/s2 and x~/y~ - x/y = (d1 - x(s2 - 1))/y~
    Yt = computed(Y, e2)
    if 0 in Yt:
        return FAIL
    k = power_of_two(Y, e2)
    D = (e1.D - X * (e2.S - ONE)) / Yt
    return rounded(I, e1.S / e2.S, D, computed(X, e1) / Yt, scale=None if k is None else -k)


def sqrt(I, X, e1):
    # sqrt(x~)/sqrt(x) = sqrt(s1) and sqrt(x~) - sqrt(x) = d1/(sqrt(x~) + sqrt(x))
    Xt = computed(X, e1)
    if Xt.lo < 0:
        return FAIL
    return rounded(I, e1.S.sqrt(), e1.D / (Xt.sqrt() + X.sqrt()), Xt.sqrt())


def err(op, I, Xs, es):
    """err_s of an e-node, from its e-class's range I and its children's ranges Xs and pairs es."""
    if FAIL in es:
        return FAIL
    match op:
        case "Var" | "Num":
            return EXACT
        case "Lit":  # a constant the format cannot represent, so it rounds
            return rounded(I, ONE, ZERO)
        case "Neg":
            return Pair(es[0].S, -es[0].D)
        case "Ctx":  # the child on a smaller subdomain, so its class may have a narrower range
            return reduce(I, *es[0])
        case "If":
            return reduce(I, es[0].S | es[1].S, es[0].D | es[1].D)
        case "Add":
            return add(I, *Xs, *es)
        case "Mul":
            return mul(I, *Xs, *es)
        case "Div":
            return div(I, *Xs, *es)
        case "Sqrt":
            return sqrt(I, *Xs, *es)


def abs_bound(e):
    """The readout g_a: |x~ - x| <= max(|lo D|, |hi D|)."""
    return e.D.mag


def rel_bound(e):
    """The readout g_r: |x~ - x|/|x| <= max(|lo S - 1|, |hi S - 1|) when x != 0."""
    return (e.S - ONE).mag


BOUNDS = {"abs": abs_bound, "rel": rel_bound}
