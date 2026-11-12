"""FPCore in and out.

A program is a tuple: ("Var", name), ("Num", Fraction), ("Add", a, b), ("Mul", a, b),
("Div", a, b), ("Neg", a), ("Sqrt", a), or ("If", (op, name, Fraction), then, else).
"""

import decimal
import math
import re
from collections import namedtuple
from fractions import Fraction
from itertools import pairwise

Core = namedtuple("Core", "name args box precision body")
SYMBOLS = {"Add": "+", "Mul": "*", "Div": "/", "Neg": "-", "Sqrt": "sqrt"}
FLIP = {">": "<", ">=": "<=", "<": ">", "<=": ">="}


def read(text):
    stack = [[]]
    for token in re.findall(r'"[^"]*"|;[^\n]*|[()\[\]]|[^\s()\[\]";]+', text):
        if token in "([":
            stack.append([])
        elif token in ")]":
            done = stack.pop()
            stack[-1].append(done)
        elif token[0] != ";":
            stack[-1].append(token)
    return stack[0]


def number(token):
    try:
        return Fraction(token)
    except (TypeError, ValueError):
        return None


def expr(s):
    if isinstance(s, str):
        return ("Var", s) if number(s) is None else ("Num", number(s))
    op, *args = s
    if op == "if" and len(args) == 3:
        return ("If", condition(args[0]), expr(args[1]), expr(args[2]))
    args = [expr(a) for a in args]
    if op == "-" and len(args) == 1:
        return ("Neg", args[0])
    if op == "-" and len(args) == 2:
        return ("Add", args[0], ("Neg", args[1]))
    if op == "/" and len(args) == 2:
        return ("Div", *args)
    if op == "sqrt" and len(args) == 1:
        return ("Sqrt", args[0])
    if op in ("+", "*") and len(args) >= 2:
        e = args[0]
        for a in args[1:]:
            e = ("Add" if op == "+" else "Mul", e, a)
        return e
    raise SyntaxError(f"unsupported expression {s}")


def condition(c):
    match c:
        case [op, str(a), str(b)] if op in FLIP and number(a) is None and number(b) is not None:
            return (op, a, number(b))
        case [op, str(a), str(b)] if op in FLIP and number(a) is not None and number(b) is None:
            return (FLIP[op], b, number(a))
    raise SyntaxError(f"unsupported condition {c}")


def variables(e):
    if e[0] == "Var":
        return {e[1]}
    if e[0] == "Num":
        return set()
    if e[0] == "If":
        return {e[1][1]} | variables(e[2]) | variables(e[3])
    return set().union(*map(variables, e[1:]))


def bounds(pre):
    """(name, lo, hi) for each bound in :pre, widened outward to floats."""
    if pre[0] == "and":
        for p in pre[1:]:
            yield from bounds(p)
        return
    if pre[0] != "<=":
        raise SyntaxError(f"unsupported precondition {pre}")
    for a, b in pairwise(pre[1:]):
        if number(a) is not None and number(b) is None:
            q, f = number(a), float(number(a))
            yield b, f if Fraction(f) <= q else math.nextafter(f, -math.inf), math.inf
        elif number(a) is None and number(b) is not None:
            q, f = number(b), float(number(b))
            yield a, -math.inf, f if Fraction(f) >= q else math.nextafter(f, math.inf)
        else:
            raise SyntaxError(f"unsupported precondition {pre}")


def parse(text):
    form = next(f for f in read(text) if isinstance(f, list) and f[:1] == ["FPCore"])
    rest = form[1:]
    name = rest.pop(0) if isinstance(rest[0], str) else None
    args, *props, body = rest
    props = dict(zip(props[::2], props[1::2]))
    box = {x: (-math.inf, math.inf) for x in args}
    for x, lo, hi in bounds(props.get(":pre", ["and"])):
        box[x] = (max(box[x][0], lo), min(box[x][1], hi))
    body = expr(body)
    if not variables(body) <= set(args):
        raise SyntaxError(f"free variables {variables(body) - set(args)}")
    name = props.get(":name", name or "").strip('"')
    return Core(name, args, box, props.get(":precision", "binary64"), body)


def numeral(q):
    if q.denominator == 1:
        return str(q.numerator)
    s = repr(float(q))
    return s if Fraction(s) == q else f"{q.numerator}/{q.denominator}"


def show(e):
    match e:
        case ("Var", x):
            return x
        case ("Num", q):
            return numeral(q)
        case ("Add", a, ("Neg", b)):
            return f"(- {show(a)} {show(b)})"
        case ("If", (op, x, k), a, b):
            return f"(if ({op} {x} {numeral(Fraction(repr(float(k))))}) {show(a)} {show(b)})"
        case (op, *args):
            return f"({SYMBOLS[op]} {' '.join(map(show, args))})"


def fpcore(core, body=None):
    if not all(map(math.isfinite, sum(core.box.values(), ()))):
        raise ValueError(f"unbounded box {core.box}")
    pre = [f"(<= {decimal.Decimal(lo):f} {x} {decimal.Decimal(hi):f})" for x, (lo, hi) in core.box.items()]
    pre = "" if not pre else f" :pre {pre[0]}" if len(pre) == 1 else f" :pre (and {' '.join(pre)})"
    return (f'(FPCore ({" ".join(core.args)}) :name "{core.name}" :precision {core.precision}{pre}'
            f" {show(core.body if body is None else body)})")
