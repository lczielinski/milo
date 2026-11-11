"""Build the sound e-graph in egglog, and read back its e-nodes and the range I_c that the range
analysis gives each e-class."""

import os
import subprocess
import tempfile
from collections import namedtuple
from fractions import Fraction
from pathlib import Path

import gmpy2
from gmpy2 import mpfr, mpz

from .fpcore import read
from .interval import DOWN, INF, UP, Iv

EGGLOG = os.environ.get("EGGLOG", os.path.expanduser("~/.cargo/bin/egglog-experimental"))
RULES = os.path.join(os.path.dirname(__file__), "egg")
CONSTRUCTORS = ("Num", "Var", "Add", "Mul", "Div", "Neg", "Sqrt", "If", "Ctx")
TABLES = CONSTRUCTORS + ("lo", "hi", "Bad", "Root")

Node = namedtuple("Node", "op payload children")
EGraph = namedtuple("EGraph", "nodes interval root")


class BadBox(Exception):
    """The range analysis found an empty range: some subexpression is undefined everywhere on the
    input intervals."""


def context(precision, round=gmpy2.RoundToNearest):
    return gmpy2.context(gmpy2.ieee(64 if precision == "binary64" else 32), round=round)


def to_float(q, precision, round=gmpy2.RoundToNearest):
    return context(precision, round).div(mpz(q.numerator), mpz(q.denominator))


def representable(q, precision):
    f = to_float(q, precision)
    return gmpy2.is_finite(f) and Fraction(*f.as_integer_ratio()) == q


def threshold(op, c, precision):
    """Turns x op c into x > k or not x > k."""
    k = to_float(c, precision)
    return float(context(precision).next_below(k) if op in (">=", "<") else k)


def f64(x):
    return repr(float(x)).replace("e+", "e")


def rational(q):
    return f'(bigrat (from-string "{q.numerator}") (from-string "{q.denominator}"))'


def guarded(e, box, precision):
    """Rewrites the input program's ifs to have guards x > k."""
    match e:
        case ("If", (op, x, c), a, b):
            k = threshold(op, c, precision)
            above, below = (a, b) if op in (">", ">=") else (b, a)
            lo, hi = box[x]
            if lo > k or hi <= k:
                return guarded(above if lo > k else below, box, precision)
            arms = guarded(above, box, precision), guarded(below, box, precision)
            return ("If", (">", x, Fraction(k)), *arms)
        case ("Var" | "Num", _):
            return e
        case (op, *args):
            return (op, *(guarded(a, box, precision) for a in args))


def script(core, iters, branch):
    """Builds the egglog program."""
    consts = set()

    def term(e, guard):
        """e on a subdomain, where guard maps a variable to the side of its guard x > k there,
        which becomes a context on that variable."""
        match e:
            case ("Var", x) if x in guard:
                return f'(Ctx "{x}" {guard[x]} (Var "{x}"))'
            case ("Var", x):
                return f'(Var "{x}")'
            case ("Num", q):
                consts.add(q)
                return f"(Num {rational(q)})"
            case ("If", (_, x, k), a, b):
                if x in guard:
                    raise ValueError(f"guards on {x} are nested, which is not supported")
                arms = (term(a, guard | {x: f"{f64(k)} true"}), term(b, guard | {x: f"{f64(k)} false"}))
                return f'(If "{x}" {f64(k)} {arms[0]} {arms[1]})'
            case (op, *args):
                return f"({op} {' '.join(term(a, guard) for a in args)})"

    files = ["analysis.egg", "rules.egg"] + ["branch.egg"] * branch
    lines = [Path(RULES, f).read_text() for f in files]
    lines.append(f"(Root {term(guarded(core.body, core.box, core.precision), {})})")
    for x, (lo, hi) in core.box.items():
        lines += [f'(set (lo (Var "{x}")) {f64(lo)})', f'(set (hi (Var "{x}")) {f64(hi)})']
    for q in consts:
        lo, hi = to_float(q, "binary64", gmpy2.RoundDown), to_float(q, "binary64", gmpy2.RoundUp)
        lines += [f"(set (lo (Num {rational(q)})) {f64(lo)})", f"(set (hi (Num {rational(q)})) {f64(hi)})"]
    # Each iteration pushes contexts to the leaves, runs the range analysis, rewrites once with
    # back-off, folds constants, and branches. With no rewriting, the range analysis runs to a fixed point.
    iteration = "(saturate dist) (repeat 30 ana) (run-with bo opt) (saturate fold) (saturate branch)"
    last = "(repeat 30 ana)" if iters else "(saturate ana)"
    lines.append(f"(run-schedule (let-scheduler bo (back-off)) (repeat {iters} {iteration})"
                 f" (saturate dist) {last})")
    lines += [f"(print-function {t} 1000000000)" for t in TABLES]
    return "\n".join(lines)


def node(op, args, precision):
    """A node from its row in egglog's dump. A constant is a Num if the format represents it, and
    otherwise a Lit, which rounds."""
    match op:
        case "Num":
            _, (_, n), (_, d) = args[0]
            q = Fraction(int(n.strip('"')), int(d.strip('"')))
            return Node("Num" if representable(q, precision) else "Lit", q, ())
        case "Var":
            return Node(op, args[0].strip('"'), ())
        case "If":
            return Node(op, (args[0].strip('"'), Fraction(float(args[1]))), (str(args[2]), str(args[3])))
        case "Ctx":
            payload = (args[0].strip('"'), Fraction(float(args[1])), args[2] == "true")
            return Node(op, payload, (str(args[3]),))
        case _:
            return Node(op, None, tuple(map(str, args)))


def read_dump(text, precision):
    """The e-graph from egglog's print-function output, one table per entry of TABLES."""
    tables = []
    for line in text.splitlines():
        if line == "(":
            tables.append([])
        elif line != ")":
            tables[-1].append(line.split(" -> "))
    nodes, ends, root = {}, {"lo": {}, "hi": {}}, None
    for table, rows in zip(TABLES, tables):
        for lhs, rhs in rows:
            op, *args = read(lhs)[0]
            if table in CONSTRUCTORS:
                nodes.setdefault(str(read(rhs)[0]), []).append(node(op, args, precision))
            elif table in ends:
                ends[table][str(args[0])] = float(rhs)
            elif table == "Bad":
                raise BadBox(lhs)
            else:
                root = str(args[0])
    interval = {}
    for c, ns in nodes.items():
        q = next((n.payload for n in ns if n.op in ("Num", "Lit")), None)
        if q is None:
            interval[c] = Iv(mpfr(ends["lo"].get(c, -INF)), mpfr(ends["hi"].get(c, INF)))
        else:  # the range of a class holding a constant is the 200-bit enclosure of its exact value
            n, d = mpz(q.numerator), mpz(q.denominator)
            interval[c] = Iv(DOWN.div(n, d), UP.div(n, d))
    return EGraph(nodes, interval, root)


def build(core, iters=4, branch=True, timeout=None):
    with tempfile.NamedTemporaryFile("w", suffix=".egg") as f:
        f.write(script(core, iters, branch))
        f.flush()
        run = subprocess.run([EGGLOG, f.name], capture_output=True, text=True, timeout=timeout)
    if run.returncode:
        raise RuntimeError(f"egglog failed: {run.stderr[-1000:]}")
    return read_dump(run.stdout, core.precision)
