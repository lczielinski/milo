"""Run time: every program of a core compiled into one scalar C loop each, timed on the same
points, one core at a time on an otherwise idle machine."""

import os
import subprocess
import sys
import tempfile

from sample import rounded

TIMES = 10


def literal(q, fmt):
    return float.hex(rounded(q, fmt)) + "f" * (fmt == "binary32")


def c(e, fmt):
    match e:
        case ("Var", x):
            return x
        case ("Num", q):
            return literal(q, fmt)
        case ("Neg", a):
            return f"(- {c(a, fmt)})"
        case ("Sqrt", a):
            return f"{'sqrtf' if fmt == 'binary32' else 'sqrt'}({c(a, fmt)})"
        case ("If", (op, x, k), a, b):
            return f"({x} {op} {literal(k, fmt)} ? {c(a, fmt)} : {c(b, fmt)})"
        case (op, a, b):
            return f"({c(a, fmt)} {dict(Add='+', Mul='*', Div='/')[op]} {c(b, fmt)})"


def measure(core, programs, points):
    """Nanoseconds per evaluation of each named program, the least of a few timings."""
    fmt = "float" if core.precision == "binary32" else "double"
    names = {x: f"v{i}" for i, x in enumerate(core.args)}

    def rename(e):
        if e[0] == "Var":
            return ("Var", names[e[1]])
        if e[0] == "If":
            return ("If", (e[1][0], names[e[1][1]], e[1][2]), rename(e[2]), rename(e[3]))
        return (e[0],) + tuple(rename(a) if isinstance(a, tuple) else a for a in e[1:])

    rows = ",\n".join("{" + ", ".join(float.hex(p[x]) for x in core.args) + "}" for p in points)
    unpack = " ".join(f"{fmt} {names[x]} = X[i][{j}];" for j, x in enumerate(core.args))
    loops = "\n".join(f"""
      start = now(); reps = 0;
      do {{ for (int i = 0; i < N; i++) {{ {unpack} Y[i] = {c(rename(e), core.precision)}; }} escape(); reps++; }}
      while (now() - start < 0.02);
      t = (now() - start) / reps / N * 1e9; if (round == 0 || t < best[{j}]) best[{j}] = t;"""
                      for j, e in enumerate(programs.values()))
    source = f"""#include <math.h>
#include <stdio.h>
#include <time.h>
#define N {len(points)}
static const {fmt} X[N][{max(1, len(core.args))}] = {{{rows}}};
static {fmt} Y[N];
static void escape(void) {{ __asm__ volatile("" : : "r"(Y) : "memory"); }}
static double now(void) {{ struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }}
int main(void) {{
  double start, t, best[{len(programs)}]; long reps;
  for (int round = 0; round < {TIMES}; round++) {{{loops}
  }}
  for (int j = 0; j < {len(programs)}; j++) printf("%g\\n", best[j]);
}}
"""
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "p.c"), "w") as f:
            f.write(source)
        library = os.path.join(d, "p.dylib")
        subprocess.run(["cc", "-O2", "-ffp-contract=off", "-fno-vectorize", "-fno-slp-vectorize", "-shared", "-fPIC",
                        "-o", library, os.path.join(d, "p.c"), "-lm"], check=True, capture_output=True)
        # The loops run as a library loaded by a fresh Python, not as an executable of their own.
        out = subprocess.run([sys.executable, "-c", "import ctypes, sys; ctypes.CDLL(sys.argv[1]).main()", library],
                             check=True, capture_output=True, text=True).stdout
    return dict(zip(programs, map(float, out.split())))
