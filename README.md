# Milo

Milo rewrites a floating-point program into the equivalent program with the smallest sound error bound.

## Setup

Needs [uv] and egglog-experimental, which Milo looks for at `~/.cargo/bin/egglog-experimental` (or `$EGGLOG`):

```
cargo install egglog-experimental --version 3.0.0
```

## Usage

```
uv run milo program.fpcore [--iters N] [--no-branch] [--time-limit S] [--bound-only]
```

prints the programs with the smallest absolute (`abs`) and relative (`rel`) error bounds. `--iters` sets the rewriting
iterations (default 4), `--no-branch` adds no if-statements, and `--time-limit` caps extraction (default 120 s).
`--bound-only` runs a sound error analysis.

The input is an [FPCore] program with:

- `+`, `-`, `*`, `/`, `sqrt`, constants, and `(if (op x c) a b)` with `op` one of `<`, `<=`, `>`, `>=`, `x` an input
  and `c` a constant, not nested on the same input;
- a `:pre` that bounds every input as `(<= lo x hi)`, joined with `and`;
- `:precision binary64` (the default) or `binary32`; inputs are floats in that format.

The program must be defined in real arithmetic on its input ranges. Milo does not check this, and where it fails the
bounds don't hold.

## Evaluation

```
uv run python bench/run.py && uv run python bench/tables.py
```

The paper's Bounding benchmark is `bench/bounding_programs/` and `bench/rewriting_programs/` together, Rewriting is
`bench/rewriting_programs/`, and Branching is `bench/branching_programs/`. `tables.py` prints the numbers the
evaluation's text gives and the data of its figures. Herbie runs as `racket -l herbie`; the other tools are found at the paths in
`bench/tools/`, which environment variables override. We used:

| Tool | Version | Local changes |
|---|---|---|
| [Herbie] | 2.3 (`9297b52`) | platform `bench/tools/herbie_platform.rkt` |
| [FPBench] | `50b1227` | |
| [Daisy] | `6a6f47a` | `Taylor.scala` overrides `max`/`min` to compile on newer JDKs |
| Daisy, for Branching | branch `regimes`, `8e58ff0` | the same fix, upstream's `137379b`, sbt 1.9.9, scala-smtlib from Maven, mpfr_java 1.4 |
| [PRECiSA] | `2626ba1` (v4.2.0) | inputs have no initial error (`initErrorMark`); decimal constants are read exactly (`fpcoreNum2Expr`) |
| [FPTaylor] | `b5a77ca` | |
| [Gappa] | 1.8.2 | |

[uv]: https://docs.astral.sh/uv/
[FPCore]: https://fpbench.org/spec/fpcore-2.0.html
[Herbie]: https://github.com/herbie-fp/herbie
[FPBench]: https://github.com/FPBench/FPBench
[Daisy]: https://github.com/malyzajko/daisy
[FPTaylor]: https://github.com/soarlab/FPTaylor
[Gappa]: https://gappa.gitlabpages.inria.fr/
[PRECiSA]: https://github.com/nasa/PRECiSA
