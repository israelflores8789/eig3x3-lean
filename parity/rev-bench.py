# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""bench.py — Performance bench of Lean's Eig3x3 closed form vs numpy/LAPACK

This module measures the computational cost per matrix between Lean's closed form
eigendecomposition and numpy/LAPACK's iterative driver.

Honesty note, printed with every report: the CLI figure is end-to-end
through eig3x3_cli and therefore includes JSON serialization — the real
cost of an IPC-style deployment. A pure in-process compute figure (ns/op)
belongs to a small Lean-side micro-benchmark, which is a separate addition.

[TODO update this with an in-process computation benchmark]

Usage: python parity/bench.py            (from the repo root, after `just cli`)
"""

import argparse
import time

import numpy as np

import compare
import gen_cases
import lean_cli


def time_lean_cli(binary: str, cases: list) -> float:
    """Seconds per matrix, one batch call (process startup amortized)."""
    t0 = time.perf_counter()
    lean_cli.run(cases, binary)
    return (time.perf_counter() - t0) / len(cases)


def time_numpy_loop(cases: list) -> float:
    """Seconds per matrix, per-matrix eigh loop."""
    t0 = time.perf_counter()
    for A in cases:
        np.linalg.eigh(compare.to_dense(A))
    return (time.perf_counter() - t0) / len(cases)


def time_numpy_batched(cases: list) -> float:
    """Seconds per matrix, one vectorized eigh over the stacked batch."""
    M = np.stack([compare.to_dense(A) for A in cases])
    t0 = time.perf_counter()
    np.linalg.eigh(M)
    return (time.perf_counter() - t0) / len(cases)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, default=20000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--lean-binary", default=lean_cli.DEFAULT_BINARY)
    args = p.parse_args()
    if not lean_cli.available(args.lean_binary):
        p.error(f"eig3x3_cli not found at {args.lean_binary!r}; run `just cli` first")

    rng = np.random.default_rng(args.seed)
    cases = gen_cases.random_uniform(rng, args.n)

    print(
        f"# bench — n={args.n} seed={args.seed} "
        f"(CLI figures include the JSON boundary; see module docstring)"
    )
    t_cli = time_lean_cli(args.lean_binary, cases)
    t_loop = time_numpy_loop(cases)
    t_batch = time_numpy_batched(cases)
    print(f"eig3x3_cli, end-to-end:          {t_cli * 1e6:9.2f} us/matrix")
    print(f"numpy/LAPACK eigh, per-matrix:   {t_loop * 1e6:9.2f} us/matrix")
    print(f"numpy/LAPACK eigh, batched:      {t_batch * 1e6:9.2f} us/matrix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
