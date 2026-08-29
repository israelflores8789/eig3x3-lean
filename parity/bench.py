# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""bench.py — Performance bench of Lean's Eig3x3 closed form vs numpy/LAPACK

This module measures the computational cost per matrix between Lean's closed form
eigendecomposition and numpy/LAPACK's iterative driver.

Reports:
* in-process Lean per matrix computation cost
* in-process Lean per matrix computation cost with accuracy certificates
* end-to-end Lean per matrix computation cost, which includes JSON
  serialization as a model of IPC-style deployment
* numpy/LAPACK per matrix in-process computation cost
* numpy/LAPACK batched in-process computation cost

Usage:
    python bench.py                      # after `just build_bench` and `just build_cli`
"""

import argparse
import os
import subprocess
import time

import numpy as np

import gen_cases
import lean_cli
import parity

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BENCH_BINARY = os.path.join(PROJECT_ROOT, ".lake", "build", "bin", "eig3x3_bench")


def time_lean_ops(bench_binary: str) -> tuple[float, float]:
    """Pure in-process µs/matrix from the Lean micro-benchmark:
    (eigendecomp only, eigendecomp + certify)."""
    proc = subprocess.run([bench_binary], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"eig3x3_bench exited {proc.returncode}: {proc.stderr.strip()}"
        )
    vals = {}
    for line in proc.stdout.splitlines():
        if line.startswith("ns_per_op_"):
            key, _, v = line.partition(" ")
            vals[key] = float(v) / 1000.0  # ns -> µs
    if "ns_per_op_decomp" not in vals or "ns_per_op_full" not in vals:
        raise RuntimeError("eig3x3_bench output missing ns_per_op lines")
    return vals["ns_per_op_decomp"], vals["ns_per_op_full"]


def time_lean_cli(binary: str, cases: list) -> float:
    """Seconds per matrix, one batch call (process startup amortized)."""
    t0 = time.perf_counter()
    lean_cli.run(cases, binary)
    return (time.perf_counter() - t0) / len(cases)


def time_numpy_loop(cases: list) -> float:
    """Seconds per matrix, per-matrix eigh loop."""
    t0 = time.perf_counter()
    for A in cases:
        np.linalg.eigh(parity.to_dense(A))
    return (time.perf_counter() - t0) / len(cases)


def time_numpy_batched(cases: list) -> float:
    """Seconds per matrix, one vectorized eigh over the stacked batch."""
    M = np.stack([parity.to_dense(A) for A in cases])
    t0 = time.perf_counter()
    np.linalg.eigh(M)
    return (time.perf_counter() - t0) / len(cases)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, default=20000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--lean-binary", default=lean_cli.DEFAULT_BINARY)
    args = p.parse_args()

    if not lean_cli.available(BENCH_BINARY):
        p.error(
            f"eig3x3_bench not found at {BENCH_BINARY!r}; run `just build_bench` first",
        )

    if not lean_cli.available(args.lean_binary):
        p.error(
            f"eig3x3_cli not found at {args.lean_binary!r}; run `just build_cli` first"
        )

    rng = np.random.default_rng(args.seed)
    cases = gen_cases.random_uniform(rng, args.n)

    print(
        f"# bench — n={args.n} seed={args.seed} "
        f"(CLI figures include the JSON boundary; see module docstring)"
    )
    t_bench_decomp, t_bench_full = time_lean_ops(BENCH_BINARY)
    t_cli = time_lean_cli(args.lean_binary, cases)
    t_loop = time_numpy_loop(cases)
    t_batch = time_numpy_batched(cases)
    print(f"eig3x3_cli, in-process:            {t_bench_decomp:9.2f} us/matrix")
    print(f"eig3x3_cli, in-process, w/ certs:  {t_bench_full:9.2f} us/matrix")
    print(f"eig3x3_cli, end-to-end:            {t_cli * 1e6:9.2f} us/matrix")
    print(f"numpy/LAPACK eigh, per-matrix:     {t_loop * 1e6:9.2f} us/matrix")
    print(f"numpy/LAPACK eigh, batched:        {t_batch * 1e6:9.2f} us/matrix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
