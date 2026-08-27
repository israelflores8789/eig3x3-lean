# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Benchmarks and the naive-vs-present discriminant study.

Timing caveat: the mirror is pure Python and is not the artifact under
test — its per-matrix time says nothing about the Lean binary. The
meaningful timing comparison (eig3x3_cli, which should be flat constant
time, vs numpy/LAPACK, which is O(n^3) with an input-dependent iteration
count) lands with the CLI. The accuracy study below is the point today.

The accuracy study reproduces the paper's naive-vs-present comparison:
the same pipeline run with the Algorithm 8 discriminant vs the clamped
naive discriminant, on the adversarial paths whose exact eigenvalues are
known in closed form.
"""

import sys
import time

import numpy as np

import compare
import gen_cases
import mirror


def time_per_matrix(fn, cases: list, repeat: int = 3) -> float:
    """Best-of-repeat seconds per matrix."""
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        for A in cases:
            fn(A)
        best = min(best, (time.perf_counter() - t0) / len(cases))
    return best


def naive_vs_present() -> None:
    """Error of the two discriminants on the adversarial paths.

    The path matrices have exactly-known float64 eigenvalues, so no
    high-precision reference is needed.
    """
    print(f"{'delta':>10} {'path':<18} {'err (Alg. 8)':>14} {'err (naive)':>14}")
    for d in gen_cases.DELTAS:
        paths = [
            (
                "near-double-top",
                (-1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0),
                (-1.0, 1.0, 1.0 + d),
            ),
            ("near-triple", (1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0 + d)),
            ("offdiag-double", (1.0, 1.0, -1.0, d, 0.0, 0.0), (-1.0, 1.0 - d, 1.0 + d)),
        ]
        for label, A, ref in paths:
            e_p = mirror.eigendecomp(A)[0]
            e_n = mirror.eigendecomp(A, use_naive=True)[0]
            err_p = max(abs(e_p[i] - ref[i]) for i in range(3))
            err_n = max(abs(e_n[i] - ref[i]) for i in range(3))
            print(f"{d:>10.1e} {label:<18} {err_p:>14.3e} {err_n:>14.3e}")


def main() -> int:
    rng = np.random.default_rng(0)
    cases = gen_cases.random_uniform(rng, 1000)

    print("== timing (see module docstring for the caveat) ==")
    t_mirror = time_per_matrix(mirror.eigendecomp, cases[:200])
    print(f"mirror (pure Python): {t_mirror * 1e6:10.1f} us/matrix")
    t_numpy = time_per_matrix(lambda A: compare.numpy_decomposition(A), cases)
    print(f"numpy/LAPACK:         {t_numpy * 1e6:10.1f} us/matrix")

    print("\n== naive-vs-present discriminant ==")
    naive_vs_present()
    return 0


if __name__ == "__main__":
    sys.exit(main())
