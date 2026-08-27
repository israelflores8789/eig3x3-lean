# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Property checks: algebraic invariants that must hold on every input.

These run against the mirror (or, later, the CLI) — never against numpy.
The bit-exact scale-invariance property in particular is ours alone: the
power-of-two preconditioner makes eigendecomp(2^k A) bit-identical in
significand to 2^k eigendecomp(A), so that check is an exact equality, not
a tolerance.

Usage: python properties.py [n] [seed]
"""

import math
import sys

import numpy as np

import gen_cases
import mirror
from gates import EPS, ORTHO_TOL, cert_tol, max_abs_entry


def det3(A) -> float:
    a00, a11, a22, a01, a02, a12 = A
    return (
        a00 * (a11 * a22 - a12 * a12)
        - a01 * (a01 * a22 - a02 * a12)
        + a02 * (a01 * a12 - a02 * a11)
    )


def check_case(A, dec, failures) -> None:
    e, Q, c = dec

    # Ordering contract: l1 <= l2 <= l3, always, exactly.
    if not (e[0] <= e[1] <= e[2]):
        failures.append(("ordering", A, e))

    # Certificates under the shared gates.
    if c["maxResidual"] > cert_tol(A):
        failures.append(("residual", A, c["maxResidual"]))
    if c["orthogonality"] > ORTHO_TOL:
        failures.append(("orthogonality", A, c["orthogonality"]))
    if c["reconstruction"] > cert_tol(A):
        failures.append(("reconstruction", A, c["reconstruction"]))

    # Trace invariant: sum l_i == tr A, within the certificate budget
    # (the sum itself adds ~2 eps of rounding).
    tr = (A[0] + A[1]) + A[2]
    s = (e[0] + e[1]) + e[2]
    if abs(s - tr) > cert_tol(A):
        failures.append(("trace", A, abs(s - tr)))

    # Determinant invariant: prod l_i == det A. Scale is maxAbs^3; the
    # budget (64) covers det evaluation (~5 eps), the product (~3 eps),
    # and the eigenvalue errors (~36 eps) at that scale.
    mA = max_abs_entry(A)
    prod = (e[0] * e[1]) * e[2]
    if abs(prod - det3(A)) > 64.0 * EPS * mA * mA * mA:
        failures.append(("det", A, abs(prod - det3(A))))

    # Right-handedness: det Q == +1 within the (dimensionless) ortho budget.
    if abs(mirror.mat_det(Q) - 1.0) > ORTHO_TOL:
        failures.append(("handedness", A, mirror.mat_det(Q)))


def check_scale_invariance(A, dec, failures) -> None:
    """eigendecomp(2^k A).eigvals == 2^k * eigendecomp(A).eigvals, bitwise.

    Only meaningful away from the subnormal/overflow frontier, so restricted
    to moderate scales (the logscale suite covers the extremes separately).
    """
    mA = max_abs_entry(A)
    if not (1e-100 <= mA <= 1e100):
        return
    e = dec[0]
    for k in (-2, -1, 1, 2):
        As = tuple(math.ldexp(x, k) for x in A)
        es = mirror.eigendecomp(As)[0]
        want = tuple(math.ldexp(x, k) for x in e)
        if es != want:
            failures.append(("scale-invariance", A, k, es, want))


def run(n: int, seed: int) -> int:
    rng = np.random.default_rng(seed)
    cases = (
        list(gen_cases.ZOO.values())
        + gen_cases.random_uniform(rng, n)
        + gen_cases.adversarial_paths()
    )
    failures = []
    for A in cases:
        dec = mirror.decomposition(A)
        check_case(A, dec, failures)
        check_scale_invariance(A, dec, failures)
    by_kind = {}
    for f in failures:
        by_kind[f[0]] = by_kind.get(f[0], 0) + 1
    print(
        f"properties: {len(cases)} cases, failures by kind: "
        f"{by_kind if by_kind else 'none'}"
    )
    for f in failures[:5]:
        print("  FAIL", f)
    return 1 if failures else 0


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    return run(n, seed)


if __name__ == "__main__":
    sys.exit(main())
