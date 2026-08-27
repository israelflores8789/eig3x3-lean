# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Exact rational-arithmetic identity checks.

Acts as a guard against transcription errors (e.g. the r10 sign error in
Algorithm 8). No floats, no tolerances. Every check here is an exact equality
of Fractions with integer-valued inputs.

Notable checks:

* The headline identity

    delta(A) == 4*J2^3 - 27*J3^2

is computed with the Algorithm 8 sum-of-squares on the left. This is precisely
the identity that caught the r10 transcription error implementation of
Habera-Zilian's Algorithm 8 and has been verified against their 2021 publishing.

* The diagonal-difference forms of J2 and J3 against their direct deviatoric
definitions

    J2 = tr(dev A)^2 / 2, J3 = det(dev A).
"""

import sys
from fractions import Fraction

import numpy as np

import gen_cases


def j2_frac(A) -> Fraction:
    a00, a11, a22, a01, a02, a12 = A
    d0, d1, d2 = a00 - a11, a00 - a22, a11 - a22
    return (d0 * d0 + d1 * d1 + d2 * d2) / 6 + (a01 * a01 + a02 * a02 + a12 * a12)


def j3_frac(A) -> Fraction:
    a00, a11, a22, a01, a02, a12 = A
    d0, d1, d2 = a00 - a11, a00 - a22, a11 - a22
    t1, t2, t3 = d1 + d2, d0 - d2, -d0 - d1
    offdiag = 2 * a01 * a12 * a02
    mixed = (a01 * a01 * t1 + a02 * a02 * t2 + a12 * a12 * t3) / 3
    diag = t1 * t2 * t3 / 27
    return offdiag + mixed - diag


def delta_frac(A) -> Fraction:
    """Habera-Zilian's Algorithm 8 (2025) with r10 correction."""
    a00, a11, a22, p, q, r = A
    d0, d1, d2 = a00 - a11, a00 - a22, a11 - a22
    r1 = p * r * q - q * p * r
    r2 = -p * q * d2 + p * p * r - q * q * r
    r3 = p * r * d1 - p * p * q + q * r * r
    r4 = q * r * d0 + p * r * r - q * q * p
    r5 = p * r * d1 - p * q * p + q * r * r
    r6 = q * r * d0 - p * q * q + p * r * r
    r7 = -q * p * d2 + p * p * r - q * r * q
    r8 = r * d0 * d1 - q * p * d1 + p * p * r - r * r * r
    r9 = r * d0 * d1 - q * p * d0 + q * r * q - r * r * r
    r10 = p * d1 * d2 + q * r * d2 + p * q * q - p * p * p
    r11 = p * d1 * d2 + q * r * d1 + p * r * r - p * p * p
    r12 = -q * d0 * d2 + p * r * d0 + q * r * r - q * q * q
    r13 = q * d0 * d2 + p * r * d2 - p * q * p + q * q * q
    r14 = d0 * d1 * d2 - p * p * d0 + q * q * d1 - r * r * d2
    return (
        9 * r1 * r1
        + 6 * (r2 * r2 + r3 * r3 + r4 * r4)
        + 8 * (r5 * r5 + r6 * r6 + r7 * r7)
        + 2 * (r8 * r8 + r9 * r9 + r10 * r10 + r11 * r11 + r12 * r12 + r13 * r13)
        + r14 * r14
    )


def dev(A):
    """The deviatoric matrix A - (tr A / 3) I."""
    a00, a11, a22, a01, a02, a12 = A
    t = (a00 + a11 + a22) / 3
    return ((a00 - t, a01, a02), (a01, a11 - t, a12), (a02, a12, a22 - t))


def j2_dev(A) -> Fraction:
    """J2 = tr((dev A)^2) / 2"""
    D = dev(A)
    return sum((D[i][j] * D[j][i] for i in range(3) for j in range(3)), Fraction(0)) / 2


def j3_dev(A) -> Fraction:
    """J3 = det(dev A)."""
    (x, p, q), (_, y, r), (_, _, z) = dev(A)
    return x * (y * z - r * r) - p * (p * z - r * q) + q * (p * r - y * q)


def check(n: int, seed: int) -> int:
    rng = np.random.default_rng(seed)
    failures = 0
    for A in gen_cases.integer_matrices(rng, n):
        A = tuple(Fraction(x) for x in A)
        J2, J3 = j2_frac(A), j3_frac(A)
        lhs = 4 * J2**3 - 27 * J3**2
        rhs = delta_frac(A)
        if lhs != rhs:
            print(
                f"FAIL discriminant identity on {A}: "
                f"4J2^3-27J3^2 = {lhs}, delta = {rhs}"
            )
            failures += 1
        if J2 != j2_dev(A):
            print(f"FAIL J2 identity on {A}")
            failures += 1
        if J3 != j3_dev(A):
            print(f"FAIL J3 identity on {A}")
            failures += 1
    print(
        f"exact: {n} integer matrices, {failures} failure(s) "
        f"(discriminant identity, J2 identity, J3 identity)"
    )
    return 1 if failures else 0


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    return check(n, seed)


if __name__ == "__main__":
    sys.exit(main())
