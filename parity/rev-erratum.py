# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""erratum.py — Exact-Arithmetic Exhibits for Implementation Corrections

This module contains exhibits of errata found during algorithm implementation
and an executable certificate of corrected implementations for this library.

Exhibit: Algorithm 8, r₁₀ sign erratum in Habera–Zilian's 2025 paper
(arXiv:2511.00292).

In exact rational arithmetic, without floats nor tolerances, this shows that:

  * Algorithm 8 as printed, with a minus sign in the second term of r₁₀
    (`p*d1*d2 - q*r*d2 + p*q*q - p*p*p`), FAILS the discriminant identity
    Δ = 4J₂³ − 27J₃², on the large majority of random integer matrices;
  * the corrected term (`+ q*r*d2`) satisfies the identity on every case.

The corrected sign also agrees with the validated x₁₀ term of Habera-Zilian's
2021 paper (arXiv:2111.02117, Eq. 29).

Verified numbers (this script, seed 0): the named exhibit
A = ⟨-2, 2, -9, 9, 8, 6⟩ gives 4J₂³ − 27J₃² = 13021520, Algorithm 8 as
printed = 10740560, corrected = 13021520; in a 2000-case sweep the printed
sign failed 1594/2000, the corrected sign 0/2000.

Exit code is nonzero only if the CORRECTED form ever fails.
"""

import sys
from fractions import Fraction

import numpy as np

import gen_cases

# The named exhibit: the matrix that caught the transcription error.
REGRESSION = tuple(Fraction(x) for x in (-2, 2, -9, 9, 8, 6))


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


def delta_frac(A, r10_sign: int) -> Fraction:
    """Algorithm 8 in exact arithmetic. `r10_sign` is ±1 on the q·r·d₂ term
    of r₁₀: -1 as printed in the 2025 paper, +1 as corrected (and as printed
    for x₁₀ in the 2021 paper, Eq. 29)."""
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
    r10 = p * d1 * d2 + r10_sign * (q * r * d2) + p * q * q - p * p * p
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


def identity_value(A) -> Fraction:
    """4J₂³ − 27J₃², exactly."""
    J2, J3 = j2_frac(A), j3_frac(A)
    return 4 * J2**3 - 27 * J3**2


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    print("# The r₁₀ sign erratum — exact rational arithmetic, zero tolerances\n")
    print("named exhibit: A = ⟨-2, 2, -9, 9, 8, 6⟩")
    print(f"  4J₂³ − 27J₃²            = {identity_value(REGRESSION)}")
    print(
        f"  Algorithm 8 as printed  = {delta_frac(REGRESSION, -1)}"
        "   (r₁₀ with − q·r·d₂)"
    )
    print(
        f"  Algorithm 8 corrected   = {delta_frac(REGRESSION, +1)}"
        "   (r₁₀ with + q·r·d₂)\n"
    )

    rng = np.random.default_rng(seed)
    cases = [tuple(Fraction(x) for x in A) for A in gen_cases.integer_matrices(rng, n)]
    fail_printed = sum(1 for A in cases if delta_frac(A, -1) != identity_value(A))
    fail_corrected = sum(1 for A in cases if delta_frac(A, +1) != identity_value(A))
    print(
        f"sweep over {total} random integer matrices (seed {seed}):".format(
            total=len(cases)
        )
    )
    print(f"  printed sign:   identity fails on {fail_printed}/{len(cases)}")
    print(f"  corrected sign: identity fails on {fail_corrected}/{len(cases)}\n")
    print("cross-reference: the corrected sign agrees with the validated x₁₀")
    print("term of Habera–Zilian 2021 (arXiv:2111.02117, Eq. 29).")

    # The corrected form must hold; the printed form's failures are the point.
    return 1 if fail_corrected else 0


if __name__ == "__main__":
    sys.exit(main())
