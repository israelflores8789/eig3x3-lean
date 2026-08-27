# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Shared numerical test standard.

Implements one machine epsilon with scale-aware gates.

These formulas mirror Tests/Util.lean exactly — the Lean test suite and this
harness enforce the same standard in the same units.

Gate classes (budget x eps x the quantity's own scale):
  * residual / reconstruction / eigenvalue error: 64 eps * maxAbsEntry(A)
  * orthonormality (also det Q, Q^T Q - I):       16 eps (dimensionless)

The constants come from the algorithm's rounding budget: the eigenvalue
error (<= ~36 eps * maxAbs in validation) dominates, the certificate
evaluation itself contributes <= 5 eps, and <= 4 eps covers platform libm
variation (cos/atan2 are exempt from IEEE 754 1/2 ULP rounding constraint).
"""

EPS = 2.0**-52  # float64 machine epsilon

# A symmetric 3x3 matrix is represented as (a00, a11, a22, a01, a02, a12),
# matching the field order of Eig3x3.SymmMat3.
Symm = tuple


def max_abs_entry(A: Symm) -> float:
    """Largest |entry|; the scale the library preconditions against."""
    a00, a11, a22, a01, a02, a12 = A
    return max(
        max(abs(a00), abs(a01)), max(max(abs(a02), abs(a11)), max(abs(a12), abs(a22)))
    )


def cert_tol(A: Symm) -> float:
    """Residual/reconstruction standard: 64 eps scaled by max |entry|.

    Note: cert_tol of the zero matrix is exactly 0.0, so degenerate cases
    are gated exactly.
    """
    return 64.0 * EPS * max_abs_entry(A)


ORTHO_TOL = 16.0 * EPS  # orthonormality standard: dimensionless (unit vectors)


def eval_tol(A: Symm) -> float:
    """Eigenvalue parity standard: the same 64 eps * maxAbs budget."""
    return cert_tol(A)
