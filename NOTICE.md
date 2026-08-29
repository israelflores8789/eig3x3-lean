# NOTICE

**eig3x3-lean** — closed-form eigenvalue and eigenvector decomposition of 3×3 real
matrices, implemented in pure Lean 4 (no FFI, no external dependencies).

Copyright © 2026 Israel Flores-Arbolay (<https://orcid.org/0009-0009-6801-4103>)

This project is distributed under the terms of its own license; see `LICENSE`.
This file provides legally and scholarly relevant attribution for the works this
project builds on. Mathematical methods are not subject to copyright; this attribution 
is maintained as scholarly courtesy.

---

## 1. Algorithmic attribution

The numerical methods implemented in this library are derived from the following
sources. If you use this library in academic work, please cite them; see
`CITATION.cff` for machine-readable records.

### Eigenvalue Provenance
Michal Habera and Andreas Zilian, "Numerically stable evaluation of closed-form
expressions for eigenvalues of 3×3 matrices," arXiv:2511.00292 [math.NA], 2025.
<https://doi.org/10.48550/arXiv.2511.00292>

Michal Habera and Andreas Zilian, "Symbolic spectral decomposition of 3x3 matrices,"
arXiv:2111.02117 [math.NA], 2021. <https://doi.org/10.48550/arXiv.2111.02117>

Closed-form evaluation of the eigenvalues of 3×3 matrices via the trace and
deviatoric invariants with a stabilized discriminant as defined in the 2025 paper
(unless noted otherwise), specifically:
* invariants I₁ (Alg. 1), J₂ (Alg. 2), J₃ (Alg. 5),
* the Algorithm 8 sum-of-squares discriminant,
* verified Algorithm 8 through factorization of the original Habera-Zilian algorithm
  (Habera–Zilian 2021, Eq. 29),
* quadrant-safe angle φ = atan2(√(27Δ), 27J₃) (Eq. 4),
* and ordered eigenvalues λ₁ ≤ λ₂ ≤ λ₃ (Eq. 2).

### Eigenvector Provenance
David H. Eberly, "A Robust Eigensolver for 3x3 Symmetric Matrices," Geometric Tools, 
LLC, 2014. <https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf>

Null-space construction of eigenvectors from known eigenvalues, specifically
the non-iterative algorithm in §5:
* max-abs preconditioning,
* isolated-eigenvector from cross products (§5, Listing 4),
* robust orthogonal complement (§5, Listing 5),
* 2×2 reduction in the complement (§5, Listing 6),
* right-handed completion method.

## 2. AI assistance disclosure

This library was developed with AI coding assistants (Kimi K3, Gemini 3.7 Flash)
used for code generation under continuous human direction. Every file was
manually reviewed and edited, and the mathematics was checked against the sources
listed in §1.

Correctness does not rest on manual review. The implementation is parity-tested
in CI against an independent NumPy/LAPACK reference over a shared corpus of test
vectors — including repeated and near-repeated eigenvalues, mpmath-based "golden"
vectors, degenerate discriminant regimes, and runtime certificates on residual
error, orthonomality, and reconstruction (QᵀQ = I) — under machine-epsilon-scaled 
tolerances.

See `tests/` and `parity/`, or run `just parity` to reproduce.

## 3. Name and Affiliation

The project name references `eig3x3`, the C implementation by Habera and Zilian, 
and marks this as the Lean 4 implementation of their published work. This 
project is not affiliated with or endorsed by the original authors.
