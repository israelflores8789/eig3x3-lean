module

public import Eig3x3.Basic
public import Eig3x3.Certificates
public import Eig3x3.Eigendecomp
public import Eig3x3.Eigenvalues

/-!
# Eig3x3 — closed-form eigendecomposition of real symmetric 3×3 matrices

Pure Lean 4 (`Float` only): no Mathlib, no FFI, no Python.

## Provenance

* **Eigenvalues** (`Eig3x3.Eigenvalues`) — Habera & Zilian, "Numerically
  stable evaluation of closed-form expressions for eigenvalues of 3×3
  matrices", arXiv:2511.00292v2 (2025): invariants I₁ (Alg. 1), J₂ (Alg. 2),
  J₃ (Alg. 5), the Algorithm 8 sum-of-squares discriminant (factorization
  originating in Habera–Zilian 2021, arXiv:2111.02117), quadrant-safe angle
  φ = atan2(√(27Δ), 27J₃) (Eq. 4), ordered eigenvalues λ₁ ≤ λ₂ ≤ λ₃ (Eq. 2).
  Reference C implementation: `eig3x3` (MIT license).

* **Eigenvectors** (`Eig3x3.Eigenvectors`, internal) — D. Eberly, "A Robust
  Eigensolver for 3×3 Symmetric Matrices", Geometric Tools (documentation:
  CC BY 4.0; code: Boost Software License 1.0), specifically the
  non-iterative `NISymmetricEigensolver3x3`: max-abs preconditioning,
  isolated-first cross products, robust orthogonal complement, 2×2 reduction
  in the complement, right-handed completion.

## Deviations from the sources

* Eberly's acos-based eigenvalue evaluation is replaced by the Habera–Zilian
  invariant pipeline (H–Z prove the acos form is unstable near repeated
  eigenvalues; their arctan form is not).
* Eberly's `sign(halfDet)` branch is replaced by a direct gap comparison on
  the ordered eigenvalues — equivalent in exact arithmetic, but more direct.
* No sorting stage: H–Z returns λ₁ ≤ λ₂ ≤ λ₃, and the eigenvector assembly is
  right-handed by construction (the third vector is always a cross product).
* `eigvecIsolated` adds a defensive exact-zero fallback; Eberly relies on
  exact-arithmetic rank 2 and has no such guard.
* `deltaNaive` is retained solely for benchmarking, to reproduce the paper's
  naive-vs-present comparisons.

## Conditioning scope note

H–Z report that Algorithm 8 exceeds its (lowest-order) forward-stability bound
only for benchmarks with an ill-conditioned eigenbasis — a regime requiring a
non-orthogonal eigenvector matrix, which symmetric matrices never possess
(κ₂ = 1). This library is scoped to real symmetric matrices, so the failing
regime is out of scope by construction.

## Validation (op-for-op float64 mirror vs `numpy.linalg.eigvalsh`)

  * 20k random symmetric matrices:        max |Δλ| ≈ 8e-15
  * near-double path diag(−1,1,1+δ):      ≈ 6e-16 for δ ∈ [1e-16, 1e-4]
  * near-triple path diag(1,1,1+δ):       ≈ 2e-16 for δ ≤ 1e-8
  * double eigenvalue at small scale:     ≈ 1e-19
  * exact scaled identity / diagonal:     exact

## Scope

Real symmetric 3×3 eigendecomposition, runtime certificates
(`Eig3x3.Certificates`), and the 3×3 `Float` matrix/vector vocabulary needed
to consume the result (`Eig3x3.Basic`: basis changes, reconstruction,
componentwise eigenvalue post-processing). General linear algebra beyond
that — including clipping and projections — is downstream's business.
Complex Hermitian is out of scope (cf. Kopp 2008).

## Usage

`import Eig3x3` re-exports the public API: the types, the vector/matrix
operations, `eigvals`, `eigendecomp`, and `certify`. Unicode notation (`Qᵀ`,
`u ⬝ v`) is activated with `open scoped Eig3x3`. The eigenvector machinery
(`Eig3x3.Eigenvectors`) is deliberately not exported.
-/
