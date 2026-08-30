/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

public import Eig3x3.Basic
public import Eig3x3.Certificates
public import Eig3x3.Eigendecomp
public import Eig3x3.Eigenvalues

/-!
# Eig3x3 — closed-form eigendecomposition of real symmetric 3×3 matrices

Pure Lean 4 eigendecomposition implementation of *real symmetric* 3×3 matrices
over the `Float` type using Habera-Zilian's method for computing the eigenvalue
vector and Eberly's non-iterative method for computing the eigenvector matrix.

No Mathlib, no FFI bindings, native Lean, zero dependencies.

## Why Habera-Zilian's Method?

The algorithms in practical use for this problem are iterative. For example,
LAPACK's symmetric eigensolvers (what NumPy, SciPy, and PyTorch call)
tridiagonalize the matrix and iterate QR steps to convergence, which means
loops, convergence tests, and running time dependent on input random variables.

The pre-existing closed-form alternatives, such as Cardano and Viète’s methods,
avoid this computational cost but suffer in accuracy by evaluating the classical
trigonometric cubic formula whose arccos near ±1 and subtractive discriminant
lose precision two eigenvalues are close.

Habera–Zilian remove this trade-off by proposing a *closed-form* alternative where
numerically dangerous steps — using invariants from diagonal differences, computing
the discriminant as a sum of squares, representing the angle as the arctan — are
re-engineered. Performance is demonstrated to match iterative solvers at machine
precision while still carrying the speed and determinism of a closed-form solution.

## Provenance

* **Eigenvalues** (`Eig3x3.Eigenvalues`) — Habera & Zilian, "Numerically
  stable evaluation of closed-form expressions for eigenvalues of 3×3
  matrices", arXiv:2511.00292v2 (2025). Specifically:
  * the invariants I₁ (Alg. 1), J₂ (Alg. 2), J₃ (Alg. 5),
  * the Algorithm 8 sum-of-squares discriminant (factorization originating in
    Habera–Zilian 2021, arXiv:2111.02117),
  * quadrant-safe angle φ = atan2(√(27Δ), 27J₃) (Eq. 4),
  * ordered eigenvalues λ₀ ≤ λ₁ ≤ λ₂ (Eq. 2).

* **Eigenvectors** (`Eig3x3.Eigenvectors`, internal) — D. Eberly, "A Robust
  Eigensolver for 3×3 Symmetric Matrices", Geometric Tools, CC BY 4.0).
  Specifically, the non-iterative algorithms drescribed in §5, inlcluding:
  * max-abs preconditioning,
  * isolated-eigenvector from cross products (§5, Listing 4),
  * robust orthogonal complement (§5, Listing 5),
  * 2×2 reduction in the complement (§5, Listing 6),
  * right-handed completion.

## Deviations from the sources

* Eberly's acos-based eigenvalue evaluation is replaced by the Habera–Zilian
  invariant pipeline (Habera–Zilian prove the acos form is unstable near repeated
  eigenvalues, and their arctan form is preferred).
* Eberly's sign half-determinant method is replaced by a direct gap comparison
  on the ordered eigenvalues — equivalent in exact arithmetic and more direct.
* Adds a final 3-element sort of the computed eigenvalues: H–Z's ordering
  guarantee (Eq. 2) is exact-arithmetic. Since cosine is a transcendental
  function and exempt from IEEE 754, floating point evaluations of the cosine
  at a degenerate angle can disagree by ~1 ulp. The sort is a pure permutation,
  runs before eigenvector assembly, and the assembly remains right-handed by
  construction (the third vector is always a cross product).
* `eigvecIsolated` adds a defensive exact-zero fallback; Eberly relies on
  exact-arithmetic rank 2 and has no such guard.

## Validation (op-for-op float64 mirror vs `numpy.linalg.eigh`)

* 20k random symmetric matrices:        max |Δλ| ≈ 8e-15
* near-double path diag(−1,1,1+δ):      ≈ 6e-16 for δ ∈ [1e-16, 1e-4]
* near-triple path diag(1,1,1+δ):       ≈ 2e-16 for δ ≤ 1e-8
* double eigenvalue at small scale:     ≈ 1e-19
* exact scaled identity / diagonal:     exact

## Scope

* Real symmetric 3×3 eigendecomposition
* Basic 3×3 `Float` type matrix/vector operations (`Eig3x3.Basic`)
* Runtime certificates for computation assurances (`Eig3x3.Certificates`)

Complex Hermitian is currently *out of scope* but is in consideration for v2.

**Note**: Habera-Zilian report that Algorithm 8 for computing the discriminant
exceeds its (lowest-order) forward-stability bound only for benchmarks with an
ill-conditioned eigenbasis. However, this requires a regime with a non-orthogonal
eigenvector matrix, which symmetric matrices never possess (κ₂ = 1). This package
is scoped to real symmetric matrices, so the failing regime is *out of scope*
by construction.

## Usage

`import Eig3x3` publically offers:
* type primtives (`Vec3`, `Mat3`, `Eigval3`, `SymmMat3`, and `Decomposition`)
* vector/matrix operations through `Vec3` and `Mat3`,
* `eigvals` for calculating eigenvalues using Habera-Zilian's method,
* `eigendecomp` for performing eigendecomposition including the eigenvectors,
* and `certify` for runtime residual error feedback.

Unicode notation for common vector and matrix operations is offered
with `open scoped Eig3x3`, including:
* transpose a matrix: `Qᵀ`
* matrix inverse: `M⁻¹`
* dot product of vectors: `u ⬝ᵥ v`
* Frobenius inner product of matrices: `A ⬝ₘ B`
* outer product of vectors: `u ⊗ᵥ v`
* cross product of vectors: `u ⨯₃ v`
* scalar multiplication: `s • v`, `s • M`
* Hadamard entrywise product: `u ⊙ v`, `A ⊙ B`
* norms: Euclidean `‖v‖`, `‖v‖²` and Frobenius `‖M‖`, `‖M‖²`
* absolute value / magnitude: `|x|`, `|v|`, `|M|` (tight bars: `|v|`, not `| v |`)
* powers: `x ^ⁿ 2`, `M ^ⁿ 2` (and `M ^ 2`)

The eigenvector-only machinery (`Eig3x3.Eigenvectors`) is deliberately not public.
Use `eigendecomp`.
-/
