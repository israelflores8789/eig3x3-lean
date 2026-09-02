/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

public import Eig3x3.Basic
import all Eig3x3.Basic
import Eig3x3.Eigenvalues
import all Eig3x3.Eigenvectors

/-!
# Eig3x3.Eigendecomp — the pipeline driver

`eigendecomp` is the package's primary entry point. It owns all steps of
the eigendecomposition pipeline:
- zero fast path,
- max-abs preconditioning (Eberly's overflow guard [Eberly2014]),
- ordered eigenvalues on the scaled matrix (Habera–Zilian [HaberaZilian2025]),
- isolated-eigenvalue-first eigenvector assembly (Eberly [Eberly2014]),
- and eigenvalue rescaling.
-/

namespace Eig3x3

/-- Full eigendecomposition of a real symmetric 3×3 matrix:
    `A = QΛQᵀ = Σᵢ λᵢ cᵢcᵢᵀ` with `l₀ ≤ l₁ ≤ l₂` and `Q` right-handed
    orthonormal.

    Pipeline: zero fast path → max-abs preconditioning (Eberly [Eberly2014]) → ordered
    eigenvalues on the scaled matrix (Habera–Zilian [HaberaZilian2025]) → eigenvector
    assembly consuming those eigenvalues (Eberly [Eberly2014]) → eigenvalue rescaling. -/
public def eigendecomp (A : SymmMat3) : Decomposition :=
  let maxAbs := A.maxAbsEntry
  if maxAbs == 0.0 then
    { eigvals := ⟨0.0, 0.0, 0.0⟩, eigvecs := Mat3.id }
  else
    -- Bit-transparent preconditioning: `frExp` gives `maxAbs = m · 2^e` with
    -- m ∈ [0.5, 1) exactly, so `m / maxAbs` is exactly 2^(-e) — the true
    -- quotient is a power of two and IEEE division is correctly rounded.
    -- Scaling by a power of two moves only exponents, so the whole pipeline
    -- is bit-identical to running unpreconditioned, and the guard costs zero
    -- rounding.
    -- Contrast: `1.0 / maxAbs` only perturbs ~96% of inputs by ~1 ULP.
    let (m, _) := maxAbs.frExp
    let B := A.scale (m / maxAbs)
    let e := eigvals B
    let Q := eigvecs B e
    -- `maxAbs / m` is exactly 2^e > 0, so rescaling is exact and order-preserving.
    { eigvals := e.scale (maxAbs / m), eigvecs := Q }

end Eig3x3
