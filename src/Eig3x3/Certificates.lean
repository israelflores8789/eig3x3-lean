/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

public import Eig3x3.Basic
import all Eig3x3.Basic

/-!
# Eig3x3.Certificates — runtime instance checking

Runtime error validation of decomposition processes. Since this package
intentionally proves no theorems about `Float`, these certificates are
the trust mechanism for downstream computation assurances.

Given any claimed decomposition of `A`, `certify` computes per-instance
evidence:
* residual — how badly each claimed eigenpair violates `Av = λv`
* orthonormality — how far the eigenvectors stray from perpendicular unit vectors
* reconstruction — how closely `QΛQᵀ` reproduces the original `A`

Deliberately uses only `mulVec`, `dot`, and subtraction — no eigenvalue
machinery — so it is trustworthy independently of the solver's complexity.

## Visibility

Exposes `Certificates` and `certify`. The component metrics are package-private.
-/

namespace Eig3x3

open scoped Eig3x3

/-- Residual ‖Av − λv‖∞ for a claimed eigenpair. -/
def residual (A : SymmMat3) (lam : Float) (v : Vec3) : Float :=
  let av := A.mulVec v
  let r : Vec3 := av - lam • v
  max3 |r.x| |r.y| |r.z|

/-- Orthonormality certificate: max |cᵢ·cⱼ − δᵢⱼ| over all pairs of columns
    of `Q`. -/
def orthogonalityError (Q : Mat3) : Float :=
  let m₀ := max |Q.c₀ ⬝ᵥ Q.c₀ - 1.0| |Q.c₁ ⬝ᵥ Q.c₁ - 1.0|
  let m₁ := max |Q.c₂ ⬝ᵥ Q.c₂ - 1.0| |Q.c₀ ⬝ᵥ Q.c₁|
  let m₂ := max |Q.c₀ ⬝ᵥ Q.c₂| |Q.c₁ ⬝ᵥ Q.c₂|
  max3 m₀ m₁ m₂

/-- Reconstruct the matrix from a decomposition: A = Σᵢ λᵢ cᵢcᵢᵀ. -/
def reconstruct (d : Decomposition) : SymmMat3 :=
  let comp (m : Float) (v : Vec3) : SymmMat3 :=
    ⟨m * v.x * v.x, m * v.y * v.y, m * v.z * v.z,
     m * v.x * v.y, m * v.x * v.z, m * v.y * v.z⟩
  let b₀ := comp d.eigvals.l₀ d.eigvecs.c₀
  let b₁ := comp d.eigvals.l₁ d.eigvecs.c₁
  let b₂ := comp d.eigvals.l₂ d.eigvecs.c₂
  ⟨b₀.a₀₀ + b₁.a₀₀ + b₂.a₀₀,
   b₀.a₁₁ + b₁.a₁₁ + b₂.a₁₁,
   b₀.a₂₂ + b₁.a₂₂ + b₂.a₂₂,
   b₀.a₀₁ + b₁.a₀₁ + b₂.a₀₁,
   b₀.a₀₂ + b₁.a₀₂ + b₂.a₀₂,
   b₀.a₁₂ + b₁.a₁₂ + b₂.a₁₂⟩

/-- Reconstruction certificate: max |entry| of QΛQᵀ − A. -/
def reconstructionError (A : SymmMat3) (d : Decomposition) : Float :=
  let B := reconstruct d
  let m₀ := max |B.a₀₀ - A.a₀₀| |B.a₁₁ - A.a₁₁|
  let m₁ := max |B.a₂₂ - A.a₂₂| |B.a₀₁ - A.a₀₁|
  let m₂ := max |B.a₀₂ - A.a₀₂| |B.a₁₂ - A.a₁₂|
  max3 m₀ m₁ m₂

/-- All three certificates for a decomposition of `A`, bundled. -/
public structure Certificates where
  /-- max ‖Avᵢ − λᵢvᵢ‖∞ over the three eigenpairs. -/
  maxResidual : Float
  /-- max |cᵢ·cⱼ − δᵢⱼ| over all column pairs of Q. -/
  orthogonality : Float
  /-- max |entry| of QΛQᵀ − A. -/
  reconstruction : Float
  deriving Repr

/-- Compute all three certificates for a claimed decomposition of `A`.
    Expect ≈ 1e-16 · ‖A‖ or smaller on outputs of `eigendecomp`. -/
public def certify (A : SymmMat3) (d : Decomposition) : Certificates :=
  let r₀ := residual A d.eigvals.l₀ d.eigvecs.c₀
  let r₁ := residual A d.eigvals.l₁ d.eigvecs.c₁
  let r₂ := residual A d.eigvals.l₂ d.eigvecs.c₂
  { maxResidual := max3 r₀ r₁ r₂,
    orthogonality := orthogonalityError d.eigvecs,
    reconstruction := reconstructionError A d }

end Eig3x3
