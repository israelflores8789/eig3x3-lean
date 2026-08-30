/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

public import Eig3x3
import all Eig3x3.Basic

/-!
# Tests.Util — assertion helpers, contract predicates, and the case zoo

Shared infrastructure for the test suites: float-tolerance assertions,
approximate-equality helpers, the executable contract predicates (ordering,
right-handedness), and the matrices exercised by every suite. Keeping the
zoo here gives all suites a single source of truth for test inputs.
-/

namespace Eig3x3

/-! ## Approximate equality and contract predicates -/

/-- Largest |component| of a vector. Test-side helper for scale-aware gates. -/
public def Vec3.maxAbs (v : Vec3) : Float :=
  max3 v.x.abs v.y.abs v.z.abs

/-- Entrywise approximate equality for vectors. -/
public def Vec3.approx (u v : Vec3) (tol : Float) : Bool :=
  decide ((u.x - v.x).abs ≤ tol)
    && decide ((u.y - v.y).abs ≤ tol)
    && decide ((u.z - v.z).abs ≤ tol)

/-- Entrywise approximate equality for matrices. -/
public def Mat3.approx (M N : Mat3) (tol : Float) : Bool :=
  M.c₀.approx N.c₀ tol && M.c₁.approx N.c₁ tol && M.c₂.approx N.c₂ tol

/-- Ordering contract on `eigvals`: `l₀ ≤ l₁ ≤ l₂`. -/
public def Eigval3.isOrdered (e : Eigval3) : Bool :=
  decide (e.l₀ ≤ e.l₁) && decide (e.l₁ ≤ e.l₂)

namespace Tests

/-! ## Machine Precision Constants -/

/-- float64 machine epsilon, 2⁻⁵². -/
-- 0x1p-52 in idomatic Lean
public def Float.eps : Float := (1 : Float) / ((2 : Float) ^ 52)

/-- The certificate standard for residual and reconstruction: 64ε scaled by
    the matrix's max |entry|. See the Certificates suite's module docstring
    for the budget derivation. -/
public def certTol (A : SymmMat3) : Float := 64.0 * Float.eps * A.maxAbsEntry

/-- The orthonormality standard: 16ε, dimensionless (unit vectors). -/
public def orthoTol : Float := 16.0 * Float.eps

/-! ## Assertions -/

/-- Fail with `name` if `|x| > tol`; otherwise report the value. -/
public def assertNear (name : String) (x tol : Float) : IO Unit := do
  if x.abs > tol then
    throw (IO.userError s!"FAIL {name}: |{x}| > {tol}")
  else
    IO.println s!"ok   {name}: {x}"

/-- Fail with `name` if `x` and `y` differ by more than `tol`.
    `tol = 0.0` demands exact equality. -/
public def assertClose (name : String) (x y tol : Float) : IO Unit :=
  assertNear name (x - y) tol

/-- Fail with `name` unless `b` holds. -/
public def assertTrue (name : String) (b : Bool) : IO Unit := do
  unless b do
    throw (IO.userError s!"FAIL {name}")

/-! ## The case zoo -/

/-- [[2,1,0],[1,2,1],[0,1,2]] — exact eigenvalues 2−√2, 2, 2+√2. -/
public def workedExample : SymmMat3 := ⟨2.0, 2.0, 2.0, 1.0, 0.0, 1.0⟩

/-- Diagonal with distinct entries: eigenvalues (1, 2, 3) with axis
    eigenvectors (up to sign). -/
public def diagonalCase : SymmMat3 := ⟨1.0, 2.0, 3.0, 0.0, 0.0, 0.0⟩

/-- Scaled identity: exact triple eigenvalue 3.7; exercises the
    all-cross-products-zero fallback in `eigvecIsolated`. -/
public def scaledIdentity : SymmMat3 := ⟨3.7, 3.7, 3.7, 0.0, 0.0, 0.0⟩

/-- The zero matrix: exercises the exact fast path in `eigendecomp`. -/
public def zeroMatrix : SymmMat3 := ⟨0.0, 0.0, 0.0, 0.0, 0.0, 0.0⟩

/-- Near-double eigenvalue: eigenvalues ≈ (−1, 1−1e-8, 1+1e-8). Exercises
    the clustered-top gap branch; with the Algorithm 8 discriminant this
    achieves machine precision (the naive Δ gave ≈2.5e-9 here). -/
public def nearDouble : SymmMat3 := ⟨1.0, 1.0, -1.0, 1.0e-8, 0.0, 0.0⟩

/-- Near-triple eigenvalue path (an H–Z benchmark path): `diag(1, 1, 1 + δ)`
    with δ = 1e-8 — an exact double eigenvalue plus a close third, driving
    the angle φ → 0 (J₃ > 0). -/
public def nearTriple : SymmMat3 := ⟨1.0, 1.0, 1.0 + 1.0e-8, 0.0, 0.0, 0.0⟩

/-- The matrix that caught the r₁₀ transcription error during porting:
    exact Δ = 13,021,520 (exactly representable in float64). -/
public def regressionMatrix : SymmMat3 := ⟨-2.0, 2.0, -9.0, 9.0, 8.0, 6.0⟩

end Tests
end Eig3x3
