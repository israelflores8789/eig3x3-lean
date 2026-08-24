module

public import all Eig3x3.Basic

/-!
# Eig3x3.Certificates — runtime instance checking

The checker half of the certifying-algorithm pair: given any claimed
decomposition of `A` — produced by this library or any other solver —
`certify` computes per-instance evidence (residual, orthonormality,
reconstruction). The checker deliberately uses only `mulVec`, `dot`, and
subtraction — no eigenvalue machinery — so it is trustworthy independently
of the solver's complexity. Since this library intentionally proves no
theorems about `Float`, these certificates are the trust mechanism.

Only `Certificates` and `certify` are public API; the component metrics are
package-private.
-/

namespace Eig3x3

/-- Residual ‖Av − λv‖∞ for a claimed eigenpair. -/
def residual (A : SymmMat3) (lam : Float) (v : Vec3) : Float :=
  let av := A.mulVec v
  let r : Vec3 := ⟨av.x - lam * v.x, av.y - lam * v.y, av.z - lam * v.z⟩
  max3 r.x.abs r.y.abs r.z.abs

/-- Orthonormality certificate: max |cᵢ·cⱼ − δᵢⱼ| over all pairs of columns
    of `Q`. -/
def orthogonalityError (Q : Mat3) : Float :=
  let m1 := max (Q.c₁.dot Q.c₁ - 1.0).abs (Q.c₂.dot Q.c₂ - 1.0).abs
  let m2 := max (Q.c₃.dot Q.c₃ - 1.0).abs (Q.c₁.dot Q.c₂).abs
  let m3 := max (Q.c₁.dot Q.c₃).abs (Q.c₂.dot Q.c₃).abs
  max3 m1 m2 m3

/-- Reconstruct the matrix from a decomposition: A = Σᵢ λᵢ cᵢcᵢᵀ. -/
def reconstruct (d : Decomposition) : SymmMat3 :=
  let comp (m : Float) (v : Vec3) : SymmMat3 :=
    ⟨m * v.x * v.x, m * v.y * v.y, m * v.z * v.z,
     m * v.x * v.y, m * v.x * v.z, m * v.y * v.z⟩
  let b1 := comp d.eigvals.l₁ d.eigvecs.c₁
  let b2 := comp d.eigvals.l₂ d.eigvecs.c₂
  let b3 := comp d.eigvals.l₃ d.eigvecs.c₃
  ⟨b1.a00 + b2.a00 + b3.a00,
   b1.a11 + b2.a11 + b3.a11,
   b1.a22 + b2.a22 + b3.a22,
   b1.a01 + b2.a01 + b3.a01,
   b1.a02 + b2.a02 + b3.a02,
   b1.a12 + b2.a12 + b3.a12⟩

/-- Reconstruction certificate: max |entry| of QΛQᵀ − A. -/
def reconstructionError (A : SymmMat3) (d : Decomposition) : Float :=
  let B := reconstruct d
  let m1 := max (B.a00 - A.a00).abs (B.a11 - A.a11).abs
  let m2 := max (B.a22 - A.a22).abs (B.a01 - A.a01).abs
  let m3 := max (B.a02 - A.a02).abs (B.a12 - A.a12).abs
  max3 m1 m2 m3

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
  let r1 := residual A d.eigvals.l₁ d.eigvecs.c₁
  let r2 := residual A d.eigvals.l₂ d.eigvecs.c₂
  let r3 := residual A d.eigvals.l₃ d.eigvecs.c₃
  { maxResidual := max3 r1 r2 r3,
    orthogonality := orthogonalityError d.eigvecs,
    reconstruction := reconstructionError A d }

end Eig3x3
