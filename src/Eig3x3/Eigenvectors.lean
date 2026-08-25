module

import all Eig3x3.Basic

/-!
# Eig3x3.Eigenvectors — the Eberly eigenvector assembly

D. Eberly, "A Robust Eigensolver for 3×3 Symmetric Matrices", Geometric
Tools (documentation: CC BY 4.0), specifically the non-iterative
algorithm (§5):
- isolated-first cross products (§5, Listing 4),
- robust orthogonal complement (§5, Listing 5),
- 2×2 reduction in the complement (§5, Listing 6),
- right-handed completion.

This module is package-private and not re-exported.
-/

namespace Eig3x3

/-- Eigenvector for the *isolated* eigenvalue λ: the largest cross product of
    rows of (A − λI) (Eberly §5, Listing 4). Variables are kept almost
    verbatim for ease of review against the paper.

    Defensive addition: if all crosses are exactly zero (A = λI up to
    rounding, i.e. a triple eigenvalue), returns e₁ — any unit vector is an
    eigenvector then, and the residual certificate stays tiny because
    ‖A − λI‖ does. -/
def eigvecIsolated (A : SymmMat3) (lam : Float) : Vec3 :=
  let r0 : Vec3 := ⟨A.a00 - lam, A.a01, A.a02⟩
  let r1 : Vec3 := ⟨A.a01, A.a11 - lam, A.a12⟩
  let r2 : Vec3 := ⟨A.a02, A.a12, A.a22 - lam⟩
  let c0 := r0.cross r1
  let c1 := r0.cross r2
  let c2 := r1.cross r2
  let d0 := c0.normSq
  let d1 := c1.normSq
  let d2 := c2.normSq
  let (best, dmax) :=
    if d1 > d0 then
      if d2 > d1 then (c2, d2) else (c1, d1)
    else
      if d2 > d0 then (c2, d2) else (c0, d0)
  if dmax == 0.0 then ⟨1.0, 0.0, 0.0⟩
  else best.scale (1.0 / dmax.sqrt)

/-- Robustly compute U, V so that {U, V, w} is a right-handed orthonormal
    set. Requires `w` unit-length (Eberly §5, Listing 5). -/
def orthonormalComplement (w : Vec3) : Vec3 × Vec3 :=
  let u :=
    if w.y.abs < w.x.abs then
      let inv := 1.0 / (w.x * w.x + w.z * w.z).sqrt
      ⟨-w.z * inv, 0.0, w.x * inv⟩
    else
      let inv := 1.0 / (w.y * w.y + w.z * w.z).sqrt
      ⟨0.0, w.z * inv, -w.y * inv⟩
  (u, w.cross u)

/-- Eigenvector for λ in the plane ⊥ `v0`, where `v0` is the unit eigenvector
    of an adjacent, well-separated eigenvalue. Restricts (A − λI) to the
    plane via J = [U V], then solves the 2×2 null system M X = 0 by
    largest-row selection with division-free normalization. If M ≈ 0
    (repeated eigenvalue), any vector in the plane works, so U is returned
    (Eberly §5, Listing 6). -/
def eigvecInPlane (A : SymmMat3) (v0 : Vec3) (lam : Float) : Vec3 :=
  let (u, v) := orthonormalComplement v0
  let au := A.mulVec u
  let av := A.mulVec v
  let m00 := u.dot au - lam
  let m01 := u.dot av
  let m11 := v.dot av - lam
  if m00.abs < m11.abs then
    -- Solve using row 1: m01·x0 + m11·x1 = 0
    let maxAbs := max m11.abs m01.abs
    if maxAbs == 0.0 then u
    else if m11.abs < m01.abs then
      let t := m11 / m01
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale (t * n)).sub (v.scale n)        -- X = (t, −1)·n
    else
      let t := m01 / m11
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale n).sub (v.scale (t * n))        -- X = (1, −t)·n
  else
    -- Solve using row 0: m00·x0 + m01·x1 = 0
    let maxAbs := max m00.abs m01.abs
    if maxAbs == 0.0 then u
    else if m00.abs < m01.abs then
      let t := m00 / m01
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale n).sub (v.scale (t * n))        -- X = (1, −t)·n
    else
      let t := m01 / m00
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale (t * n)).sub (v.scale n)        -- X = (t, −1)·n

/-- Assemble a right-handed orthonormal eigenbasis for the preconditioned
    matrix `B` from `e`, its ordered eigenvalues as produced by `eigvals B`.
    Column `cᵢ` of the result is a unit eigenvector for `lᵢ`.

    Contract (precondition): `e.l₁ ≤ e.l₂ ≤ e.l₃` must be the spectrum of
    `B` (H–Z Eq. 2 ordering). Isolated-eigenvalue-first (Eberly); the gap
    comparison on the ordered eigenvalues replaces Eberly's `sign(halfDet)`
    branch. -/
def eigvecs (B : SymmMat3) (e : Eigval3) : Mat3 :=
  let gapLo := e.l₂ - e.l₁
  let gapHi := e.l₃ - e.l₂
  if gapLo < gapHi then
    -- λ₃ is isolated (Eberly's halfDet ≥ 0 case)
    let v3 := eigvecIsolated B e.l₃
    let v2 := eigvecInPlane B v3 e.l₂
    ⟨v2.cross v3, v2, v3⟩
  else
    -- λ₁ is isolated (Eberly's halfDet < 0 case)
    let v1 := eigvecIsolated B e.l₁
    let v2 := eigvecInPlane B v1 e.l₂
    ⟨v1, v2, v1.cross v2⟩

end Eig3x3
