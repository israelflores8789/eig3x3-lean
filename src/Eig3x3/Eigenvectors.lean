/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import all Eig3x3.Basic

/-!
# Eig3x3.Eigenvectors — the Eberly eigenvector assembly

Given the eigenvalues from Habera-Zilian's method (see `Eig3x3.Eigenvalues`),
this module performs Eberly's non-iterative method to compute the eigenvectors.

## Algorithm

For an eigenvalue λ that is well-separated from its neighbors, the matrix A − λI
has rank 2: its three rows lie in a plane, and the eigenvector is the one
direction perpendicular to that plane computed directly as a cross product
of two rows (`eigvecIsolated`). Cross products lose accuracy as eigenvalues
cluster, so we first compare the gaps λ₂ − λ₁ and λ₃ − λ₂ and apply that
construction only at the isolated end of the spectrum.

Symmetry then does the rest. Eigenvectors of distinct eigenvalues of a symmetric
matrix are orthogonal, so the second eigenvector lies in the plane perpendicular
to the first, and the problem reduces to solving a 2×2 system (`eigvecInPlane`).

Only one perpendicular direction remains, so it is automatically an eigenvector,
unit-length, and right-handed. The third eigenvector is then just the cross
product of the other two. When eigenvalues repeat, any perpendicular direction
is valid, and the code falls back to the unit vector of the eigenbasis.

## Provenance

D. Eberly, "A Robust Eigensolver for 3×3 Symmetric Matrices", Geometric
Tools, (2014). CC BY 4.0.

Specifically, the non-iterative algorithm (§5):
* isolated-eigenvector from cross products (§5, Listing 4),
* robust orthogonal complement (§5, Listing 5),
* 2×2 reduction in the complement (§5, Listing 6),
* right-handed completion.

The gap comparison method replaces Eberly's sign method, which is equivalent
in exact arithmetic but more direct.

## Visbility

This is an internal, package-private module and not intended to be used directly.
-/

namespace Eig3x3

/-- Eigenvector for an *isolated* eigenvalue λ. When λ is well-separated from
    the other eigenvalues, A − λI has rank 2: its three rows lie in a plane,
    and the eigenvector is the one direction perpendicular to that plane. The
    cross product of any two rows points along that perpendicular, so we
    compute all three and keep the longest — in floating point, the longest
    cross carries the most accuracy.

    Reference: Eberly §5, Listing 4.

    Algorithmic Addition: If all crosses are exactly zero then A = λI, every
    direction is an eigenvector, and we return the unit eigenbasis vector e₁. -/
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

/-- Complete a unit vector `w` to a right-handed orthonormal frame {u, v, w}.
    A perpendicular to `w` can be considered `w × eᵢ` for some coordinate axis
    `eᵢ`, and we cross with the axis `w` is least aligned with. The longest
    cross product has the most accurate direction. The third vector is then
    just `w × u`: automatically unit-length and perpendicular to both.

    Reference: Eberly §5, Listing 5. -/
def orthonormalComplement (w : Vec3) : Vec3 × Vec3 :=
  let u :=
    if w.y.abs < w.x.abs then
      let inv := 1.0 / (w.x * w.x + w.z * w.z).sqrt
      ⟨-w.z * inv, 0.0, w.x * inv⟩
    else
      let inv := 1.0 / (w.y * w.y + w.z * w.z).sqrt
      ⟨0.0, w.z * inv, -w.y * inv⟩
  (u, w.cross u)

/-- Eigenvector for λ, exploiting that it must lie in the plane perpendicular
    to `v0`, where `v0` is the unit eigenvector of an adjacent, well-separated
    eigenvalue. Eigenvectors of a symmetric matrix for distinct eigenvalues are
    orthogonal, and symmetry keeps that plane closed under A. Restricting A − λI
    to the plane via [u v] turns the problem into a 2×2 null system M X = 0,
    whose solution is just the perpendicular of the better-conditioned row, with
    the division arranged so we always divide by the larger coefficient. If the
    whole 2×2 vanishes, such that M ≈ 0, then λ is a repeated eigenvalue, any
    vector in the plane is an eigenvector, and we return `u`.

    Reference: Eberly §5, Listing 6. -/
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

    The cross-product construction is only accurate for an eigenvalue far from
    its neighbors, so we first compare the gaps `l₂ − l₁` and `l₃ − l₂` to find
    the isolated end of the spectrum and compute that eigenvector directly; the
    middle eigenvector then comes from the plane perpendicular to it. The third
    is simply the cross product of the other two since only one perpendicular
    direction remains; thus, it must be the remaining eigenvector, and the cross
    product picks the sign that makes the basis right-handed.

    Precondition (contract): `e.l₁ ≤ e.l₂ ≤ e.l₃` must be the spectrum of `B`.

    References: Eberly §3; Habera-Zilian Eq. 2 for ordering. -/
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
