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
  let r₀ : Vec3 := ⟨A.a₀₀ - lam, A.a₀₁, A.a₀₂⟩
  let r₁ : Vec3 := ⟨A.a₀₁, A.a₁₁ - lam, A.a₁₂⟩
  let r₂ : Vec3 := ⟨A.a₀₂, A.a₁₂, A.a₂₂ - lam⟩
  let c₀ := r₀.cross r₁
  let c₁ := r₀.cross r₂
  let c₂ := r₁.cross r₂
  let d₀ := c₀.normSq
  let d₁ := c₁.normSq
  let d₂ := c₂.normSq
  let (best, dmax) :=
    if d₁ > d₀ then
      if d₂ > d₁ then (c₂, d₂) else (c₁, d₁)
    else
      if d₂ > d₀ then (c₂, d₂) else (c₀, d₀)
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
    to `v₀`, where `v₀` is the unit eigenvector of an adjacent, well-separated
    eigenvalue. Eigenvectors of a symmetric matrix for distinct eigenvalues are
    orthogonal, and symmetry keeps that plane closed under A. Restricting A − λI
    to the plane via [u v] turns the problem into a 2×2 null system M X = 0,
    whose solution is just the perpendicular of the better-conditioned row, with
    the division arranged so we always divide by the larger coefficient. If the
    whole 2×2 vanishes, such that M ≈ 0, then λ is a repeated eigenvalue, any
    vector in the plane is an eigenvector, and we return `u`.

    Reference: Eberly §5, Listing 6. -/
def eigvecInPlane (A : SymmMat3) (v₀ : Vec3) (lam : Float) : Vec3 :=
  let (u, v) := orthonormalComplement v₀
  let au := A.mulVec u
  let av := A.mulVec v
  let m₀₀ := u.dot au - lam
  let m₀₁ := u.dot av
  let m₁₁ := v.dot av - lam
  if m₀₀.abs < m₁₁.abs then
    -- Solve using row 1: m₀₁·x₀ + m₁₁·x₁ = 0
    let maxAbs := max m₁₁.abs m₀₁.abs
    if maxAbs == 0.0 then u
    else if m₁₁.abs < m₀₁.abs then
      let t := m₁₁ / m₀₁
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale (t * n)).sub (v.scale n)        -- X = (t, −1)·n
    else
      let t := m₀₁ / m₁₁
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale n).sub (v.scale (t * n))        -- X = (1, −t)·n
  else
    -- Solve using row 0: m₀₀·x₀ + m₀₁·x₁ = 0
    let maxAbs := max m₀₀.abs m₀₁.abs
    if maxAbs == 0.0 then u
    else if m₀₀.abs < m₀₁.abs then
      let t := m₀₀ / m₀₁
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale n).sub (v.scale (t * n))        -- X = (1, −t)·n
    else
      let t := m₀₁ / m₀₀
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale (t * n)).sub (v.scale n)        -- X = (t, −1)·n

/-- Assemble a right-handed orthonormal eigenbasis for the preconditioned
    matrix `B` from `e`, its ordered eigenvalues as produced by `eigvals B`.
    Column `cᵢ` of the result is a unit eigenvector for `lᵢ`.

    The cross-product construction is only accurate for an eigenvalue far from
    its neighbors, so we first compare the gaps `l₁ − l₀` and `l₂ − l₁` to find
    the isolated end of the spectrum and compute that eigenvector directly; the
    middle eigenvector then comes from the plane perpendicular to it. The third
    is simply the cross product of the other two since only one perpendicular
    direction remains; thus, it must be the remaining eigenvector, and the cross
    product picks the sign that makes the basis right-handed.

    Precondition (contract): `e.l₀ ≤ e.l₁ ≤ e.l₂` must be the spectrum of `B`.

    References: Eberly §3; Habera-Zilian Eq. 2 for ordering. -/
def eigvecs (B : SymmMat3) (e : Eigval3) : Mat3 :=
  let gapLo := e.l₁ - e.l₀
  let gapHi := e.l₂ - e.l₁
  if gapLo < gapHi then
    -- λ₂ is isolated (Eberly's halfDet ≥ 0 case)
    let v₂ := eigvecIsolated B e.l₂
    let v₁ := eigvecInPlane B v₂ e.l₁
    ⟨v₁.cross v₂, v₁, v₂⟩
  else
    -- λ₀ is isolated (Eberly's halfDet < 0 case)
    let v₀ := eigvecIsolated B e.l₀
    let v₁ := eigvecInPlane B v₀ e.l₁
    ⟨v₀, v₁, v₀.cross v₁⟩

end Eig3x3
