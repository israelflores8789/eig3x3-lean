/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

public import Eig3x3.Basic

/-!
# Eig3x3.Eigenvalues — the Habera–Zilian invariant pipeline

Given a *real symmetric matrix* `A` to decompose, performs Habera-Zilian's
method to compute an ordered vector of eigenvalues, `Eigval3`.

## Algorithm

The eigenvalues of `A` are the roots of its cubic characteristic polynomial.
While cubics have a closed-form trigonometric solution, evaluating it naively
loses accuracy exactly when two eigenvalues are close.

The solution is to sidestep approaching the polynomial directly, and instead
recenter by the mean of the eigenvalues, which the trace `i1` provides trivially,
so the remaining unknowns sum to zero. Then, measure the spread of the recentered
eigenvalues with two moments — `j2`, roughly their variance, and `j3`,
roughly their skewness — computed from *differences* of diagonal entries so
that nearly-equal eigenvalues do not cancel away.

The discriminant, Δ = 4J₂³ − 27J₃², (`delta`) is evaluated as a sum of squares so
no subtraction of nearly-equal eigenvalues can occur and cause floating-point
instability. While the classical formula uses the arccos, the required angle can
be expressed as the arctan2, leveraging its numerical stability near zero. Finally,
the eigenvalues are computed λₖ = (I₁ + 2√(3J₂)·cos(φ/3 + 2πk/3))/3 for
k = { 1, 2, 3 } (`eigvals`).

## Provenance

Habera & Zilian, "Numerically stable evaluation of closed-form expressions
for eigenvalues of 3×3 matrices", arXiv:2511.00292v2 (2025). CC BY 4.0.

Specifically:
* invariants I₁ (Alg. 1), J₂ (Alg. 2), J₃ (Alg. 5),
* the Algorithm 8 sum-of-squares discriminant,
* verified Algorithm 8 through factorization of the
  original Habera-Zilian algorithm
  (Habera–Zilian 2021, Eq. 29, arXiv:2111.02117),
* quadrant-safe angle φ = atan2(√(27Δ), 27J₃) (Eq. 4),
* and ordered eigenvalues λ₁ ≤ λ₂ ≤ λ₃ (Eq. 2).

Deviations:
* Eq. 4's arctan of the ratio is realized as `atan2` of numerator and denominator,
  which is identical when J₃ > 0, quadrant-correct when J₃ < 0, and NaN-free at
  J₃ = Δ = 0.

## Visibility

Exposes `eigvals`.

The rest of this module is internal, package-private.
-/

namespace Eig3x3

open scoped Eig3x3

/-- Algorithm 1: I₁ = tr(A). -/
def i₁ (A : SymmMat3) : Float := A.a₀₀ + A.a₁₁ + A.a₂₂

/-- Algorithm 2, symmetric case: J₂ = ½ tr(dev A)² from diagonal differences
    and off-diagonal squares. Exactly zero for scaled identities (H–Z Eq. 62),
    which is what keeps the near-triple-eigenvalue case stable. -/
def j₂ (A : SymmMat3) : Float :=
  let d₀ := A.a₀₀ - A.a₁₁
  let d₁ := A.a₀₀ - A.a₂₂
  let d₂ := A.a₁₁ - A.a₂₂
  let diag := (d₀ * d₀ + d₁ * d₁ + d₂ * d₂) / 6.0
  let offdiag := (A.a₀₁ * A.a₀₁ + A.a₀₂ * A.a₀₂ + A.a₁₂ * A.a₁₂)
  diag + offdiag

/-- Algorithm 5, symmetric case: J₃ = det(dev A) via diagonal differences. -/
def j₃ (A : SymmMat3) : Float :=
  let d₀ := A.a₀₀ - A.a₁₁
  let d₁ := A.a₀₀ - A.a₂₂
  let d₂ := A.a₁₁ - A.a₂₂
  let t₀ := d₁ + d₂
  let t₁ := d₀ - d₂
  let t₂ := -d₀ - d₁
  let offdiag := 2.0 * A.a₀₁ * A.a₁₂ * A.a₀₂
  let mixed := (A.a₀₁ * A.a₀₁ * t₀ + A.a₀₂ * A.a₀₂ * t₁ + A.a₁₂ * A.a₁₂ * t₂) / 3.0
  let diag := t₀ * t₁ * t₂ / 27.0
  offdiag + mixed - diag

/-- Discriminant Δ = 4J₂³ − 27J₃² = ∏_{i<j}(λᵢ − λⱼ)².

    Habera–Zilian 2025, Algorithm 8, specialized to symmetric A: the auxiliary
    vectors u = DX(A) and v = DX(Aᵀ) coincide, so the sum-of-products
    Δ = Σᵢ wᵢ uᵢ vᵢ becomes a sum of squares with the published weights
    w = (9,6,6,6,8,8,8,2,2,2,2,2,2,1). All 14 DX terms are kept verbatim for
    ease of review against the paper; for symmetric A, r₁ = 0, r₅ = r₃ and
    r₇ = r₂ in exact arithmetic.

    TRANSCRIPTION NOTE: the second term of r₁₀ is `+ q·r·d₂` here. With a
    minus sign there, the formula fails the exact identity Δ = 4J₂³ − 27J₃²
    (verified in exact rational arithmetic over 200 random symmetric
    matrices); the plus sign also agrees with the validated x₁₀ term of the
    2021 paper (arXiv:2111.02117, Eq. 29).

    Validated: exact identity with 4J₂³ − 27J₃² in exact rational arithmetic;
    machine precision on both 2025 benchmark paths (D1, D2), on adversarial
    double-eigenvalue-at-small-scale cases, and on 20k random symmetric
    matrices; exact for scaled identities. -/
def delta (A : SymmMat3) : Float :=
  let p := A.a₀₁
  let q := A.a₀₂
  let r := A.a₁₂
  let d₀ := A.a₀₀ - A.a₁₁
  let d₁ := A.a₀₀ - A.a₂₂
  let d₂ := A.a₁₁ - A.a₂₂
  let r₁  := p*r*q - q*p*r
  let r₂  := -p*q*d₂ + p*p*r - q*q*r
  let r₃  := p*r*d₁ - p*p*q + q*r*r
  let r₄  := q*r*d₀ + p*r*r - q*q*p
  let r₅  := p*r*d₁ - p*q*p + q*r*r
  let r₆  := q*r*d₀ - p*q*q + p*r*r
  let r₇  := -q*p*d₂ + p*p*r - q*r*q
  let r₈  := r*d₀*d₁ - q*p*d₁ + p*p*r - r*r*r
  let r₉  := r*d₀*d₁ - q*p*d₀ + q*r*q - r*r*r
  let r₁₀ := p*d₁*d₂ + q*r*d₂ + p*q*q - p*p*p
  let r₁₁ := p*d₁*d₂ + q*r*d₁ + p*r*r - p*p*p
  let r₁₂ := -q*d₀*d₂ + p*r*d₀ + q*r*r - q*q*q
  let r₁₃ := q*d₀*d₂ + p*r*d₂ - p*q*p + q*q*q
  let r₁₄ := d₀*d₁*d₂ - p*p*d₀ + q*q*d₁ - r*r*d₂
  9.0*r₁*r₁ + 6.0*(r₂*r₂ + r₃*r₃ + r₄*r₄) + 8.0*(r₅*r₅ + r₆*r₆ + r₇*r₇)
    + 2.0*(r₈*r₈ + r₉*r₉ + r₁₀*r₁₀ + r₁₁*r₁₁ + r₁₂*r₁₂ + r₁₃*r₁₃)
    + r₁₄*r₁₄

/-- 2π/3 to full double precision. -/
def twoPiOver3 : Float := 2.0943951023931953

/-- Eigenvalues in increasing order (H–Z Eq. 2 with the Eq. 4 arctan angle).
    `atan2` yields φ ∈ [0, π]; k = 1, 2, 3 then gives λ₁ ≤ λ₂ ≤ λ₃.
    For a scaled identity, J₂ = J₃ = Δ = 0, φ = atan2(0,0) = 0, and all three
    eigenvalues come out as exactly I₁/3.

    Postcondition (contract): the result satisfies `l₀ ≤ l₁ ≤ l₂` enforced by
    the final sort. -/
public def eigvals (A : SymmMat3) : Eigval3 :=
  let I₁ := i₁ A
  let J₂ := j₂ A
  -- Exact scaled identity (H–Z Eq. 62): J₂ is a sum of squares, so it is
  -- exactly zero iff A = cI, and then all three eigenvalues are exactly the
  -- common diagonal entry.
  if J₂ == 0.0 then
    ⟨A.a₀₀, A.a₀₀, A.a₀₀⟩
  else
    let J₃ := j₃ A
    let d  := delta A
    let phi := Float.atan2 (Float.sqrt (27.0 * d)) (27.0 * J₃)
    let c  := 2.0 * Float.sqrt (3.0 * J₂)
    let lam := fun k : Float => (I₁ + c * Float.cos (phi / 3.0 + twoPiOver3 * k)) / 3.0
    let x := lam 1.0
    let y := lam 2.0
    let z := lam 3.0
    -- Exact-arithmetic ordering of Habera-Zilian's Eq. 2 holds only up to
    -- float cosine evaluation at degenerate angles.
    -- Since cosine is a transcendental function, it is exempt from IEEE 754,
    -- and when near Δ = 0 a tied eigenvalue pair can come out ~1 ulp apart
    -- and inverted (caught by the property parity test suite).
    -- A 3-element compare-exchange sort guards the postcondition.
    let (a, b) := if y < x then (y, x) else (x, y)
    let (b, c) := if z < b then (z, b) else (b, z)
    let (a, b) := if b < a then (b, a) else (a, b)
    { l₀ := a, l₁ := b, l₂ := c }

end Eig3x3
