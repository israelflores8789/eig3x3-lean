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
no subtraction of nearly-equal eigenvalues can occur and cause float-point
instability. While the classical formula uses the arccos, the required angle can
be expressed as the arctan2 leveraging its numerical stability near zero. Finally,
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

/-- Algorithm 1: I₁ = tr(A). -/
def i1 (A : SymmMat3) : Float := A.a00 + A.a11 + A.a22

/-- Algorithm 2, symmetric case: J₂ = ½ tr(dev A)² from diagonal differences
    and off-diagonal squares. Exactly zero for scaled identities (H–Z Eq. 62),
    which is what keeps the near-triple-eigenvalue case stable. -/
def j2 (A : SymmMat3) : Float :=
  let d0 := A.a00 - A.a11
  let d1 := A.a00 - A.a22
  let d2 := A.a11 - A.a22
  let diag := (d0 * d0 + d1 * d1 + d2 * d2) / 6.0
  let offdiag := (A.a01 * A.a01 + A.a02 * A.a02 + A.a12 * A.a12)
  diag + offdiag

/-- Algorithm 5, symmetric case: J₃ = det(dev A) via diagonal differences. -/
def j3 (A : SymmMat3) : Float :=
  let d0 := A.a00 - A.a11
  let d1 := A.a00 - A.a22
  let d2 := A.a11 - A.a22
  let t1 := d1 + d2
  let t2 := d0 - d2
  let t3 := -d0 - d1
  let offdiag := 2.0 * A.a01 * A.a12 * A.a02
  let mixed := (A.a01 * A.a01 * t1 + A.a02 * A.a02 * t2 + A.a12 * A.a12 * t3) / 3.0
  let diag := t1 * t2 * t3 / 27.0
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
  let p := A.a01
  let q := A.a02
  let r := A.a12
  let d0 := A.a00 - A.a11
  let d1 := A.a00 - A.a22
  let d2 := A.a11 - A.a22
  let r1  := p*r*q - q*p*r
  let r2  := -p*q*d2 + p*p*r - q*q*r
  let r3  := p*r*d1 - p*p*q + q*r*r
  let r4  := q*r*d0 + p*r*r - q*q*p
  let r5  := p*r*d1 - p*q*p + q*r*r
  let r6  := q*r*d0 - p*q*q + p*r*r
  let r7  := -q*p*d2 + p*p*r - q*r*q
  let r8  := r*d0*d1 - q*p*d1 + p*p*r - r*r*r
  let r9  := r*d0*d1 - q*p*d0 + q*r*q - r*r*r
  let r10 := p*d1*d2 + q*r*d2 + p*q*q - p*p*p
  let r11 := p*d1*d2 + q*r*d1 + p*r*r - p*p*p
  let r12 := -q*d0*d2 + p*r*d0 + q*r*r - q*q*q
  let r13 := q*d0*d2 + p*r*d2 - p*q*p + q*q*q
  let r14 := d0*d1*d2 - p*p*d0 + q*q*d1 - r*r*d2
  9.0*r1*r1 + 6.0*(r2*r2 + r3*r3 + r4*r4) + 8.0*(r5*r5 + r6*r6 + r7*r7)
    + 2.0*(r8*r8 + r9*r9 + r10*r10 + r11*r11 + r12*r12 + r13*r13)
    + r14*r14

/-- 2π/3 to full double precision. -/
def twoPiOver3 : Float := 2.0943951023931953

/-- Eigenvalues in increasing order (H–Z Eq. 2 with the Eq. 4 arctan angle).
    `atan2` yields φ ∈ [0, π]; k = 1, 2, 3 then gives λ₁ ≤ λ₂ ≤ λ₃.
    For a scaled identity, J₂ = J₃ = Δ = 0, φ = atan2(0,0) = 0, and all three
    eigenvalues come out as exactly I₁/3.

    Postcondition (contract): the result satisfies `l₁ ≤ l₂ ≤ l₃`. -/
public def eigvals (A : SymmMat3) : Eigval3 :=
  let I1 := i1 A
  let J2 := j2 A
  -- Exact scaled identity (H–Z Eq. 62): J₂ is a sum of squares, so it is
  -- exactly zero iff A = cI, and then all three eigenvalues are exactly the
  -- common diagonal entry.
  if J2 == 0.0 then
    ⟨A.a00, A.a00, A.a00⟩
  else
    let J3 := j3 A
    let d  := delta A
    let phi := Float.atan2 (Float.sqrt (27.0 * d)) (27.0 * J3)
    let c  := 2.0 * Float.sqrt (3.0 * J2)
    let lam := fun k : Float => (I1 + c * Float.cos (phi / 3.0 + twoPiOver3 * k)) / 3.0
    { l₁ := lam 1.0, l₂ := lam 2.0, l₃ := lam 3.0 }

end Eig3x3
