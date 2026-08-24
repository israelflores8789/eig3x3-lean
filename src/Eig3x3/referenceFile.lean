/-!
# Eig3x3 — closed-form eigendecomposition of real symmetric 3×3 matrices

Pure Lean 4 (`Float` only): no Mathlib, no FFI, no Python.

## Provenance

* **Eigenvalues** — Habera & Zilian, "Numerically stable evaluation of
  closed-form expressions for eigenvalues of 3×3 matrices",
  arXiv:2511.00292v2 (2025). Invariants I₁ (Alg. 1), J₂ (Alg. 2), J₃ (Alg. 5),
  quadrant-safe angle φ = atan2(√(27Δ), 27J₃) (Eq. 4), ordered eigenvalues
  λ₁ ≤ λ₂ ≤ λ₃ (Eq. 2). Reference C implementation: `eig3x3` (MIT license).

* **Discriminant** — Habera & Zilian 2025, Algorithm 8: the Cauchy–Binet
  sum-of-products form Δ = Σᵢ wᵢ uᵢ vᵢ of Parlett's discriminant
  factorization (formula originating in Habera–Zilian 2021,
  arXiv:2111.02117), specialized to symmetric A where u = v = DX(A) and the
  sum becomes a sum of squares. The auxiliary terms contain the diagonal
  entries of A only through their differences — the deviatoric
  backward-stability property that keeps the evaluation accurate near J₂ = 0 —
  while the sum-of-squares form eliminates catastrophic cancellation near
  Δ = 0 (near-multiple eigenvalues).

* **Eigenvectors** — D. Eberly, "A Robust Eigensolver for 3×3 Symmetric
  Matrices", Geometric Tools (documentation: CC BY 4.0), specifically the
  non-iterative `NISymmetricEigensolver3x3` in `SymmetricEigensolver3x3.h`
  (code: Boost Software License 1.0). Max-abs preconditioning, isolated-first
  cross products (`ComputeEigenvector0`), robust orthogonal complement
  (`ComputeOrthogonalComplement`), 2×2 reduction in the complement
  (`ComputeEigenvector1`), right-handed completion.

## Deviations from the sources

* Eberly's acos-based eigenvalue evaluation is replaced by the Habera–Zilian
  invariant pipeline (H–Z prove the acos form is unstable near repeated
  eigenvalues; their arctan form is not).
* Eberly's `sign(halfDet)` branch is replaced by a direct gap comparison on the
  ordered eigenvalues — equivalent in exact arithmetic, but more direct.
* No sorting stage: H–Z returns λ₁ ≤ λ₂ ≤ λ₃, and the eigenvector assembly is
  right-handed by construction (the third vector is always a cross product).
* `eigvecIsolated` adds a defensive exact-zero fallback; Eberly relies on
  exact-arithmetic rank 2 and has no such guard.
* `deltaNaive` (4J₂³ − 27J₃²) is retained solely for benchmarking, to
  reproduce the paper's naive-vs-present comparisons.

## Conditioning scope note

H–Z report that Algorithm 8 exceeds its (lowest-order) forward-stability bound
for benchmarks with an *ill-conditioned eigenbasis* (their Figs. 5e/5f, the U₂
congruences). That regime requires a non-orthogonal eigenvector matrix, which
only non-symmetric matrices possess; symmetric matrices always have an
orthogonal eigenbasis (κ₂ = 1), the regime in which Algorithm 8 meets the
bound (their Figs. 5a–5d). This library is scoped to real symmetric matrices,
so the failing regime is out of scope by construction.

## Validation (op-for-op float64 mirror vs `numpy.linalg.eigvalsh`)

  * 20k random symmetric matrices:        max |Δλ| ≈ 8e-15
  * near-double path diag(−1,1,1+δ):      ≈ 6e-16 for δ ∈ [1e-16, 1e-4]
  * near-triple path diag(1,1,1+δ):       ≈ 2e-16 for δ ≤ 1e-8
  * double eigenvalue at small scale:     ≈ 1e-19 (boundary adversarial cases)
  * exact scaled identity / diagonal:     exact

## Scope

Real symmetric matrices (i.e. Hermitian over ℝ). Complex Hermitian is out of
scope (cf. Kopp 2008 for that case). Strictly eigendecomposition plus
certificates: no clipping, no projections, no downstream numerics.
-/

namespace Eig3x3

/-! ## Types -/

/-- Symmetric 3×3 matrix stored as its six independent entries. -/
structure SymmMat3 where
  a00 : Float
  a11 : Float
  a22 : Float
  a01 : Float
  a02 : Float
  a12 : Float
  deriving Repr

/-- Vector in ℝ³. -/
structure Vec3 where
  x : Float
  y : Float
  z : Float
  deriving Repr

/-- Ordered eigenvalues `l1 ≤ l2 ≤ l3`. -/
structure Eigen3 where
  l1 : Float
  l2 : Float
  l3 : Float
  deriving Repr

/-- Full eigendecomposition: `A = Σᵢ λᵢ vᵢvᵢᵀ` with `{v1, v2, v3}` a
    right-handed orthonormal set and `vᵢ` a unit eigenvector for `λᵢ`. -/
structure Decomposition where
  evals : Eigen3
  v1 : Vec3
  v2 : Vec3
  v3 : Vec3
  deriving Repr

/-! ## Vec3 and SymmMat3 operations -/

@[inline] def Vec3.cross (u v : Vec3) : Vec3 :=
  ⟨u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩

@[inline] def Vec3.dot (u v : Vec3) : Float :=
  u.x * v.x + u.y * v.y + u.z * v.z

@[inline] def Vec3.normSq (u : Vec3) : Float :=
  u.x * u.x + u.y * u.y + u.z * u.z

@[inline] def Vec3.scale (v : Vec3) (s : Float) : Vec3 :=
  ⟨s * v.x, s * v.y, s * v.z⟩

@[inline] def Vec3.sub (u v : Vec3) : Vec3 :=
  ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩

@[inline] def SymmMat3.mulVec (A : SymmMat3) (v : Vec3) : Vec3 :=
  ⟨A.a00 * v.x + A.a01 * v.y + A.a02 * v.z,
   A.a01 * v.x + A.a11 * v.y + A.a12 * v.z,
   A.a02 * v.x + A.a12 * v.y + A.a22 * v.z⟩

@[inline] def SymmMat3.scale (A : SymmMat3) (s : Float) : SymmMat3 :=
  ⟨s * A.a00, s * A.a11, s * A.a22, s * A.a01, s * A.a02, s * A.a12⟩

/-- Largest |entry|, used for preconditioning (Eberly's overflow guard). -/
def SymmMat3.maxAbsEntry (A : SymmMat3) : Float :=
  let m1 := Float.max A.a00.abs A.a01.abs
  let m2 := Float.max A.a02.abs A.a11.abs
  let m3 := Float.max A.a12.abs A.a22.abs
  Float.max m1 (Float.max m2 m3)

/-- Scale eigenvalues by `s`. Requires `s > 0` to preserve ordering. -/
@[inline] def Eigen3.scale (e : Eigen3) (s : Float) : Eigen3 :=
  ⟨s * e.l1, s * e.l2, s * e.l3⟩

/-! ## Eigenvalues (Habera–Zilian 2025) -/

/-- Algorithm 1: I₁ = tr(A). -/
def i1 (A : SymmMat3) : Float := A.a00 + A.a11 + A.a22

/-- Algorithm 2, symmetric case: J₂ = ½ tr(dev A)² from diagonal differences
    and off-diagonal squares. Exactly zero for scaled identities (H–Z Eq. 62),
    which is what keeps the near-triple-eigenvalue case stable. -/
def j2 (A : SymmMat3) : Float :=
  let d0 := A.a00 - A.a11
  let d1 := A.a00 - A.a22
  let d2 := A.a11 - A.a22
  (d0 * d0 + d1 * d1 + d2 * d2) / 6.0
    + (A.a01 * A.a01 + A.a02 * A.a02 + A.a12 * A.a12)

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

/-- Naive discriminant Δ = 4J₂³ − 27J₃², clamped to [0, ∞). Kept for
    benchmarking only (reproduces the paper's naive-vs-present comparisons);
    suffers catastrophic cancellation near double eigenvalues with finite J₂
    (observed ≈5e-9 absolute eigenvalue error on the D2 path). -/
def deltaNaive (J2 J3 : Float) : Float :=
  let d := 4.0 * J2 * J2 * J2 - 27.0 * J3 * J3
  if d < 0.0 then 0.0 else d

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
    2021 paper (arXiv:2111.02117, Eq. 29). Double-check this sign against the
    paper PDF before release.

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
    eigenvalues come out as exactly I₁/3. -/
def eigvals (A : SymmMat3) : Eigen3 :=
  let I1 := i1 A
  let J2 := j2 A
  let J3 := j3 A
  let d  := delta A
  let phi := Float.atan2 (Float.sqrt (27.0 * d)) (27.0 * J3)
  let c  := 2.0 * Float.sqrt (3.0 * J2)
  let lam := fun k : Float => (I1 + c * Float.cos (phi / 3.0 + twoPiOver3 * k)) / 3.0
  { l1 := lam 1.0, l2 := lam 2.0, l3 := lam 3.0 }

/-! ## Eigenvectors (Eberly) -/

/-- Eigenvector for the *isolated* eigenvalue λ: the largest cross product of
    rows of (A − λI). GTE `ComputeEigenvector0`.

    Defensive addition: if all crosses are exactly zero (A = λI up to rounding,
    i.e. a triple eigenvalue), returns e₁ — any unit vector is an eigenvector
    then, and the residual certificate stays tiny because ‖A − λI‖ does. -/
def eigvecIsolated (A : SymmMat3) (lam : Float) : Vec3 :=
  let r0 : Vec3 := ⟨A.a00 - lam, A.a01, A.a02⟩
  let r1 : Vec3 := ⟨A.a01, A.a11 - lam, A.a12⟩
  let r2 : Vec3 := ⟨A.a02, A.a12, A.a22 - lam⟩
  let c0 := r1.cross r2
  let c1 := r0.cross r2
  let c2 := r0.cross r1
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

/-- Robustly compute U, V so that {U, V, w} is a right-handed orthonormal set.
    Requires `w` unit-length. GTE `ComputeOrthogonalComplement`. -/
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
    of an adjacent, well-separated eigenvalue. Restricts (A − λI) to the plane
    via J = [U V], then solves the 2×2 null system M X = 0 by largest-row
    selection with division-free normalization. If M ≈ 0 (repeated eigenvalue),
    any vector in the plane works, so U is returned.
    GTE `ComputeEigenvector1`. -/
def eigvecInPlane (A : SymmMat3) (v0 : Vec3) (lam : Float) : Vec3 :=
  let (u, v) := orthonormalComplement v0
  let au := A.mulVec u
  let av := A.mulVec v
  let m00 := u.dot au - lam
  let m01 := u.dot av
  let m11 := v.dot av - lam
  if m00.abs < m11.abs then
    -- Solve using row 1: m01·x0 + m11·x1 = 0
    let maxAbs := Float.max m11.abs m01.abs
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
    let maxAbs := Float.max m00.abs m01.abs
    if maxAbs == 0.0 then u
    else if m00.abs < m01.abs then
      let t := m00 / m01
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale n).sub (v.scale (t * n))        -- X = (1, −t)·n
    else
      let t := m01 / m00
      let n := 1.0 / (1.0 + t * t).sqrt
      (u.scale (t * n)).sub (v.scale n)        -- X = (t, −1)·n

/-! ## Assembly -/

/-- Full eigendecomposition of a real symmetric 3×3 matrix.

    Pipeline: zero fast path → max-abs preconditioning (Eberly) → ordered
    eigenvalues on the scaled matrix (Habera–Zilian) → isolated-eigenvalue-first
    eigenvector assembly (Eberly), where the gap comparison on the ordered
    eigenvalues replaces Eberly's `sign(halfDet)` branch → right-handed
    completion by cross product → eigenvalue rescaling. -/
def eigendecomp (A : SymmMat3) : Decomposition :=
  let maxAbs := A.maxAbsEntry
  if maxAbs == 0.0 then
    { evals := ⟨0.0, 0.0, 0.0⟩,
      v1 := ⟨1.0, 0.0, 0.0⟩,
      v2 := ⟨0.0, 1.0, 0.0⟩,
      v3 := ⟨0.0, 0.0, 1.0⟩ }
  else
    let s := 1.0 / maxAbs
    let B := A.scale s
    let e := eigvals B
    let gapLo := e.l2 - e.l1
    let gapHi := e.l3 - e.l2
    let (v1, v2, v3) :=
      if gapLo < gapHi then
        -- λ₃ is isolated (Eberly's halfDet ≥ 0 case)
        let v3 := eigvecIsolated B e.l3
        let v2 := eigvecInPlane B v3 e.l2
        (v2.cross v3, v2, v3)
      else
        -- λ₁ is isolated (Eberly's halfDet < 0 case)
        let v1 := eigvecIsolated B e.l1
        let v2 := eigvecInPlane B v1 e.l2
        (v1, v2, v1.cross v2)
    { evals := e.scale maxAbs, v1 := v1, v2 := v2, v3 := v3 }

/-! ## Certificates -/

/-- Residual ‖Av − λv‖∞ for a claimed eigenpair. -/
def residual (A : SymmMat3) (lam : Float) (v : Vec3) : Float :=
  let av := A.mulVec v
  let r : Vec3 := ⟨av.x - lam * v.x, av.y - lam * v.y, av.z - lam * v.z⟩
  Float.max r.x.abs (Float.max r.y.abs r.z.abs)

/-- Orthonormality certificate: max |vᵢ·vⱼ − δᵢⱼ| over all pairs. -/
def orthogonalityError (v1 v2 v3 : Vec3) : Float :=
  let m1 := Float.max (v1.dot v1 - 1.0).abs (v2.dot v2 - 1.0).abs
  let m2 := Float.max (v3.dot v3 - 1.0).abs (v1.dot v2).abs
  let m3 := Float.max (v1.dot v3).abs (v2.dot v3).abs
  Float.max m1 (Float.max m2 m3)

/-- Reconstruct the matrix from a decomposition: A = Σᵢ λᵢ vᵢvᵢᵀ. -/
def reconstruct (d : Decomposition) : SymmMat3 :=
  let comp (m : Float) (v : Vec3) : SymmMat3 :=
    ⟨m * v.x * v.x, m * v.y * v.y, m * v.z * v.z,
     m * v.x * v.y, m * v.x * v.z, m * v.y * v.z⟩
  let b1 := comp d.evals.l1 d.v1
  let b2 := comp d.evals.l2 d.v2
  let b3 := comp d.evals.l3 d.v3
  ⟨b1.a00 + b2.a00 + b3.a00,
   b1.a11 + b2.a11 + b3.a11,
   b1.a22 + b2.a22 + b3.a22,
   b1.a01 + b2.a01 + b3.a01,
   b1.a02 + b2.a02 + b3.a02,
   b1.a12 + b2.a12 + b3.a12⟩

/-- Reconstruction certificate: max |entry| of QΛQᵀ − A. -/
def reconstructionError (A : SymmMat3) (d : Decomposition) : Float :=
  let B := reconstruct d
  let m1 := Float.max (B.a00 - A.a00).abs (B.a11 - A.a11).abs
  let m2 := Float.max (B.a22 - A.a22).abs (B.a01 - A.a01).abs
  let m3 := Float.max (B.a02 - A.a02).abs (B.a12 - A.a12).abs
  Float.max m1 (Float.max m2 m3)

/-- All three certificates for a decomposition of `A`, bundled. -/
structure Certificates where
  maxResidual : Float
  orthogonality : Float
  reconstruction : Float
  deriving Repr

def certify (A : SymmMat3) (d : Decomposition) : Certificates :=
  let r1 := residual A d.evals.l1 d.v1
  let r2 := residual A d.evals.l2 d.v2
  let r3 := residual A d.evals.l3 d.v3
  { maxResidual := Float.max r1 (Float.max r2 r3),
    orthogonality := orthogonalityError d.v1 d.v2 d.v3,
    reconstruction := reconstructionError A d }

/-! ## Smoke tests -/

section SmokeTests

/-- [[2,1,0],[1,2,1],[0,1,2]] — exact eigenvalues 2−√2, 2, 2+√2. -/
def workedExample : SymmMat3 := ⟨2.0, 2.0, 2.0, 1.0, 0.0, 1.0⟩

-- Expect l1 ≈ 0.5857864376, l2 = 2.0, l3 ≈ 3.4142135624
#eval (eigendecomp workedExample).evals

-- Expect all three certificates ≈ 1e-16 or smaller
#eval certify workedExample (eigendecomp workedExample)

-- Zero matrix: exact fast path, identity basis
#eval eigendecomp ⟨0.0, 0.0, 0.0, 0.0, 0.0, 0.0⟩

-- Scaled identity: exact triple eigenvalue 3.7, exercises the dmax == 0 fallback
#eval eigendecomp ⟨3.7, 3.7, 3.7, 0.0, 0.0, 0.0⟩

-- Diagonal with distinct entries: expect evals (1, 2, 3) and axis eigenvectors
#eval eigendecomp ⟨1.0, 2.0, 3.0, 0.0, 0.0, 0.0⟩

/-- Near-double eigenvalue: eigenvalues ≈ (−1, 1−1e-8, 1+1e-8).
    Exercises the clustered-top gap branch; with the Algorithm 8 discriminant
    this achieves machine precision (the naive Δ gave ≈2.5e-9 here). -/
def nearDouble : SymmMat3 := ⟨1.0, 1.0, -1.0, 1.0e-8, 0.0, 0.0⟩

-- Expect l1 ≈ -1.0, l2 ≈ 0.99999999, l3 ≈ 1.00000001
#eval (eigendecomp nearDouble).evals

-- Expect certificates ≈ 1e-16
#eval certify nearDouble (eigendecomp nearDouble)

/-- The matrix that caught the r₁₀ transcription error during porting:
    exact Δ = 13,021,520. Kept as a regression test. -/
def regressionMatrix : SymmMat3 := ⟨-2.0, 2.0, -9.0, 9.0, 8.0, 6.0⟩

-- Expect delta = 13021520.0 (exactly representable)
#eval delta regressionMatrix

end SmokeTests

end Eig3x3
