module

/-!
# Eig3x3.Basic — types and the 3×3 vector/matrix vocabulary

Core types (`Vec3`, `Mat3`, `SymmMat3`, `Eigval3`, `Decomposition`), the
public `Float` arithmetic vocabulary, typeclass instances for arithmetic
notation, and scoped unicode notation.
-/

namespace Eig3x3

public section

/-- Vector over `Float` in real (ℝ³) space. -/
structure Vec3 where
  x : Float
  y : Float
  z : Float
  deriving Repr

/-- General 3×3 matrix over `Float`, stored by columns. The primary matrix
    this library produces is the eigenvector matrix Q, whose columns are
    eigenvectors matches conventional mathematics and NumPy/SciPy/PyTorch. -/
structure Mat3 where
  c₁ : Vec3
  c₂ : Vec3
  c₃ : Vec3
  deriving Repr

/-- Symmetric 3×3 matrix stored as its six independent entries. -/
structure SymmMat3 where
  a00 : Float
  a11 : Float
  a22 : Float
  a01 : Float
  a02 : Float
  a12 : Float
  deriving Repr

/-- Ordered eigenvalues `l₁ ≤ l₂ ≤ l₃`.

    A distinct type rather than a bare `Vec3`, so ordering guarantees required
    by Eberly's algorithm are carried by the type. Coerces to `Vec3` for
    downstream arithmetic. The coercion is one-way, and ordering guarantees
    do not hold after coercion. -/
structure Eigval3 where
  l₁ : Float
  l₂ : Float
  l₃ : Float
  deriving Repr

/-- Full eigendecomposition: `A = QΛQᵀ = Σᵢ λᵢ cᵢcᵢᵀ`.

    `eigvals`:
      Eigenvalues in increasing order: `l₁ ≤ l₂ ≤ l₃`.

    `eigvecs`:
      Eigenvector matrix Q; column `cᵢ` is a unit eigenvector for `lᵢ`.
      Right-handed: `Q.det = 1`.
-/
structure Decomposition where
  eigvals : Eigval3
  eigvecs : Mat3
  deriving Repr

/-! ## Vec3 operations -/

@[inline] def Vec3.add (u v : Vec3) : Vec3 :=
  ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩

@[inline] def Vec3.sub (u v : Vec3) : Vec3 :=
  ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩

@[inline] def Vec3.scale (v : Vec3) (s : Float) : Vec3 :=
  ⟨s * v.x, s * v.y, s * v.z⟩

@[inline] def Vec3.dot (u v : Vec3) : Float :=
  u.x * v.x + u.y * v.y + u.z * v.z

@[inline] def Vec3.cross (u v : Vec3) : Vec3 :=
  ⟨u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩

@[inline] def Vec3.normSq (u : Vec3) : Float :=
  u.x * u.x + u.y * u.y + u.z * u.z

/-- Componentwise map — the workhorse for eigenvalue post-processing
    (clipping, shifting, reciprocals). -/
@[inline] def Vec3.map (f : Float → Float) (v : Vec3) : Vec3 :=
  ⟨f v.x, f v.y, f v.z⟩

/-- Componentwise zip. -/
@[inline] def Vec3.map2 (f : Float → Float → Float) (u v : Vec3) : Vec3 :=
  ⟨f u.x v.x, f u.y v.y, f u.z v.z⟩

@[inline] def Vec3.abs (v : Vec3) : Vec3 := v.map Float.abs

/-! ## Mat3 operations -/

/-- Rows become columns: nine field moves, no arithmetic. -/
@[inline] def Mat3.transpose (M : Mat3) : Mat3 :=
  ⟨⟨M.c₁.x, M.c₂.x, M.c₃.x⟩,
   ⟨M.c₁.y, M.c₂.y, M.c₃.y⟩,
   ⟨M.c₁.z, M.c₂.z, M.c₃.z⟩⟩

/-- M·v: linear combination of the columns. -/
@[inline] def Mat3.mulVec (M : Mat3) (v : Vec3) : Vec3 :=
  (M.c₁.scale v.x).add (M.c₂.scale v.y) |>.add (M.c₃.scale v.z)

/-- Mᵀ·v without materializing the transpose: dot with each column. -/
@[inline] def Mat3.transposeMulVec (M : Mat3) (v : Vec3) : Vec3 :=
  ⟨M.c₁.dot v, M.c₂.dot v, M.c₃.dot v⟩

/-- M·N: column `j` of the product is M applied to column `j` of N. -/
@[inline] def Mat3.mul (M N : Mat3) : Mat3 :=
  ⟨M.mulVec N.c₁, M.mulVec N.c₂, M.mulVec N.c₃⟩

@[inline] def Mat3.scale (M : Mat3) (s : Float) : Mat3 :=
  ⟨M.c₁.scale s, M.c₂.scale s, M.c₃.scale s⟩

/-- det M: the triple product `(c₁ × c₂) ⬝ c₃` of the columns. -/
@[inline] def Mat3.det (M : Mat3) : Float :=
  (M.c₁.cross M.c₂).dot M.c₃

@[inline] def Mat3.trace (M : Mat3) : Float :=
  M.c₁.x + M.c₂.y + M.c₃.z

/-- The identity matrix. -/
def Mat3.id : Mat3 :=
  ⟨⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨0.0, 0.0, 1.0⟩⟩

/-- Row-minded constructor: build from three rows. -/
@[inline] def Mat3.ofRows (r1 r2 r3 : Vec3) : Mat3 :=
  (⟨r1, r2, r3⟩ : Mat3).transpose

/-- Embed a symmetric matrix. -/
def SymmMat3.toMat3 (A : SymmMat3) : Mat3 :=
  ⟨⟨A.a00, A.a01, A.a02⟩, ⟨A.a01, A.a11, A.a12⟩, ⟨A.a02, A.a12, A.a22⟩⟩

/-- Eigenvalues as a plain vector, for downstream arithmetic. -/
@[inline] def Eigval3.toVec3 (e : Eigval3) : Vec3 := ⟨e.l₁, e.l₂, e.l₃⟩

/-- One-way coercion: ordered eigenvalues are a vector, but an arbitrary
    vector is not ordered eigenvalues. -/
instance : Coe Eigval3 Vec3 := ⟨Eigval3.toVec3⟩

/-! ## Arithmetic instances -/

instance : Add Vec3 := ⟨Vec3.add⟩
instance : Sub Vec3 := ⟨Vec3.sub⟩
instance : HSMul Float Vec3 Vec3 := ⟨fun s v => v.scale s⟩
instance : Mul Mat3 := ⟨Mat3.mul⟩
instance : HMul Mat3 Vec3 Vec3 := ⟨Mat3.mulVec⟩
instance : HSMul Float Mat3 Mat3 := ⟨fun s M => M.scale s⟩

end  -- public section

/-- Transpose notation, matching Mathlib's `Matrix.transpose` and postfix
    precedence convention. Activate with `open scoped Eig3x3`. -/
public scoped postfix:max "ᵀ" => Mat3.transpose

/-- Dot-product notation at Mathlib's infixl precedence.
    Activate with `open scoped Eig3x3`. -/
public scoped infixl:72 " ⬝ " => Vec3.dot

/-! ## Internal pipeline helpers (package-private) -/

/-- Simple inline helper for computing the max of 3 values. -/
@[inline] def max3 {α : Type} [Max α] (a b c : α) : α :=
  max a (max b c)

/-- Apply the symmetric matrix to a vector. -/
@[inline] def SymmMat3.mulVec (A : SymmMat3) (v : Vec3) : Vec3 :=
  ⟨A.a00 * v.x + A.a01 * v.y + A.a02 * v.z,
   A.a01 * v.x + A.a11 * v.y + A.a12 * v.z,
   A.a02 * v.x + A.a12 * v.y + A.a22 * v.z⟩

@[inline] def SymmMat3.scale (A : SymmMat3) (s : Float) : SymmMat3 :=
  ⟨s * A.a00, s * A.a11, s * A.a22, s * A.a01, s * A.a02, s * A.a12⟩

/-- Largest |entry|, used for preconditioning (Eberly's overflow guard). -/
def SymmMat3.maxAbsEntry (A : SymmMat3) : Float :=
  let m1 := max A.a00.abs A.a01.abs
  let m2 := max A.a02.abs A.a11.abs
  let m3 := max A.a12.abs A.a22.abs
  max3 m1 m2 m3

/-- Scale eigenvalues by `s`. Requires `s > 0` to preserve ordering. -/
@[inline] def Eigval3.scale (e : Eigval3) (s : Float) : Eigval3 :=
  ⟨s * e.l₁, s * e.l₂, s * e.l₃⟩

end Eig3x3
