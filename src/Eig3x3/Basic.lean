/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module


/-!
# Eig3x3.Basic — types and the 3×3 vector/matrix vocabulary

Core types (`Vec3`, `Mat3`, `SymmMat3`, `Eigval3`, `Decomposition`), the
public `Float` arithmetic vocabulary, typeclass instances for arithmetic
notation, and scoped unicode notation.

## Usage

Importing this module activates every typeclass instance (`Add`, `Mul`,
`HSMul`, `Inv`, …), so arithmetic on `Vec3` and `Mat3` works immediately
through Lean's core syntax. The Unicode operator glyphs are opt-in:

```
open scoped Eig3x3
```

The single exception is `Pow Float Nat` (the `x ^ 2` spelling on scalars):
it is a `scoped instance`, because an instance on core's `Float` is foreign
territory and must not enter synthesis for consumers who have not opted in.
It activates with the same `open scoped Eig3x3`.

## Notation and Operators

A Mathlib-consistent operator layer for vector/matrix algebra on `Float`,
with no Mathlib dependency. Precedences follow Lean core and Mathlib:

| Notation | Operation              | Declaration   | Convention source               |
|----------|------------------------|---------------|---------------------------------|
| `a + b`  | addition               | `infixl:65`   | Lean core (`HAdd.hAdd`)         |
| `a - b`  | subtraction            | `infixl:65`   | Lean core (`HSub.hSub`)         |
| `-a`     | negation               | prefix        | Lean core (`Neg.neg`)           |
| `a * b`  | multiplication/action  | `infixl:70`   | Lean core (`HMul.hMul`)         |
| `a / s`  | scalar division        | `infixl:70`   | Lean core (`HDiv.hDiv`)         |
| `u ⊗ᵥ v` | vector outer product   | `infixl:70`   | (new; subscript family of `⬝ᵥ`) |
| `u ⬝ᵥ v` | vector dot product     | `infixl:72`   | Mathlib `Matrix.dotProduct`     |
| `A ⬝ₘ B` | matrix Frobenius dot   | `infixl:72`   | (family extension of `⬝ᵥ`)      |
| `s • v`  | scalar multiplication  | `infixr:73`   | Lean core (`HSMul.hSMul`)       |
| `u ⨯₃ v` | cross product          | `infixl:74`   | Mathlib `crossProduct`          |
| `x ^ n`  | power                  | `infixr:80`   | Lean core (`HPow.hPow`)         |
| `A ⊙ B`  | Hadamard product       | `infixl:100`  | Mathlib `Matrix.hadamard`       |
| `\|x\|`  | absolute value         | delimited     | Mathlib `\|a\|` for `abs`       |
| `‖v‖`    | norm                   | delimited     | Mathlib `Norm.norm`             |
| `‖v‖²`   | squared norm           | delimited     | (Mathlib writes `‖x‖ ^ 2`)      |
| `Aᵀ`     | transpose              | `postfix:max` | Mathlib `Matrix.transpose`      |
| `A⁻¹`    | inverse                | `postfix:max` | Lean core (`Inv.inv`)           |

Binding strength: `+ -` (65) < `* /` (70) < `⬝ᵥ ⬝ₘ` (72) < `•` (73) < `⨯₃` (74)
< `^` (80) < `⊙` (100); the postfixes `ᵀ` and `⁻¹` (1024) bind tightest.

Useful consequences:

* `u ⬝ᵥ v ⨯₃ w` parses as `u ⬝ᵥ (v ⨯₃ w)` — the scalar triple product
  needs no parentheses.
* `s • u ⨯₃ v` parses as `s • (u ⨯₃ v)`.
* `u ⬝ᵥ u + v ⬝ᵥ v` parses as `(u ⬝ᵥ u) + (v ⬝ᵥ v)`.
* `M⁻¹ᵀ` parses as `(M⁻¹)ᵀ`.

Watch out: `⬝ᵥ` (72) binds tighter than `*` (70), so `M * v ⬝ᵥ w` parses as
`M * (v ⬝ᵥ w)` — write `(M * v) ⬝ᵥ w` instead.
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
    eigenvectors, matching conventional mathematics and NumPy/SciPy/PyTorch. -/
structure Mat3 where
  c₀ : Vec3
  c₁ : Vec3
  c₂ : Vec3
  deriving Repr


/-- Symmetric 3×3 matrix stored as its six independent entries. -/
structure SymmMat3 where
  a₀₀ : Float
  a₁₁ : Float
  a₂₂ : Float
  a₀₁ : Float
  a₀₂ : Float
  a₁₂ : Float
  deriving Repr


/-- Ordered eigenvalues `l₀ ≤ l₁ ≤ l₂`.

    A distinct type rather than a bare `Vec3`, so ordering guarantees required
    by Eberly's algorithm are carried by the type. Coerces to `Vec3` for
    downstream arithmetic. The coercion is one-way, and ordering guarantees
    do not hold after coercion. -/
structure Eigval3 where
  l₀ : Float
  l₁ : Float
  l₂ : Float
  deriving Repr


/-- Full eigendecomposition: `A = QΛQᵀ = Σᵢ λᵢ cᵢcᵢᵀ`. -/
structure Decomposition where
  /- Eigenvalues in increasing order: `l₀ ≤ l₁ ≤ l₂`. -/
  eigvals : Eigval3
  /- Eigenvector matrix Q; column `cᵢ` is a unit eigenvector for `lᵢ`.
     Right-handed: `Q.det = 1`. -/
  eigvecs : Mat3
  deriving Repr


/-! ## Vec3 operations -/


@[inline] def Vec3.add (u v : Vec3) : Vec3 :=
  ⟨u.x + v.x, u.y + v.y, u.z + v.z⟩


@[inline] def Vec3.neg (v : Vec3) : Vec3 :=
  ⟨-v.x, -v.y, -v.z⟩


@[inline] def Vec3.sub (u v : Vec3) : Vec3 :=
  ⟨u.x - v.x, u.y - v.y, u.z - v.z⟩


@[inline] def Vec3.scale (v : Vec3) (s : Float) : Vec3 :=
  ⟨s * v.x, s * v.y, s * v.z⟩


@[inline] def Vec3.div (v : Vec3) (s : Float) : Vec3 :=
  let invS := 1.0 / s
  ⟨invS * v.x, invS * v.y, invS * v.z⟩


@[inline] def Vec3.dot (u v : Vec3) : Float :=
  u.x * v.x + u.y * v.y + u.z * v.z


@[inline] def Vec3.cross (u v : Vec3) : Vec3 :=
  ⟨u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x⟩


@[inline] def Vec3.normSq (u : Vec3) : Float :=
  u.x * u.x + u.y * u.y + u.z * u.z


@[inline] def Vec3.norm (u : Vec3) : Float :=
  Float.sqrt u.normSq


/-- Outer product of two vectors, producing the rank-1 matrix `u vᵀ`.
    Exposed as `⊗ᵥ`. The columns are `u` scaled by the entries of `v`,
    matching `Mat3`'s column-major layout. -/
@[inline] def Vec3.outer (u v : Vec3) : Mat3 :=
  ⟨u.scale v.x, u.scale v.y, u.scale v.z⟩


/-- Componentwise product, exposed as `⊙` via the `Hadamard` class.
    `Vec3` deliberately has no `Mul` instance: the vector products of this
    library are `⬝ᵥ` and `⨯₃`, and `⊙` is the single entrywise spelling. -/
@[inline] def Vec3.hadamard (u v : Vec3) : Vec3 :=
  ⟨u.x * v.x, u.y * v.y, u.z * v.z⟩


/-- Componentwise map — workhorse for eigenvalue post-processing
    (clipping, shifting, reciprocals). -/
@[inline] def Vec3.map (f : Float → Float) (v : Vec3) : Vec3 :=
  ⟨f v.x, f v.y, f v.z⟩


/-- Componentwise zip. -/
@[inline] def Vec3.map2 (f : Float → Float → Float) (u v : Vec3) : Vec3 :=
  ⟨f u.x v.x, f u.y v.y, f u.z v.z⟩


@[inline] def Vec3.abs (v : Vec3) : Vec3 := v.map Float.abs


/-- Eigenvalues as a plain vector, for downstream arithmetic. -/
@[inline] def Eigval3.toVec3 (e : Eigval3) : Vec3 := ⟨e.l₀, e.l₁, e.l₂⟩


/-- One-way coercion: ordered eigenvalues are a vector, but an arbitrary
    vector is not ordered eigenvalues. -/
instance : Coe Eigval3 Vec3 := ⟨Eigval3.toVec3⟩


/-! ## Mat3 operations -/


/-- Rows become columns: nine field moves, no arithmetic. -/
@[inline] def Mat3.transpose (M : Mat3) : Mat3 :=
  ⟨⟨M.c₀.x, M.c₁.x, M.c₂.x⟩,
   ⟨M.c₀.y, M.c₁.y, M.c₂.y⟩,
   ⟨M.c₀.z, M.c₁.z, M.c₂.z⟩⟩


/-- M·v: linear combination of the columns. -/
@[inline] def Mat3.mulVec (M : Mat3) (v : Vec3) : Vec3 :=
  (M.c₀.scale v.x).add (M.c₁.scale v.y) |>.add (M.c₂.scale v.z)


/-- Mᵀ·v without materializing the transpose: dot with each column. -/
@[inline] def Mat3.transposeMulVec (M : Mat3) (v : Vec3) : Vec3 :=
  ⟨M.c₀.dot v, M.c₁.dot v, M.c₂.dot v⟩


/-- M·N: column `j` of the product is M applied to column `j` of N. -/
@[inline] def Mat3.mul (M N : Mat3) : Mat3 :=
  ⟨M.mulVec N.c₀, M.mulVec N.c₁, M.mulVec N.c₂⟩


@[inline] def Mat3.scale (M : Mat3) (s : Float) : Mat3 :=
  ⟨M.c₀.scale s, M.c₁.scale s, M.c₂.scale s⟩


@[inline] def Mat3.div (M : Mat3) (s : Float) : Mat3 :=
  let invS := 1.0 / s
  ⟨M.c₀.scale invS, M.c₁.scale invS, M.c₂.scale invS⟩


@[inline] def Mat3.add (M N : Mat3) : Mat3 :=
  ⟨M.c₀.add N.c₀, M.c₁.add N.c₁, M.c₂.add N.c₂⟩


@[inline] def Mat3.neg (M : Mat3) : Mat3 :=
  ⟨M.c₀.neg, M.c₁.neg, M.c₂.neg⟩


@[inline] def Mat3.sub (M N : Mat3) : Mat3 :=
  ⟨M.c₀.sub N.c₀, M.c₁.sub N.c₁, M.c₂.sub N.c₂⟩


/-- Componentwise (Hadamard) product, exposed as `⊙`. -/
@[inline] def Mat3.hadamard (M N : Mat3) : Mat3 :=
  ⟨M.c₀.hadamard N.c₀, M.c₁.hadamard N.c₁, M.c₂.hadamard N.c₂⟩


/-- Componentwise absolute value. -/
@[inline] def Mat3.abs (M : Mat3) : Mat3 :=
  ⟨M.c₀.abs, M.c₁.abs, M.c₂.abs⟩


/-- Matrix Frobenius inner product `∑ᵢⱼ Mᵢⱼ Nᵢⱼ = tr(MᵀN)`. -/
@[inline] def Mat3.dot (M N : Mat3) : Float :=
  M.c₀.dot N.c₀ + M.c₁.dot N.c₁ + M.c₂.dot N.c₂


/-- Squared Frobenius norm: sum of squares of all entries. -/
@[inline] def Mat3.normSq (M : Mat3) : Float :=
  M.c₀.normSq + M.c₁.normSq + M.c₂.normSq


/-- Frobenius norm: square root of sum of squares of all entries.
    (Not the operator norm; chosen because it is exactly computable with
    nine multiplications.) -/
@[inline] def Mat3.norm (M : Mat3) : Float :=
  Float.sqrt M.normSq


/-- det M: the triple product `(c₀ × c₁) ⬝ c₂` of the columns. -/
@[inline] def Mat3.det (M : Mat3) : Float :=
  (M.c₀.cross M.c₁).dot M.c₂


/-- tr M: sum of the diagonal entries. -/
@[inline] def Mat3.trace (M : Mat3) : Float :=
  M.c₀.x + M.c₁.y + M.c₂.z


/-- The identity matrix. -/
def Mat3.id : Mat3 :=
  ⟨⟨1.0, 0.0, 0.0⟩, ⟨0.0, 1.0, 0.0⟩, ⟨0.0, 0.0, 1.0⟩⟩


/-- Row-minded constructor: build from three rows. -/
@[inline] def Mat3.ofRows (r₀ r₁ r₂ : Vec3) : Mat3 :=
  (⟨r₀, r₁, r₂⟩ : Mat3).transpose


/-- Invert a 3×3 matrix via cofactor cross products: `M⁻¹ = (1/det M) adj(M)`,
    with `adj(M)` computed from cross products of the rows. Totalized like the
    rest of the `Float` vocabulary: a singular `M` yields `inf`/`NaN` entries
    rather than an error. -/
def Mat3.inv (M : Mat3) : Mat3 :=
  let r₀ : Vec3 := ⟨M.c₀.x, M.c₁.x, M.c₂.x⟩
  let r₁ : Vec3 := ⟨M.c₀.y, M.c₁.y, M.c₂.y⟩
  let r₂ : Vec3 := ⟨M.c₀.z, M.c₁.z, M.c₂.z⟩
  let c₀' := r₁.cross r₂
  let c₁' := r₂.cross r₀
  let c₂' := r₀.cross r₁
  let det := r₀.dot c₀'
  let invDet := 1.0 / det
  ⟨c₀'.scale invDet, c₁'.scale invDet, c₂'.scale invDet⟩


/-- Embed a symmetric matrix. -/
def SymmMat3.toMat3 (A : SymmMat3) : Mat3 :=
  ⟨⟨A.a₀₀, A.a₀₁, A.a₀₂⟩, ⟨A.a₀₁, A.a₁₁, A.a₁₂⟩, ⟨A.a₀₂, A.a₁₂, A.a₂₂⟩⟩


/-- One-way coercion: a symmetric matrix is a matrix, but an arbitrary
    matrix is not symmetric. -/
instance : Coe SymmMat3 Mat3 := ⟨SymmMat3.toMat3⟩


/-! ## Powers on `Float` and `Mat3` -/


/-- Power by repeated multiplication — adequate for the small exponents
(2 and 3) that appear in characteristic-polynomial computations, and cheaper
than a C `pow` call for them. Lean core already provides `Float.pow` with a
`HomogeneousPow Float` instance, so `x ^ (2.0 : Float)` works without this
instance; this one exists for the `x ^ 2` (`Nat` exponent) spelling. -/
def powNat (x : Float) : Nat → Float
  | 0     => 1.0
  | n + 1 => x * powNat x n


/-- Scoped because this is an orphan instance: `Pow` and `Float` both belong
to Lean core, so it must not enter synthesis for consumers who have not opted
in. Activates with `open scoped Eig3x3`. -/
scoped instance : Pow Float Nat where
  pow := powNat


/-- Matrix power by repeated multiplication, with `M ^ 0 = Mat3.id`. -/
def Mat3.powNat (M : Mat3) : Nat → Mat3
  | 0     => Mat3.id
  | n + 1 => M.mul (Mat3.powNat M n)


/-- Global, unlike `Pow Float Nat`: `Mat3` is this package's own type, so no
downstream collision is possible. -/
instance : Pow Mat3 Nat where
  pow := Mat3.powNat


/-! ## Custom typeclasses behind the delimited/Hadamard notation -/


/-- Absolute value behind a typeclass, so the same bars work on
`Float` (magnitude), `Vec3` (componentwise), and `Mat3` (componentwise).
The class is this package's own, so instances on core types are safe. -/
class Abs (α : Type) where
  abs : α → α


instance : Abs Float where
  abs := Float.abs


instance : Abs Vec3 where
  abs := Vec3.abs


instance : Abs Mat3 where
  abs := Mat3.abs


/-- Euclidean norm for vectors, Frobenius norm for matrices. -/
class Norm (α : Type) where
  norm : α → Float


/-- Squared norm: Euclidean for vectors, Frobenius for matrices. -/
class NormSq (α : Type) where
  normSq : α → Float


instance : Norm Vec3 where
  norm := Vec3.norm


instance : NormSq Vec3 where
  normSq := Vec3.normSq


instance : Norm Mat3 where
  norm := Mat3.norm


instance : NormSq Mat3 where
  normSq := Mat3.normSq


/-- Hadamard (componentwise) product. -/
class Hadamard (α : Type) where
  hadamard : α → α → α


instance : Hadamard Vec3 where
  hadamard := Vec3.hadamard


instance : Hadamard Mat3 where
  hadamard := Mat3.hadamard


/-! ## Arithmetic instances (global: these are this package's own types) -/


instance : Add Vec3 := ⟨Vec3.add⟩
instance : Neg Vec3 := ⟨Vec3.neg⟩
instance : Sub Vec3 := ⟨Vec3.sub⟩
-- Deliberately no `Mul Vec3`: the entrywise product is `⊙` (see `Hadamard`),
-- and `*` on this package's types always means a genuine product or action.
instance : HSMul Float Vec3 Vec3 := ⟨fun s v => v.scale s⟩
instance : HDiv Vec3 Float Vec3 := ⟨Vec3.div⟩


instance : Add Mat3 := ⟨Mat3.add⟩
instance : Neg Mat3 := ⟨Mat3.neg⟩
instance : Sub Mat3 := ⟨Mat3.sub⟩
instance : Mul Mat3 := ⟨Mat3.mul⟩
instance : HMul Mat3 Vec3 Vec3 := ⟨Mat3.mulVec⟩
/-- Row-vector times matrix product `vᵀ M` without materializing the transpose. -/
instance : HMul Vec3 Mat3 Vec3 := ⟨fun v M => M.transposeMulVec v⟩
instance : HSMul Float Mat3 Mat3 := ⟨fun s M => M.scale s⟩
instance : HDiv Mat3 Float Mat3 := ⟨Mat3.div⟩


/-- Matrix inverse through Lean core's `Inv` class, so `M⁻¹` works with no new
syntax at all: core defines `postfix:max "⁻¹" => Inv.inv`. -/
instance : Inv Mat3 := ⟨Mat3.inv⟩


/-! ## Scoped operators (`open scoped Eig3x3` to activate) -/


/-- Outer product producing the rank-1 matrix `u vᵀ`, at multiplicative
precedence 70. Type `\otimes` then `\_v`. -/
scoped infixl:70 " ⊗ᵥ " => Vec3.outer


/-- Vector dot product at precedence 72, matching Mathlib's
`Matrix.dotProduct`. Type `\cdot` then `\_v` (or `\dot\_v`). -/
scoped infixl:72 " ⬝ᵥ " => Vec3.dot


/-- Matrix Frobenius inner product at precedence 72. Type `\cdot\_m`. -/
scoped infixl:72 " ⬝ₘ " => Mat3.dot


/-- Cross product at precedence 74, matching Mathlib's `crossProduct`.
`⨯₃` (U+2A2F) is Mathlib's current glyph; `×₃` (U+00D7, typed `\times\_3`) is
the legacy alias, kept because it is substantially easier to type. -/
scoped infixl:74 " ⨯₃ " => Vec3.cross
scoped infixl:74 " ×₃ " => Vec3.cross


/-- Absolute value: magnitude on `Float`, entrywise on `Vec3`/`Mat3`.
Mirrors Mathlib's `|a|` notation for `abs`. -/
scoped notation:max "|" x "|" => Abs.abs x


/-- Norm: Euclidean on `Vec3`, Frobenius on `Mat3`. Type `\Vert`. -/
scoped notation:max "‖" v "‖" => Norm.norm v


/-- Squared norm. A bare postfix `x²` is not possible in Lean (`²` is an
identifier-continuation character, so `x²` lexes as a single name). `‖²` works
because after the closing `‖` the lexer is reading a symbol token and absorbs
the superscript — the same mechanism that makes `⨯₃` lexable. Write it tight:
`‖v‖²`, no space before the `²`. -/
scoped notation:max "‖" v "‖²" => NormSq.normSq v


/-- Hadamard (componentwise) product at precedence 100, matching Mathlib's
`Matrix.hadamard`. Type `\odot`. -/
scoped infixl:100 " ⊙ " => Hadamard.hadamard


/-- Transpose at `postfix:max` (= 1024), exactly Mathlib's `Matrix.transpose`
convention. Type `\^T`. -/
scoped postfix:max "ᵀ" => Mat3.transpose


end  -- public section


/-! ## Internal pipeline helpers (package-private) -/


/-- Simple inline helper for computing the max of 3 values. -/
@[inline] def max3 {α : Type} [Max α] (a b c : α) : α :=
  max a (max b c)


/-- Apply the symmetric matrix to a vector. -/
@[inline] def SymmMat3.mulVec (A : SymmMat3) (v : Vec3) : Vec3 :=
  ⟨A.a₀₀ * v.x + A.a₀₁ * v.y + A.a₀₂ * v.z,
   A.a₀₁ * v.x + A.a₁₁ * v.y + A.a₁₂ * v.z,
   A.a₀₂ * v.x + A.a₁₂ * v.y + A.a₂₂ * v.z⟩


@[inline] def SymmMat3.scale (A : SymmMat3) (s : Float) : SymmMat3 :=
  ⟨s * A.a₀₀, s * A.a₁₁, s * A.a₂₂, s * A.a₀₁, s * A.a₀₂, s * A.a₁₂⟩


/-- Largest |entry|, used for preconditioning (Eberly's overflow guard). -/
def SymmMat3.maxAbsEntry (A : SymmMat3) : Float :=
  let m₀ := max A.a₀₀.abs A.a₀₁.abs
  let m₁ := max A.a₀₂.abs A.a₁₁.abs
  let m₂ := max A.a₁₂.abs A.a₂₂.abs
  max3 m₀ m₁ m₂


/-- Scale eigenvalues by `s`. Requires `s > 0` to preserve ordering. -/
@[inline] def Eigval3.scale (e : Eigval3) (s : Float) : Eigval3 :=
  ⟨s * e.l₀, s * e.l₁, s * e.l₂⟩


end Eig3x3
