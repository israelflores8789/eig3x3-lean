/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3
import Tests.Util

/-!
# Tests.KnownAnswer — fixed inputs with known outputs

Exact-in-float64 checks for the vector/matrix vocabulary (small integers
only), known-answer eigenvalue checks, the exact fast paths (zero matrix,
scaled identity), and end-to-end exercise of the notation (`ᵀ`, `⁻¹`, `⬝ᵥ`,
`⬝ₘ`, `⊗ᵥ`, `⨯₃`, `•`, `⊙`, `|·|`, `‖·‖`, `‖·‖²`, `^ⁿ`, `^`, `*`) and the
`Eigval3 → Vec3` and `SymmMat3 → Mat3` coercions.
-/

namespace Eig3x3.Tests

open scoped Eig3x3

public def runKnownAnswer : IO Unit := do
  let u : Vec3 := ⟨1.0, 2.0, 3.0⟩
  let v : Vec3 := ⟨4.0, 5.0, 6.0⟩
  let w : Vec3 := ⟨2.0, -1.0, 1.0⟩
  let s : Float := 2.5

  -- 1. Vec3 arithmetic & operators (all exact)
  assertTrue "vec add"
    ((u + v).approx ⟨5.0, 7.0, 9.0⟩ 0.0)
  assertTrue "vec sub"
    ((v - u).approx ⟨3.0, 3.0, 3.0⟩ 0.0)
  assertTrue "vec neg"
    ((-u).approx ⟨-1.0, -2.0, -3.0⟩ 0.0)
  assertTrue "vec div"
    (((⟨2.0, 4.0, 6.0⟩ : Vec3) / 2.0).approx ⟨1.0, 2.0, 3.0⟩ 0.0)
  assertTrue "vec smul"
    ((2.0 • u).approx ⟨2.0, 4.0, 6.0⟩ 0.0)
  assertClose "vec dot ⬝ᵥ" (u ⬝ᵥ v) 32.0 0.0
  assertTrue "vec cross ⨯₃"
    (((⟨1.0, 0.0, 0.0⟩ : Vec3) ⨯₃ ⟨0.0, 1.0, 0.0⟩).approx ⟨0.0, 0.0, 1.0⟩ 0.0)
  assertTrue "vec cross self is zero"
    ((u ⨯₃ u).approx ⟨0.0, 0.0, 0.0⟩ 0.0)
  assertTrue "vec outer ⊗ᵥ"
    ((u ⊗ᵥ v).approx ⟨⟨4.0, 8.0, 12.0⟩, ⟨5.0, 10.0, 15.0⟩, ⟨6.0, 12.0, 18.0⟩⟩ 0.0)
  assertTrue "vec hadamard ⊙"
    ((u ⊙ v).approx ⟨4.0, 10.0, 18.0⟩ 0.0)
  assertClose "vec normSq ‖v‖²" ‖u‖² 14.0 0.0
  assertClose "vec norm ‖v‖" ‖(⟨3.0, 4.0, 0.0⟩ : Vec3)‖ 5.0 0.0
  assertTrue "vec abs |v|"
    (|(⟨-1.0, 2.0, -3.0⟩ : Vec3)|.approx ⟨1.0, 2.0, 3.0⟩ 0.0)
  assertTrue "vec map"
    (((⟨-1.0, 2.0, -3.0⟩ : Vec3).map Float.abs).approx ⟨1.0, 2.0, 3.0⟩ 0.0)
  assertTrue "vec map2"
    ((Vec3.map2 (· * ·) ⟨1.0, 2.0, 3.0⟩ ⟨4.0, 5.0, 6.0⟩).approx ⟨4.0, 10.0, 18.0⟩ 0.0)

  -- 2. Mat3 arithmetic & operators; A has columns (1,2,3), (4,5,6), (7,8,10)
  let A : Mat3 := ⟨⟨1.0, 2.0, 3.0⟩, ⟨4.0, 5.0, 6.0⟩, ⟨7.0, 8.0, 10.0⟩⟩
  let B : Mat3 := ⟨⟨2.0, 0.0, 1.0⟩, ⟨1.0, 3.0, 2.0⟩, ⟨0.0, 1.0, 4.0⟩⟩

  assertClose "mat det" A.det (-3.0) 0.0
  assertClose "mat trace" A.trace 16.0 0.0
  assertTrue "mat add"
    ((A + B).approx
      ⟨⟨3.0, 2.0, 4.0⟩, ⟨5.0, 8.0, 8.0⟩, ⟨7.0, 9.0, 14.0⟩⟩ 0.0)
  assertTrue "mat sub"
    ((A - B).approx
      ⟨⟨-1.0, 2.0, 2.0⟩, ⟨3.0, 2.0, 4.0⟩, ⟨7.0, 7.0, 6.0⟩⟩ 0.0)
  assertTrue "mat neg"
    ((-A).approx
      ⟨⟨-1.0, -2.0, -3.0⟩, ⟨-4.0, -5.0, -6.0⟩, ⟨-7.0, -8.0, -10.0⟩⟩ 0.0)
  assertTrue "mat smul"
    ((2.0 • A).approx
      ⟨⟨2.0, 4.0, 6.0⟩, ⟨8.0, 10.0, 12.0⟩, ⟨14.0, 16.0, 20.0⟩⟩ 0.0)
  assertTrue "mat div"
    ((A / 2.0).approx (0.5 • A) 0.0)
  assertTrue "mat transpose ᵀ"
    (Aᵀ.approx
      (Mat3.ofRows (⟨1.0, 2.0, 3.0⟩ : Vec3) ⟨4.0, 5.0, 6.0⟩ ⟨7.0, 8.0, 10.0⟩) 0.0)
  assertTrue "mat transpose involution"
    (Aᵀᵀ.approx A 0.0)
  assertTrue "mat mulVec M * v"
    ((A * (⟨1.0, 1.0, 1.0⟩ : Vec3)).approx ⟨12.0, 15.0, 19.0⟩ 0.0)
  assertTrue "vec mulMat v * M"
    (((⟨1.0, 1.0, 1.0⟩ : Vec3) * A).approx ⟨6.0, 15.0, 25.0⟩ 0.0)
  assertTrue "vec mulMat consistency with transpose"
    (((⟨1.0, 1.0, 1.0⟩ : Vec3) * A).approx (Aᵀ * (⟨1.0, 1.0, 1.0⟩ : Vec3)) 0.0)
  assertTrue "mat mul M * N"
    ((A * B).approx
      ⟨⟨9.0, 12.0, 16.0⟩, ⟨27.0, 33.0, 41.0⟩, ⟨32.0, 37.0, 46.0⟩⟩ 0.0)
  assertTrue "mat mul self"
    ((A * A).approx
      ⟨⟨30.0, 36.0, 45.0⟩, ⟨66.0, 81.0, 102.0⟩, ⟨109.0, 134.0, 169.0⟩⟩ 0.0)
  assertTrue "mat id"
    ((Mat3.id * A).approx A 0.0)
  assertClose "mat Frobenius dot ⬝ₘ" (A ⬝ₘ B) 84.0 0.0
  assertClose "mat Frobenius dot self is normSq" (A ⬝ₘ A) ‖A‖² 0.0
  assertClose "mat normSq ‖M‖²" ‖A‖² 304.0 0.0
  assertClose "mat norm ‖M‖" ‖A‖ (Float.sqrt 304.0) 0.0
  assertTrue "mat hadamard ⊙"
    ((A ⊙ B).approx
      ⟨⟨2.0, 0.0, 3.0⟩, ⟨4.0, 15.0, 12.0⟩, ⟨0.0, 8.0, 40.0⟩⟩ 0.0)
  assertTrue "mat abs |M|"
    (|(⟨⟨-1.0, 2.0, -3.0⟩, ⟨4.0, -5.0, 6.0⟩, ⟨-7.0, 8.0, -10.0⟩⟩ : Mat3)|.approx A 0.0)

  let invA := A⁻¹
  assertTrue "mat inv A * A⁻¹ = I"
    ((A * invA).approx Mat3.id 1e-14)
  assertTrue "mat inv A⁻¹ * A = I"
    ((invA * A).approx Mat3.id 1e-14)
  assertTrue "mat inv transpose parses M⁻¹ᵀ"
    (A⁻¹ᵀ.approx invAᵀ 0.0)

  assertTrue "symm toMat3"
    ((SymmMat3.toMat3 workedExample).approx
      ⟨⟨2.0, 1.0, 0.0⟩, ⟨1.0, 2.0, 1.0⟩, ⟨0.0, 1.0, 2.0⟩⟩ 0.0)

  -- 3. Powers: x ^ⁿ n, M ^ⁿ n, and M ^ n
  assertClose "float abs |x|" |-3.5| 3.5 0.0
  assertClose "float powNat 0" (3.0 ^ⁿ 0) 1.0 0.0
  assertClose "float powNat 1" (3.0 ^ⁿ 1) 3.0 0.0
  assertClose "float powNat 2" (3.0 ^ⁿ 2) 9.0 0.0
  assertClose "float powNat 3" (3.0 ^ⁿ 3) 27.0 0.0
  assertClose "float powNat 8 (small boundary)" (2.0 ^ⁿ 8) 256.0 0.0
  assertClose "float powNat 9 (large branch)" (2.0 ^ⁿ 9) 512.0 0.0
  assertClose "float powNat 10" (2.0 ^ⁿ 10) 1024.0 0.0
  assertClose "float powNat negative odd" ((-2.0) ^ⁿ 3) (-8.0) 0.0
  assertClose "float powNat negative even" ((-2.0) ^ⁿ 4) 16.0 0.0

  -- Mat3 powers via PowNat (^ⁿ) — exponentiation by squaring
  assertTrue "mat powNat 0"
    ((A ^ⁿ 0).approx Mat3.id 0.0)
  assertTrue "mat powNat 1"
    ((A ^ⁿ 1).approx A 0.0)
  assertTrue "mat powNat 2"
    ((A ^ⁿ 2).approx (A * A) 0.0)
  assertTrue "mat powNat 3"
    ((A ^ⁿ 3).approx (A * A * A) 0.0)
  assertTrue "mat powNat 4 (even squaring)"
    ((A ^ⁿ 4).approx ((A * A) * (A * A)) 0.0)
  assertTrue "mat powNat 5 (odd squaring)"
    ((A ^ⁿ 5).approx (A * (A ^ⁿ 4)) 0.0)
  assertTrue "mat powNat 8"
    ((A ^ⁿ 8).approx ((A ^ⁿ 4) * (A ^ⁿ 4)) 0.0)

  -- Mat3 powers via standard Pow (^)
  assertTrue "mat pow 0"
    ((A ^ 0).approx Mat3.id 0.0)
  assertTrue "mat pow 1"
    ((A ^ 1).approx A 0.0)
  assertTrue "mat pow 2"
    ((A ^ 2).approx (A * A) 0.0)
  assertTrue "mat pow 3"
    ((A ^ 3).approx (A * A * A) 0.0)
  assertTrue "mat pow 4"
    ((A ^ 4).approx ((A * A) * (A * A)) 0.0)

  -- 4. Precedence & identity rules
  assertClose "powNat mul precedence"
    (2.0 * 3.0 ^ⁿ 2)
    18.0
    0.0
  assertClose "powNat add precedence"
    (1.0 + 3.0 ^ⁿ 2)
    10.0
    0.0
  assertTrue "powNat smul precedence"
    ((2.0 • A ^ⁿ 2).approx (2.0 • (A * A)) 0.0)
  assertClose "powNat frobenius dot precedence"
    (A ^ⁿ 2 ⬝ₘ B)
    ((A * A) ⬝ₘ B)
    0.0
  assertTrue "powNat transpose rule"
    ((A ^ⁿ 2)ᵀ.approx (Aᵀ ^ⁿ 2) 0.0)
  assertClose "triple product precedence"
    (u ⬝ᵥ v ⨯₃ w)
    (u ⬝ᵥ (v ⨯₃ w))
    0.0
  assertTrue "smul cross precedence"
    ((s • u ⨯₃ v).approx (s • (u ⨯₃ v)) 0.0)
  assertClose "dot add precedence"
    (u ⬝ᵥ u + v ⬝ᵥ v)
    ((u ⬝ᵥ u) + (v ⬝ᵥ v))
    0.0
  assertTrue "transpose product identity"
    (((A * B)ᵀ).approx (Bᵀ * Aᵀ) 0.0)
  assertClose "matVec dot"
    ((A * u) ⬝ᵥ v)
    (u ⬝ᵥ (Aᵀ * v))
    0.0

  -- 5. Worked example: eigenvalues 2−√2, 2, 2+√2
  let d := eigendecomp workedExample
  assertClose "worked λ₀" d.eigvals.l₀ 0.5857864376269049 (certTol workedExample)
  assertClose "worked λ₁" d.eigvals.l₁ 2.0 (certTol workedExample)
  assertClose "worked λ₂" d.eigvals.l₂ 3.4142135623730951 (certTol workedExample)
  assertTrue "worked ordered" d.eigvals.isOrdered

  -- 6. Diagonal: eigenvalues (1, 2, 3) with axis eigenvectors (up to sign)
  let dd := eigendecomp diagonalCase
  assertClose "diag λ₀" dd.eigvals.l₀ 1.0 (certTol diagonalCase)
  assertClose "diag λ₁" dd.eigvals.l₁ 2.0 (certTol diagonalCase)
  assertClose "diag λ₂" dd.eigvals.l₂ 3.0 (certTol diagonalCase)
  assertClose "diag axis 0" dd.eigvecs.c₀.x.abs 1.0 (certTol diagonalCase)
  assertClose "diag axis 1" dd.eigvecs.c₁.y.abs 1.0 (certTol diagonalCase)
  assertClose "diag axis 2" dd.eigvecs.c₂.z.abs 1.0 (certTol diagonalCase)

  -- 7. Zero matrix: exact fast path, identity basis
  let dz := eigendecomp zeroMatrix
  assertClose "zero λ₀" dz.eigvals.l₀ 0.0 0.0
  assertClose "zero λ₁" dz.eigvals.l₁ 0.0 0.0
  assertClose "zero λ₂" dz.eigvals.l₂ 0.0 0.0
  assertTrue "zero basis" (dz.eigvecs.approx Mat3.id 0.0)

  -- 8. Scaled identity: exact triple eigenvalue 3.7
  let ds := eigendecomp scaledIdentity
  assertClose "scaled id λ₀" ds.eigvals.l₀ 3.7 0.0
  assertClose "scaled id λ₁" ds.eigvals.l₁ 3.7 0.0
  assertClose "scaled id λ₂" ds.eigvals.l₂ 3.7 0.0
  assertTrue "scaled id ordered" ds.eigvals.isOrdered

  -- 9. Notation end-to-end: change of basis and back is the identity on g
  -- Roundtrip budget: two mat-vecs (≈6ε) + Q's non-orthogonality (≈16ε),
  -- scaled by ‖g‖∞ — rounded to 32ε·maxAbs(g).
  let g : Vec3 := ⟨1.0, 2.0, -1.0⟩
  assertTrue "basis roundtrip"
    ((d.eigvecs * (d.eigvecsᵀ * g)).approx g (32.0 * Float.eps * g.maxAbs))

  -- 10. Coercion: Eigval3 elaborates as Vec3 where a Vec3 is expected
  let wEig : Vec3 := d.eigvals
  assertClose "coe l₀" wEig.x d.eigvals.l₀ 0.0
  assertClose "coe l₁" wEig.y d.eigvals.l₁ 0.0
  assertClose "coe l₂" wEig.z d.eigvals.l₂ 0.0

  -- 11. Coercion: SymmMat3 elaborates as Mat3 where a Mat3 is expected
  let wMat : Mat3 := workedExample
  assertClose "coe M₀₀" wMat.c₀.x workedExample.a₀₀ 0.0
  assertClose "coe M₁₀" wMat.c₀.y workedExample.a₀₁ 0.0
  assertClose "coe M₂₀" wMat.c₀.z workedExample.a₀₂ 0.0
  assertClose "coe M₀₁" wMat.c₁.x workedExample.a₀₁ 0.0
  assertClose "coe M₁₁" wMat.c₁.y workedExample.a₁₁ 0.0
  assertClose "coe M₂₁" wMat.c₁.z workedExample.a₁₂ 0.0
  assertClose "coe M₀₂" wMat.c₂.x workedExample.a₀₂ 0.0
  assertClose "coe M₁₂" wMat.c₂.y workedExample.a₁₂ 0.0
  assertClose "coe M₂₂" wMat.c₂.z workedExample.a₂₂ 0.0

end Eig3x3.Tests
