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
scaled identity), and end-to-end exercise of the notation (`ᵀ`, `⬝`, `*`,
`•`) and the `Eigval3 → Vec3` coercion as a downstream consumer would use
them.
-/

namespace Eig3x3.Tests

open scoped Eig3x3

public def runKnownAnswer : IO Unit := do
  -- Vec3 arithmetic (all exact)
  assertTrue "vec add"
    (((⟨1.0, 2.0, 3.0⟩ : Vec3) + ⟨4.0, 5.0, 6.0⟩).approx ⟨5.0, 7.0, 9.0⟩ 0.0)
  assertTrue "vec sub"
    (((⟨4.0, 5.0, 6.0⟩ : Vec3) - ⟨1.0, 2.0, 3.0⟩).approx ⟨3.0, 3.0, 3.0⟩ 0.0)
  assertClose "vec dot" ((⟨1.0, 2.0, 3.0⟩ : Vec3) ⬝ ⟨4.0, 5.0, 6.0⟩) 32.0 0.0
  assertTrue "vec cross"
    (((⟨1.0, 0.0, 0.0⟩ : Vec3).cross ⟨0.0, 1.0, 0.0⟩).approx ⟨0.0, 0.0, 1.0⟩ 0.0)
  assertClose "vec normSq" (Vec3.normSq ⟨1.0, 2.0, 3.0⟩) 14.0 0.0
  assertTrue "vec smul"
    ((2.0 • (⟨1.0, 2.0, 3.0⟩ : Vec3)).approx ⟨2.0, 4.0, 6.0⟩ 0.0)
  assertTrue "vec map"
    (((⟨-1.0, 2.0, -3.0⟩ : Vec3).map Float.abs).approx ⟨1.0, 2.0, 3.0⟩ 0.0)
  assertTrue "vec map2"
    ((Vec3.map2 (· * ·) ⟨1.0, 2.0, 3.0⟩ ⟨4.0, 5.0, 6.0⟩).approx ⟨4.0, 10.0, 18.0⟩ 0.0)

  -- Mat3 arithmetic (all exact); M has columns (1,2,3), (4,5,6), (7,8,10)
  let M : Mat3 := ⟨⟨1.0, 2.0, 3.0⟩, ⟨4.0, 5.0, 6.0⟩, ⟨7.0, 8.0, 10.0⟩⟩
  assertClose "mat det" M.det (-3.0) 0.0
  assertClose "mat trace" M.trace 16.0 0.0
  assertTrue "mat transpose" (Mᵀ.approx
    (Mat3.ofRows ⟨1.0, 2.0, 3.0⟩ ⟨4.0, 5.0, 6.0⟩ ⟨7.0, 8.0, 10.0⟩) 0.0)
  assertTrue "mat mulVec" ((M * (⟨1.0, 1.0, 1.0⟩ : Vec3)).approx ⟨12.0, 15.0, 19.0⟩ 0.0)
  assertTrue "mat transposeMulVec" ((Mᵀ * (⟨1.0, 1.0, 1.0⟩ : Vec3)).approx ⟨6.0, 15.0, 25.0⟩ 0.0)
  assertTrue "mat mul" ((M * M).approx
    ⟨⟨30.0, 36.0, 45.0⟩, ⟨66.0, 81.0, 102.0⟩, ⟨109.0, 134.0, 169.0⟩⟩ 0.0)
  assertTrue "mat id" ((Mat3.id * M).approx M 0.0)
  assertTrue "symm toMat3" ((SymmMat3.toMat3 workedExample).approx
    ⟨⟨2.0, 1.0, 0.0⟩, ⟨1.0, 2.0, 1.0⟩, ⟨0.0, 1.0, 2.0⟩⟩ 0.0)

  -- Worked example: eigenvalues 2−√2, 2, 2+√2
  let d := eigendecomp workedExample
  assertClose "worked λ₁" d.eigvals.l₁ 0.5857864376269049 (certTol workedExample)
  assertClose "worked λ₂" d.eigvals.l₂ 2.0 (certTol workedExample)
  assertClose "worked λ₃" d.eigvals.l₃ 3.4142135623730951 (certTol workedExample)
  assertTrue "worked ordered" d.eigvals.isOrdered

  -- Diagonal: eigenvalues (1, 2, 3) with axis eigenvectors (up to sign)
  let dd := eigendecomp diagonalCase
  assertClose "diag λ₁" dd.eigvals.l₁ 1.0 (certTol diagonalCase)
  assertClose "diag λ₂" dd.eigvals.l₂ 2.0 (certTol diagonalCase)
  assertClose "diag λ₃" dd.eigvals.l₃ 3.0 (certTol diagonalCase)
  assertClose "diag axis 1" dd.eigvecs.c₁.x.abs 1.0 (certTol diagonalCase)
  assertClose "diag axis 2" dd.eigvecs.c₂.y.abs 1.0 (certTol diagonalCase)
  assertClose "diag axis 3" dd.eigvecs.c₃.z.abs 1.0 (certTol diagonalCase)

  -- Zero matrix: exact fast path, identity basis
  let dz := eigendecomp zeroMatrix
  assertClose "zero λ₁" dz.eigvals.l₁ 0.0 0.0
  assertClose "zero λ₂" dz.eigvals.l₂ 0.0 0.0
  assertClose "zero λ₃" dz.eigvals.l₃ 0.0 0.0
  assertTrue "zero basis" (dz.eigvecs.approx Mat3.id 0.0)

  -- Scaled identity: exact triple eigenvalue 3.7
  let ds := eigendecomp scaledIdentity
  assertClose "scaled id λ₁" ds.eigvals.l₁ 3.7 0.0
  assertClose "scaled id λ₂" ds.eigvals.l₂ 3.7 0.0
  assertClose "scaled id λ₃" ds.eigvals.l₃ 3.7 0.0
  assertTrue "scaled id ordered" ds.eigvals.isOrdered

  -- Notation end-to-end: change of basis and back is the identity on g
  -- Roundtrip budget: two mat-vecs (≈6ε) + Q's non-orthogonality (≈16ε),
  -- scaled by ‖g‖∞ — rounded to 32ε·maxAbs(g).
  let g : Vec3 := ⟨1.0, 2.0, -1.0⟩
  assertTrue "basis roundtrip"
    ((d.eigvecs * (d.eigvecsᵀ * g)).approx g (32.0 * Float.eps * g.maxAbs))

  -- Coercion: Eigval3 elaborates as Vec3 where a Vec3 is expected
  let w : Vec3 := d.eigvals
  assertClose "coe l₁" w.x d.eigvals.l₁ 0.0
  assertClose "coe l₂" w.y d.eigvals.l₂ 0.0
  assertClose "coe l₃" w.z d.eigvals.l₃ 0.0

end Eig3x3.Tests
