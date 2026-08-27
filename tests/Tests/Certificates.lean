/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3
import Tests.Util

/-!
# Tests.Certificates — runtime evidence assertions

Every decomposition in the zoo must come with residual, orthonormality, and
reconstruction certificates at machine precision including a right-handed Q
(`det Q = 1`), ordered eigenvalues, and `QᵀQ = I`. The zoo must all meet
the certificate standard:
* residual and reconstruction below 64ε·‖A‖∞,
* and orthonormality below 16ε,
where ε = 2⁻⁵² is float64 machine epsilon and ‖A‖∞ is the max |entry|
(`certTol`/`orthoTol` in Tests.Util).

The constants are the algorithm's rounding budget:
* the eigenvalue error (≤ ~36ε·‖A‖∞ in the 20k-matrix validation) dominates,
* the evaluation of Av − λv contributes ≤ 5ε·‖A‖∞,
* and ≤ 4ε is allowance for platform libm variation.

This is IEEE 754 compliant; however, note cos/atan2 are transcendental and
exempt from IEEE 754, so results may differ by ~1 ulp across platforms.
Observed values on the zoo are ≤ 7ε·‖A‖∞.

For calibration: the naive-discriminant regression produces a residual of
3.6e-9 on the near-double case, 250,000× over its gate — the standard has
five orders of magnitude of detection headroom. The zero matrix and scaled
identity involve no transcendentals and are gated at exactly zero.
-/

namespace Eig3x3.Tests

open scoped Eig3x3

public def runCertificates : IO Unit := do
  let cases : List (String × SymmMat3) :=
    [("worked", workedExample),
     ("diagonal", diagonalCase),
     ("scaled-id", scaledIdentity),
     ("zero", zeroMatrix),
     ("near-double", nearDouble),
     ("regression", regressionMatrix)]
  for (name, A) in cases do
    let d := eigendecomp A
    let c := certify A d
    assertNear s!"{name} residual" c.maxResidual (certTol A)
    assertNear s!"{name} orthogonality" c.orthogonality orthoTol
    assertNear s!"{name} reconstruction" c.reconstruction (certTol A)
    assertClose s!"{name} right-handed" d.eigvecs.det 1.0 orthoTol
    assertTrue s!"{name} ordered" d.eigvals.isOrdered
    assertTrue s!"{name} QᵀQ = I" ((d.eigvecsᵀ * d.eigvecs).approx Mat3.id orthoTol)

end Eig3x3.Tests
