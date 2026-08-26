module

import Eig3x3
import Tests.Util

/-!
# Tests.Certificates — runtime evidence assertions

Every decomposition in the zoo must come with residual, orthonormality, and
reconstruction certificates at machine precision, a right-handed Q
(`det Q = 1`), ordered eigenvalues, and `QᵀQ = I`. Observed certificate
values are ≈ 1e-16 · ‖A‖; the 1e-12 gates are deliberately loose so the
suite is robust, while remaining two orders of magnitude below anything a
regression could sneak past.
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
    assertNear s!"{name} residual" c.maxResidual 1e-12
    assertNear s!"{name} orthogonality" c.orthogonality 1e-12
    assertNear s!"{name} reconstruction" c.reconstruction 1e-12
    assertClose s!"{name} right-handed" d.eigvecs.det 1.0 1e-12
    assertTrue s!"{name} ordered" d.eigvals.isOrdered
    assertTrue s!"{name} QᵀQ = I" ((d.eigvecsᵀ * d.eigvecs).approx Mat3.id 1e-12)

end Eig3x3.Tests
