/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3
import all Eig3x3.Eigenvalues
import Tests.Util

/-!
# Tests.Regression — pinned historical failures

The r₁₀ transcription regression (exact discriminant value) and the
near-double path that distinguishes Habera-Zilian's Algorithm 8
discriminant from the naive representation.
-/

namespace Eig3x3.Tests

/-- Naive discriminant Δ = 4J₂³ − 27J₃², clamped to [0, ∞). Kept for
    benchmarking only (reproduces the paper's naive-vs-present comparisons);
    suffers catastrophic cancellation near double eigenvalues with finite J₂
    (observed ≈5e-9 absolute eigenvalue error on the D2 path). -/
def deltaNaive (J2 J3 : Float) : Float :=
  let d := 4.0 * J2 * J2 * J2 - 27.0 * J3 * J3
  if d < 0.0 then 0.0 else d

public def runRegression : IO Unit := do
  -- Exact discriminant on the matrix that caught the r₁₀ sign error
  -- in Algorithm 8 (Habera-Zilian 2025).
  -- 13,021,520 < 2^53, so it is exactly representable in float64.
  assertClose "regression Δ" (delta regressionMatrix) 13021520.0 0.0

  -- Near-double eigenvalues at machine precision.
  -- (the naive Δ gave ≈2.5e-9 absolute eigenvalue error here)
  let d := eigendecomp nearDouble
  assertClose "nearDouble λ₁" d.eigvals.l₁ (-1.0) 1e-12
  assertClose "nearDouble λ₂" d.eigvals.l₂ 0.99999999 1e-12
  assertClose "nearDouble λ₃" d.eigvals.l₃ 1.00000001 1e-12
  assertTrue "nearDouble ordered" d.eigvals.isOrdered

  -- Informational: present vs naive discriminant on the near-double path.
  IO.println s!"info nearDouble Δ (Alg. 8): {delta nearDouble}"
  IO.println s!"info nearDouble Δ (naive):  {deltaNaive (j2 nearDouble) (j3 nearDouble)}"

end Eig3x3.Tests
