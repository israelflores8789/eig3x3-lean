/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3
import Tests.Util
public import Eig3x3.Basic

/-!
# Tests.Golden — high-precision golden vectors

References computed at 50-digit precision with mpmath (`parity/golden.py`)
and embedded here as 17-significant-digit literals, which round-trip to the
correct float64. This module is the Lean suite's high-precision anchor: it
validates `eigendecomp` against essentially-exact references with no Python
at test time. Gates are the shared `certTol` standard (64ε · max |entry|).

Two notes on the literals:
* simple decimals (`3.7`, `1.0e-8`) parse exactly, but long mantissas and
  extreme exponents rely on the toolchain's float parsing, which may differ
  from an ideal rounding by ~1 ulp — four orders of magnitude below the gates.
* And `golden.py`'s curated set **must** stay in sync with this file;
  regenerate with `just golden`, **never** edit literals by hand.
-/

namespace Eig3x3.Tests

/-- (name, matrix, correctly-rounded eigenvalues `l₁ ≤ l₂ ≤ l₃`). -/
public def goldenCases : List (String × SymmMat3 × Eigval3) :=
  [("worked", ⟨2.0, 2.0, 2.0, 1.0, 0.0, 1.0⟩,
    ⟨0.58578643762690497, 2.0, 3.4142135623730949⟩),
   ("diagonal", ⟨1.0, 2.0, 3.0, 0.0, 0.0, 0.0⟩,
    ⟨1.0, 2.0, 3.0⟩),
   ("scaled-id", ⟨3.7, 3.7, 3.7, 0.0, 0.0, 0.0⟩,
    ⟨3.7000000000000002, 3.7000000000000002, 3.7000000000000002⟩),
   ("zero", ⟨0.0, 0.0, 0.0, 0.0, 0.0, 0.0⟩,
    ⟨0.0, 0.0, 0.0⟩),
   ("near-double", ⟨1.0, 1.0, -1.0, 1.0e-8, 0.0, 0.0⟩,
    ⟨-1.0, 0.99999998999999995, 1.0000000099999999⟩),
   ("regression", ⟨-2.0, 2.0, -9.0, 9.0, 8.0, 6.0⟩,
    ⟨-14.233076938854385, -8.2169045134999497, 13.449981452354335⟩),
   ("path-near-double", ⟨-1.0, 1.0, 1.0000000099999999, 0.0, 0.0, 0.0⟩,
    ⟨-1.0, 1.0, 1.0000000099999999⟩),
   ("path-near-triple", ⟨1.0, 1.0, 1.0000000099999999, 0.0, 0.0, 0.0⟩,
    ⟨1.0, 1.0, 1.0000000099999999⟩),
   ("frontier-huge-id", ⟨1.0e300, 1.0e300, 1.0e300, 0.0, 0.0, 0.0⟩,
    ⟨1.0000000000000001e300, 1.0000000000000001e300, 1.0000000000000001e300⟩),
   ("frontier-tiny-id", ⟨1.0e-300, 1.0e-300, 1.0e-300, 0.0, 0.0, 0.0⟩,
    ⟨1.0e-300, 1.0e-300, 1.0e-300⟩)]

public def runGolden : IO Unit := do
  for (name, A, ref) in goldenCases do
    let d := eigendecomp A
    let tol := certTol A
    assertClose s!"{name} golden l₁" d.eigvals.l₁ ref.l₁ tol
    assertClose s!"{name} golden l₂" d.eigvals.l₂ ref.l₂ tol
    assertClose s!"{name} golden l₃" d.eigvals.l₃ ref.l₃ tol

end Eig3x3.Tests
