/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
import Eig3x3

/-!
# Bench — in-process operation timing

Benchmarks pure computation cost of `eigendecomp` + `certify` per matrix, in ns/op.

Methodology: monotonic clock (`IO.monoNanosNow`), a warm-up pass, a fixed mixed
batch cycled round-robin, and an accumulator (`sink`) that depends on every result
so the compiler cannot eliminate the computation as dead code.
-/

namespace Eig3x3.Bench

/-- A fixed mixed batch: worked example, diagonal, near-double, regression,
    scaled identity. -/
def cases : Array SymmMat3 :=
  #[⟨2.0, 2.0, 2.0, 1.0, 0.0, 1.0⟩,
    ⟨1.0, 2.0, 3.0, 0.0, 0.0, 0.0⟩,
    ⟨1.0, 1.0, -1.0, 1.0e-8, 0.0, 0.0⟩,
    ⟨-2.0, 2.0, -9.0, 9.0, 8.0, 6.0⟩,
    ⟨3.7, 3.7, 3.7, 0.0, 0.0, 0.0⟩]

def run (reps : Nat) : IO (Float × Float × Float) := do
  have hcases : 0 < cases.size := by decide
  let mut sink := 0.0
  for i in [0 : reps / 10] do                      -- warm-up
    let d := eigendecomp (cases[i % cases.size]'(Nat.mod_lt i hcases))
    sink := sink + d.eigvals.l₀
  -- Workload 1: eigendecomp alone — the fair LAPACK comparison
  let t0 ← IO.monoNanosNow
  for i in [0:reps] do
    let d := eigendecomp (cases[i % cases.size]'(Nat.mod_lt i hcases))
    sink := sink + d.eigvals.l₀
  let t1 ← IO.monoNanosNow
  -- Workload 2: eigendecomp + certify — the CLI's per-matrix workload.
  -- (end-to-end CLI) − (workload 2) isolates the JSON boundary;
  -- (workload 2) − (workload 1) is what certification costs.
  let t2 ← IO.monoNanosNow
  for i in [0:reps] do
    let A := cases[i % cases.size]'(Nat.mod_lt i hcases)
    let d := eigendecomp A
    let c := certify A d
    sink := sink + d.eigvals.l₀ + c.maxResidual
  let t3 ← IO.monoNanosNow
  let per := Float.ofNat reps
  let decomp := Float.ofNat (t1 - t0) / per
  let full  := Float.ofNat (t3 - t2) / per
  pure (decomp, full, sink)

end Eig3x3.Bench

def main : IO UInt32 := do
  let reps := 200000
  let (decomp, full, sink) ← Eig3x3.Bench.run reps
  IO.println s!"ns_per_op_decomp {decomp}"
  IO.println s!"ns_per_op_full {full}"
  IO.println s!"reps {reps}"
  IO.println s!"sink {sink}"
  return 0
