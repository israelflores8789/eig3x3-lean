/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3
import all Eig3x3.Basic
import Tests.Util

/-!
# Tests.Properties — invariant battery at N deterministic draws

Reference-free invariants of a correct decomposition, checked on thousands
of random symmetric matrices from a deterministic PRNG (splitmix64 — pure
UInt64 arithmetic, bit-identical on every platform since no libm is
involved), including:
* eigenvalue ordering contract,
* three certificates under the shared gates,
* trace and determinant identities,
* right-handedness,
* and bit-exact power-of-two scale invariance.
-/

namespace Eig3x3.Tests

open scoped Eig3x3

/-- splitmix64: a tiny deterministic PRNG. -/
structure Prng where
  s : UInt64

def Prng.next (p : Prng) : UInt64 × Prng :=
  let s := p.s + 0x9E3779B97F4A7C15
  let z := (s ^^^ (s >>> (30 : UInt64))) * 0xBF58476D1CE4E5B9
  let z := (z ^^^ (z >>> (27 : UInt64))) * 0x94D049BB133111EB
  (z ^^^ (z >>> (31 : UInt64)), ⟨s⟩)

/-- Uniform in [0, 1): top 53 bits over 2⁵³ (ofNat is exact below 2⁵³, and
    division by a power of two is exact). -/
def Prng.float01 (p : Prng) : Float × Prng :=
  let (z, p') := p.next
  (Float.ofNat (z >>> (11 : UInt64)).toNat / 9007199254740992.0, p')

/-- Uniform in [−1, 1). -/
def Prng.float11 (p : Prng) : Float × Prng :=
  let (u, p') := p.float01
  (2.0 * u - 1.0, p')

/-- A random symmetric matrix with entries in [−1, 1). -/
def Prng.symmMat3 (p : Prng) : SymmMat3 × Prng :=
  let (x₀, p) := p.float11
  let (x₁, p) := p.float11
  let (x₂, p) := p.float11
  let (x₃, p) := p.float11
  let (x₄, p) := p.float11
  let (x₅, p) := p.float11
  (⟨x₀, x₁, x₂, x₃, x₄, x₅⟩, p)

public def runProperties : IO Unit := do
  let n := 5000
  let mut prng : Prng := ⟨0x9E3779B97F4A7C15⟩
  let mut failures : Array String := #[]
  for _ in [0:n] do
    let (A, p') := prng.symmMat3
    prng := p'
    let d := eigendecomp A
    let e := d.eigvals
    let c := certify A d
    let tol := certTol A

    -- The ordering contract (enforced by the sort in `eigvals`).
    if !e.isOrdered then
      failures := failures.push s!"ordering: ({e.l₀}, {e.l₁}, {e.l₂})"

    -- Certificates under the shared gates.
    if c.maxResidual > tol then
      failures := failures.push s!"residual {c.maxResidual} > {tol}"
    if c.orthogonality > orthoTol then
      failures := failures.push s!"orthogonality {c.orthogonality} > {orthoTol}"
    if c.reconstruction > tol then
      failures := failures.push s!"reconstruction {c.reconstruction} > {tol}"

    -- Trace invariant (the sum adds ~2ε of its own rounding, inside the gate).
    let tr := (A.a₀₀ + A.a₁₁) + A.a₂₂
    let se := (e.l₀ + e.l₁) + e.l₂
    if (se - tr).abs > tol then
      failures := failures.push s!"trace: |{se} − {tr}| > {tol}"

    -- Determinant invariant; the scale is maxAbs³ (budget 64 covers det
    -- evaluation ≈5ε, the product ≈3ε, eigenvalue error ≈36ε at that scale).
    let mA := A.maxAbsEntry
    let prod := (e.l₀ * e.l₁) * e.l₂
    if (prod - A.toMat3.det).abs > 64.0 * Float.eps * mA ^ⁿ 3 then
      failures := failures.push s!"det: |{prod} − {A.toMat3.det}|"

    -- Right-handedness.
    if (d.eigvecs.det - 1.0).abs > orthoTol then
      failures := failures.push s!"handedness: det Q = {d.eigvecs.det}"

    -- Bit-exact power-of-two scale invariance: the preconditioner makes
    -- scaling by 2^k significand-transparent, so this is `==`, not a gate.
    for s in ([0.5, 2.0] : List Float) do
      let eB := (eigendecomp (A.scale s)).eigvals
      if !(eB.l₀ == s * e.l₀ && eB.l₁ == s * e.l₁ && eB.l₂ == s * e.l₂) then
        failures := failures.push s!"scale-invariance ×{s}"

  if failures.size > 0 then
    IO.eprintln s!"FAIL properties: {failures.size} failure(s) over {n} draws"
    for f in failures.take 5 do
      IO.eprintln s!"  {f}"
    throw (IO.userError "properties suite failed")
  else
    IO.println s!"ok   properties: {n} draws, 0 failures"

end Eig3x3.Tests
