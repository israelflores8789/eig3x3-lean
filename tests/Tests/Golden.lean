/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3
import Tests.Util
import JsonMiniReader

/-!
# Tests.Golden — high-precision golden test vectors

Reads `generated/golden.json` (generated eigenvalue "golden vectors" at
50-digit mpmath precision by `parity/golden.py`, git versioned, regenerated
with `just golden`). Each reference eigenvalue travels as an exact dyadic pair
`[sig, exp]` with value `sig · 2^exp`, so the transfer from Python to Lean is
bit-exact without decimal parser rounding.
-/

namespace Eig3x3.Tests

open JsonMiniReader

/-- Reconstruct a float from its exact dyadic pair. `ofInt` is exact below
    2⁵³ and `scaleB` (C `scalbn`) is exact by construction, so this transfer
    is bit-exact and immune to decimal-parser rounding. If a toolchain
    lacks `Float.scaleB`, `Float.ofInt sig * Float.pow 2.0 exp.toFloat` is
    exact for integer powers of two in every real libm (the test gates would
    catch a violation regardless). -/
def ofDyadic (sig : Int) (exp : Int) : Float :=
  Float.scaleB (Float.ofInt sig) exp

/-- One golden case: name, matrix, three reference eigenvalues. -/
structure GoldenCase where
  name : String
  matrix : SymmMat3
  refs : Eigval3

def parseDyadic (cs : List Char) : Option (Float × List Char) := do
  let cs ← expectChar '[' cs
  let (s, cs) ← parseInt cs
  let cs ← expectChar ',' cs
  let (e, cs) ← parseInt cs
  let cs ← expectChar ']' cs
  pure (ofDyadic s e, cs)

def parseRefs (cs : List Char) : Option (Eigval3 × List Char) := do
  let cs ← expectChar '[' cs
  let (a, cs) ← parseDyadic cs
  let cs ← expectChar ',' cs
  let (b, cs) ← parseDyadic cs
  let cs ← expectChar ',' cs
  let (c, cs) ← parseDyadic cs
  let cs ← expectChar ']' cs
  pure (⟨a, b, c⟩, cs)

/-- Skip human-readable display strings. -/
def skipStringArray (cs : List Char) : Option (List Char) := do
  let cs ← expectChar '[' cs
  let (_, cs) ← parseStringLit cs
  let cs ← expectChar ',' cs
  let (_, cs) ← parseStringLit cs
  let cs ← expectChar ',' cs
  let (_, cs) ← parseStringLit cs
  let cs ← expectChar ']' cs
  pure cs

def parseCase (cs : List Char) : Option (GoldenCase × List Char) := do
  let cs ← expectChar '{' cs
  let cs ← expectString "name" cs
  let cs ← expectChar ':' cs
  let (name, cs) ← parseStringLit cs
  let cs ← expectChar ',' cs
  let cs ← expectString "matrix" cs
  let cs ← expectChar ':' cs
  let (m, cs) ← parseSymmMat3 cs
  let cs ← expectChar ',' cs
  let cs ← expectString "eigvals_display" cs
  let cs ← expectChar ':' cs
  let cs ← skipStringArray cs
  let cs ← expectChar ',' cs
  let cs ← expectString "eigvals" cs
  let cs ← expectChar ':' cs
  let (refs, cs) ← parseRefs cs
  let cs ← expectChar '}' cs
  pure ({ name, matrix := m, refs }, cs)

/-- `partial`: a hand-rolled recursive-descent loop whose termination the
    compiler cannot see; the input is machine-generated, so failure returns
    `none`. -/
partial def parseCases (acc : List GoldenCase) (cs : List Char)
    : Option (List GoldenCase) :=
  match skipWs cs with
  | ']' :: _ => some acc.reverse
  | _ =>
    match parseCase cs with
    | none => none
    | some (gc, cs) =>
      match skipWs cs with
      | ',' :: cs => parseCases (gc :: acc) cs
      | ']' :: _ => some (gc :: acc).reverse
      | _ => none

/-- The generator writes "cases" first and "provenance" metadata second;
    we parse the former and ignore the latter. -/
def parseGolden (s : String) : Option (List GoldenCase) := do
  let cs ← expectChar '{' s.toList
  let cs ← expectString "cases" cs
  let cs ← expectChar ':' cs
  let cs ← expectChar '[' cs
  parseCases [] cs

def loadGolden (path : String) : IO (List GoldenCase) := do
  let contents ←
    try
      IO.FS.readFile path
    catch _ =>
      throw (IO.userError
        s!"golden vectors not found at {path}; run `just golden` from the repo root")
  match parseGolden contents with
  | some cs => pure cs
  | none =>
    throw (IO.userError
      s!"could not parse {path}; regenerate with `just golden`")

/-- Each case must meet the shared standard (64ε · max |entry|) against the
    essentially-exact references. -/
public def runGolden : IO Unit := do
  let cases ← loadGolden "generated/golden.json"
  for gc in cases do
    let d := eigendecomp gc.matrix
    let tol := certTol gc.matrix
    assertClose s!"{gc.name} golden l₀" d.eigvals.l₀ gc.refs.l₀ tol
    assertClose s!"{gc.name} golden l₁" d.eigvals.l₁ gc.refs.l₁ tol
    assertClose s!"{gc.name} golden l₂" d.eigvals.l₂ gc.refs.l₂ tol

end Eig3x3.Tests
