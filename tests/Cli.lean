/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
import Eig3x3
import JsonMiniReader

/-!
# Cli.Main — `eig3x3_cli`: a JSON batch interface to the library

Contract (consumed by `parity/lean_cli.py`; one process per batch, never
per matrix):

  ```text
  in : {"matrices": [[a00, a11, a22, a01, a02, a12], ...]}
  out: {"results": [{"matrix": [a00, ..., a12],
                     "eigvals": [l1, l2, l3],
                     "eigvecs": [c1x, c1y, c1z, c2x, ..., c3z],
                     "certificates": {"maxResidual": r,
                                      "orthogonality": o,
                                      "reconstruction": rc}},
                    ...]}
  ```

The six input entries are `Eig3x3.SymmMat3` field order. `eigvecs` is
column-major: column `i` is the unit eigenvector for `l_i`.
-/

namespace Eig3x3.Cli

open JsonMiniReader

/-! ## Exact decimal printing -/

/-- Drop trailing zeros. Apply to fraction digits only, never to an integer
    part. -/
def dropTrailingZeros (s : String) : String :=
  String.ofList (s.toList.reverse.dropWhile (· == '0')).reverse

/-- Exact decimal expansion of a finite float. Every double is a dyadic
    rational `sig · 2^j`; we print it exactly with big-integer arithmetic,
    so the string round-trips to exactly the same bits. Non-finite values
    print as `null` — JSON has no NaN/Infinity, and the harness treats
    `null` as the failure it is. -/
def exactDecimal (x : Float) : String :=
  if x == 0.0 then "0.0"
  else if !(x.abs ≤ 1.7976931348623157e308) then "null"  -- NaN or ±∞
  else
    let (m, e) := x.frExp          -- x = m · 2^e, |m| ∈ [0.5, 1)
    let neg := m < 0.0
    -- |m| has at most 53 fraction bits, so |m| · 2^53 is an exact integer.
    let sig := (m.abs * 9007199254740992.0).toUInt64.toNat
    let j := e - 53                -- |x| = sig · 2^j
    let body :=
      if j ≥ 0 then
        toString (sig * 2 ^ j.toNat) ++ ".0"
      else
        -- sig / 2^k = sig · 5^k / 10^k, with k = 53 - e
        let k := (53 - e).toNat
        let num := toString (sig * 5 ^ k)
        let len := num.length
        if len > k then
          let frac := dropTrailingZeros (num.drop (len - k)).toString
          (num.take (len - k)).toString ++ "." ++ (if frac.isEmpty then "0" else frac)
        else
          "0." ++ String.ofList (List.replicate (k - len) '0')
                ++ dropTrailingZeros num
    (if neg then "-" else "") ++ body

/-! ## Output -/

def parseInput (s : String) : Option (List SymmMat3) := do
  let cs ← expectChar '{' s.toList
  let cs ← expectString "matrices" cs
  let cs ← expectChar ':' cs
  let cs ← expectChar '[' cs
  let (ms, cs) ← parseSymmMat3Rows [] cs
  let _ ← expectChar '}' cs
  pure ms

/-- One result object, exactly matching the `lean_cli.py` schema. -/
def resultJson (A : SymmMat3) : String :=
  let d := eigendecomp A
  let c := certify A d
  let e := d.eigvals
  let q := d.eigvecs
  let af := [A.a00, A.a11, A.a22, A.a01, A.a02, A.a12].map exactDecimal
  let evs := [e.l₁, e.l₂, e.l₃].map exactDecimal
  let qf := [q.c₁.x, q.c₁.y, q.c₁.z,
             q.c₂.x, q.c₂.y, q.c₂.z,
             q.c₃.x, q.c₃.y, q.c₃.z].map exactDecimal
  "{\"matrix\":[" ++ String.intercalate "," af
    ++ "],\"eigvals\":[" ++ String.intercalate "," evs
    ++ "],\"eigvecs\":[" ++ String.intercalate "," qf
    ++ "],\"certificates\":{\"maxResidual\":" ++ exactDecimal c.maxResidual
    ++ ",\"orthogonality\":" ++ exactDecimal c.orthogonality
    ++ ",\"reconstruction\":" ++ exactDecimal c.reconstruction ++ "}}"

def run : IO UInt32 := do
  let stdin ← IO.getStdin
  let input ← stdin.readToEnd
  match parseInput input with
  | none =>
    IO.eprintln "eig3x3_cli: expected {\"matrices\": [[a00,a11,a22,a01,a02,a12], ...]} on stdin"
    return 1
  | some matrices =>
    IO.print "{\"results\":["
    let mut first := true
    for A in matrices do
      if first then first := false else IO.print ","
      IO.print (resultJson A)
    IO.println "]}"
    return 0

end Eig3x3.Cli

def main : IO UInt32 := Eig3x3.Cli.run
