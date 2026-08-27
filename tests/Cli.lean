/-
Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Israel Flores-Arbolay
-/
module

import Eig3x3

/-!
# Cli.Main — `eig3x3-cli`: a JSON batch interface to the library

Contract (consumed by `parity/lean_cli.py`; one process per batch, never
per matrix):

    in : {"matrices": [[a00, a11, a22, a01, a02, a12], ...]}
    out: {"results": [{"eigvals": [l1, l2, l3],
                       "eigvecs": [c1x, c1y, c1z, c2x, ..., c3z],
                       "certificates": {"maxResidual": r,
                                        "orthogonality": o,
                                        "reconstruction": rc}}, ...]}

The six input entries are `Eig3x3.SymmMat3` field order. `eigvecs` is
column-major: column `i` is the unit eigenvector for `l_i`.

Two pieces of infrastructure live here because Lean core lacks them: a
minimal JSON reader for the fixed input schema (core has no string→Float
parser yet — leanprover/lean4#14659), and an exact decimal printer for
`Float` (the default display is fixed-precision and lossy; every double has
a finite exact decimal expansion, computed here with big-integer
arithmetic).

This executable is also the library's first true downstream consumer: it
imports only the public facade.
-/

namespace Eig3x3.Cli

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

/-! ## Minimal JSON reader (fixed schema) -/

def isWs (c : Char) : Bool :=
  c == ' ' || c == '\n' || c == '\t' || c == '\r'

def skipWs : List Char → List Char
  | c :: cs => if isWs c then skipWs cs else c :: cs
  | [] => []

def expectChar (c : Char) (cs : List Char) : Option (List Char) :=
  match skipWs cs with
  | c' :: cs' => if c == c' then some cs' else none
  | [] => none

def expectString (s : String) (cs : List Char) : Option (List Char) :=
  match skipWs cs with
  | '"' :: cs => go s.toList cs
  | _ => none
where
  go : List Char → List Char → Option (List Char)
    | [], '"' :: cs => some cs
    | c :: pat, c' :: cs => if c == c' then go pat cs else none
    | _, _ => none

def isNumChar (c : Char) : Bool :=
  c.isDigit || c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E'

def takeNumTok : List Char → List Char → List Char × List Char
  | acc, c :: cs =>
    if isNumChar c then takeNumTok (c :: acc) cs else (acc.reverse, c :: cs)
  | acc, [] => (acc.reverse, [])

/-- Parse a decimal/scientific literal via `OfScientific` (mantissa and
    decimal exponent as exact `Nat`s, then one Float scaling). May differ
    from ideal rounding by ~1 ulp — four orders of magnitude below the
    parity gates, and harmless because symmetric eigenvalues are perfectly
    conditioned. -/
def parseDecimal (s : String) : Option Float := do
  let (mant, expStr) ←
    match (s.replace "E" "e").splitOn "e" with
    | [m] => some (m, "0")
    | [m, e] => some (m, e)
    | _ => none
  let expStr := if expStr.startsWith "+" then (expStr.drop 1).toString else expStr
  let neg := mant.startsWith "-"
  let mantAbs := if neg then (mant.drop 1).toString else mant
  let (intPart, fracPart) ←
    match mantAbs.splitOn "." with
    | [i] => some (i, "")
    | [i, f] => some (i, f)
    | _ => none
  let mantissa ← (intPart ++ fracPart).toNat?
  let e10 ← expStr.toInt?
  let eFinal := e10 - Int.ofNat fracPart.length
  let f : Float := OfScientific.ofScientific mantissa (eFinal < 0) eFinal.natAbs
  pure (if neg then -f else f)

def parseNumber (cs : List Char) : Option (Float × List Char) := do
  let (tok, rest) := takeNumTok [] (skipWs cs)
  if tok.isEmpty then none else pure (← parseDecimal (String.ofList tok), rest)

/-- One matrix row: `[a00, a11, a22, a01, a02, a12]`. -/
def parseRow (cs : List Char) : Option (SymmMat3 × List Char) := do
  let cs ← expectChar '[' cs
  let (x0, cs) ← parseNumber cs
  let cs ← expectChar ',' cs
  let (x1, cs) ← parseNumber cs
  let cs ← expectChar ',' cs
  let (x2, cs) ← parseNumber cs
  let cs ← expectChar ',' cs
  let (x3, cs) ← parseNumber cs
  let cs ← expectChar ',' cs
  let (x4, cs) ← parseNumber cs
  let cs ← expectChar ',' cs
  let (x5, cs) ← parseNumber cs
  let cs ← expectChar ']' cs
  pure (⟨x0, x1, x2, x3, x4, x5⟩, cs)

/-- The `"matrices"` array. `partial`: a hand-rolled recursive-descent loop
    whose termination the compiler cannot see; the input is
    machine-generated, so failure returns `none`, never a crash. -/
partial def parseRows (acc : List SymmMat3) (cs : List Char)
    : Option (List SymmMat3 × List Char) :=
  match skipWs cs with
  | ']' :: cs => some (acc.reverse, cs)
  | _ =>
    match parseRow cs with
    | none => none
    | some (m, cs) =>
      match skipWs cs with
      | ',' :: cs => parseRows (m :: acc) cs
      | ']' :: cs => some ((m :: acc).reverse, cs)
      | _ => none

def parseInput (s : String) : Option (List SymmMat3) := do
  let cs ← expectChar '{' s.toList
  let cs ← expectString "matrices" cs
  let cs ← expectChar ':' cs
  let cs ← expectChar '[' cs
  let (ms, cs) ← parseRows [] cs
  let _ ← expectChar '}' cs
  pure ms

/-! ## Output -/

/-- One result object, exactly matching the `lean_cli.py` schema. -/
def resultJson (A : SymmMat3) : String :=
  let d := eigendecomp A
  let c := certify A d
  let e := d.eigvals
  let q := d.eigvecs
  let evs := [e.l₁, e.l₂, e.l₃].map exactDecimal
  let qf := [q.c₁.x, q.c₁.y, q.c₁.z,
             q.c₂.x, q.c₂.y, q.c₂.z,
             q.c₃.x, q.c₃.y, q.c₃.z].map exactDecimal
  "{\"eigvals\":[" ++ String.intercalate "," evs
    ++ "],\"eigvecs\":[" ++ String.intercalate "," qf
    ++ "],\"certificates\":{\"maxResidual\":" ++ exactDecimal c.maxResidual
    ++ ",\"orthogonality\":" ++ exactDecimal c.orthogonality
    ++ ",\"reconstruction\":" ++ exactDecimal c.reconstruction ++ "}}"

def run : IO UInt32 := do
  let stdin ← IO.getStdin
  let input ← stdin.readToEnd
  match parseInput input with
  | none =>
    IO.eprintln "eig3x3-cli: expected {\"matrices\": [[a00,a11,a22,a01,a02,a12], ...]} on stdin"
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
