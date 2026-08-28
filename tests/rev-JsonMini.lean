module

import Eig3x3

/-!
# JsonMini — a minimal JSON reader for repo infrastructure

Shared by `Cli.lean` (the eig3x3_cli binary) and `Tests/Golden.lean`.
Scope: exactly what our own generated JSON needs — whitespace, string
literals without escapes, decimal/scientific numbers, and integer members
(dyadic pairs). This is not a general JSON library; if a schema ever grows
beyond that, reconsider rather than extend.

Lean core currently has no string→Float parser (leanprover/lean4#14659), so
`parseDecimal` goes through `OfScientific`: mantissa and decimal exponent as
exact `Nat`s, then one Float scaling. That can sit ~1 ulp from ideal
rounding — documented and harmless at our gates — which is why reference
*values* travel as exact dyadic pairs instead (see `Tests/Golden.lean`).
-/

namespace JsonMini

def isWs (c : Char) : Bool :=
  c == ' ' || c == '\n' || c == '\t' || c == '\r'

def skipWs : List Char → List Char
  | c :: cs => if isWs c then skipWs cs else c :: cs
  | [] => []

def expectChar (c : Char) (cs : List Char) : Option (List Char) :=
  match skipWs cs with
  | c' :: cs' => if c == c' then some cs' else none
  | [] => none

/-- Expect the string literal `"s"` exactly. -/
def expectString (s : String) (cs : List Char) : Option (List Char) :=
  match skipWs cs with
  | '"' :: cs => go s.toList cs
  | _ => none
where
  go : List Char → List Char → Option (List Char)
    | [], '"' :: cs => some cs
    | c :: pat, c' :: cs => if c == c' then go pat cs else none
    | _, _ => none

/-- Parse a string literal (no escape sequences — our generators emit plain
    ASCII names). -/
def parseStringLit : List Char → Option (String × List Char)
  | cs => match skipWs cs with
    | '"' :: cs => go [] cs
    | _ => none
where
  go (acc : List Char) : List Char → Option (String × List Char)
    | '"' :: cs => some (String.ofList acc.reverse, cs)
    | c :: cs => go (c :: acc) cs
    | [] => none

def isNumChar (c : Char) : Bool :=
  c.isDigit || c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E'

def takeNumTok : List Char → List Char → List Char × List Char
  | acc, c :: cs =>
    if isNumChar c then takeNumTok (c :: acc) cs else (acc.reverse, c :: cs)
  | acc, [] => (acc.reverse, [])

def takeIntTok : List Char → List Char → List Char × List Char
  | acc, c :: cs =>
    if c.isDigit || c == '-' || c == '+' then takeIntTok (c :: acc) cs
    else (acc.reverse, c :: cs)
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

/-- Parse a Float from the character stream. -/
def parseFloat (cs : List Char) : Option (Float × List Char) := do
  let (tok, rest) := takeNumTok [] (skipWs cs)
  if tok.isEmpty then none else pure (← parseDecimal (String.ofList tok), rest)

/-- Parse an Int from the character stream (dyadic-pair members). -/
def parseInt (cs : List Char) : Option (Int × List Char) := do
  let (tok, rest) := takeIntTok [] (skipWs cs)
  if tok.isEmpty then none
  else
    let s := String.ofList tok
    let s := if s.startsWith "+" then (s.drop 1).toString else s
    pure (← s.toInt?, rest)

/-- One matrix row: `[a00, a11, a22, a01, a02, a12]` — SymmMat3 field
    order. Shared by the CLI's input schema and Golden's cases. -/
def parseSymmMat3 (cs : List Char) : Option (SymmMat3 × List Char) := do
  let cs ← expectChar '[' cs
  let (x0, cs) ← parseFloat cs
  let cs ← expectChar ',' cs
  let (x1, cs) ← parseFloat cs
  let cs ← expectChar ',' cs
  let (x2, cs) ← parseFloat cs
  let cs ← expectChar ',' cs
  let (x3, cs) ← parseFloat cs
  let cs ← expectChar ',' cs
  let (x4, cs) ← parseFloat cs
  let cs ← expectChar ',' cs
  let (x5, cs) ← parseFloat cs
  let cs ← expectChar ']' cs
  pure (⟨x0, x1, x2, x3, x4, x5⟩, cs)

end JsonMini
