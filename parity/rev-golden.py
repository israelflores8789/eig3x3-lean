"""golden.py — generate golden.json, the high-precision reference set.

Role in the harness: pure generator. This module contains no assertions.
The Lean suite (tests/Tests/Golden.lean) consumes golden.json directly, and
compare.py may use it as a reference. The module exists because the library
is Float-only by design and manufacturing 50-digit eigenvalue references
requires arbitrary precision — mpmath does in one line what Lean
deliberately cannot.

Exactness channel: each reference eigenvalue is stored twice —
  * "eigvals_display": the 17-significant-digit decimal (for humans; 17
    digits uniquely determine any float64),
  * "eigvals": the exact dyadic pair [sig, exp] with value sig * 2^exp
    (via math.frexp; sig an exact integer below 2^53) — the bit-exact
    transfer channel that Golden.lean uses, immune to decimal-parser
    rounding.

Schema note: "cases" is written before "provenance" — the Lean mini-reader
parses the cases and ignores the rest. Run from the repo root
(`just golden`) so golden.json lands at the root.
"""

import json
import math

import mpmath as mp

import gen_cases

mp.mp.dps = 50

GOLDEN_JSON = "golden.json"


def golden_eigvals(A) -> list:
    """Eigenvalues of the *exact* float64 matrix at 50 digits, ascending."""
    a00, a11, a22, a01, a02, a12 = (mp.mpf(x) for x in A)
    M = mp.matrix(3, 3)
    M[0, 0], M[0, 1], M[0, 2] = a00, a01, a02
    M[1, 0], M[1, 1], M[1, 2] = a01, a11, a12
    M[2, 0], M[2, 1], M[2, 2] = a02, a12, a22
    E = mp.eigsy(M, eigvals_only=True)
    return sorted(float(x) for x in E)  # float() rounds correctly to float64


def dyadic(x: float) -> list:
    """Exact pair [sig, exp] with x == sig * 2^exp."""
    if x == 0.0:
        return [0, 0]
    m, e = math.frexp(x)
    return [int(m * 2 ** 53), e - 53]


def curated_cases() -> dict:
    cases = dict(gen_cases.ZOO)
    d = 1e-8
    cases["path-near-double"] = (-1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0)
    cases["path-near-triple"] = (1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0)
    cases["frontier-huge-id"] = (1e300, 1e300, 1e300, 0.0, 0.0, 0.0)
    cases["frontier-tiny-id"] = (1e-300, 1e-300, 1e-300, 0.0, 0.0, 0.0)
    return cases


def generate(path: str = GOLDEN_JSON) -> None:
    cases = []
    for name, A in curated_cases().items():
        vals = golden_eigvals(A)
        cases.append({
            "name": name,
            "matrix": [float(x) for x in A],
            "eigvals_display": [format(v, ".17g") for v in vals],
            "eigvals": [dyadic(v) for v in vals],
        })
        print(f"  {name:<18} " + "  ".join(f"{v:.17g}" for v in vals))
    doc = {
        "cases": cases,
        "provenance": {
            "generator": "parity/golden.py",
            "mpmath": mp.__version__,
            "dps": mp.mp.dps,
            "encoding": "eigvals are exact dyadic pairs [sig, exp]: "
                        "value = sig * 2^exp, the correctly-rounded float64 "
                        "of the 50-digit result.",
        },
    }
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
    print(f"wrote {path} ({len(cases)} cases)")


def main() -> int:
    generate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
