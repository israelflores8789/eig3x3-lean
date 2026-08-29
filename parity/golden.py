# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""golden.py — generates golden.json

This module is purely a generator of 50-digit eigenvalue "golden" reference
vectors using Python's mpmath, stored in generated/golden.json which has
two consumers:

    * The test Tests/Golden.lean consumes golden.json directly.
    * The Python module compare.py references golden.json for its
      parity checks.

Exactness channel: each reference eigenvalue is stored twice —
  * "eigvals_display": the 17-significant-digit decimal
    (human-readable 17-digit float64 representation),
  * "eigvals": the exact dyadic pair [sig, exp] with value sig * 2^exp
    (via math.frexp; sig an exact integer below 2^53) — the bit-exact
    transfer channel that Golden.lean uses, immune to decimal-parser
    rounding.

Schema note: "cases" is written before "provenance" — the Lean mini-reader
(JsonMini.lean) parses the cases and ignores the rest.

IMPORTANT: NEVER modify golden.json directly.

Usage:
    python golden.py                     # writes generated/golden.json
    just golden                          # equivalent
"""

import json
import math
import os

import mpmath as mp

import gen_cases

mp.mp.dps = 50

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLDEN_JSON = os.path.join(PROJECT_ROOT, "generated", "golden.json")


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
    return [int(m * 2**53), e - 53]


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
        cases.append(
            {
                "name": name,
                "matrix": [float(x) for x in A],
                "eigvals_display": [format(v, ".17g") for v in vals],
                "eigvals": [dyadic(v) for v in vals],
            }
        )
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
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
    print(f"wrote {path} ({len(cases)} cases)")


def main() -> int:
    generate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
