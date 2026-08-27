# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Golden vectors: 50-digit references via mpmath.

A small curated set (the zoo plus a few nasties) is solved at 50-digit
precision and emitted as 17-significant-digit decimals, which round-trip
exactly to the correct float64. Two consumers:

  * this harness (--check) validates an implementation against them;
  * the Golden.lean task checks a subset into Tests/Golden.lean, giving the
    Lean suite high-precision anchors with no Python at test time.

Usage:
    python golden.py                # write golden.json
    python golden.py --check        # validate the mirror against golden.json
"""

import json
import sys

import mpmath as mp

from gates import cert_tol
import gen_cases
import mirror

mp.mp.dps = 50

GOLDEN_JSON = "golden.json"


def golden_eigvals(A) -> list:
    """Eigenvalues of the *exact* float64 matrix, at 50 digits, ascending."""
    a00, a11, a22, a01, a02, a12 = (mp.mpf(x) for x in A)
    M = mp.matrix(3, 3)
    M[0, 0], M[0, 1], M[0, 2] = a00, a01, a02
    M[1, 0], M[1, 1], M[1, 2] = a01, a11, a12
    M[2, 0], M[2, 1], M[2, 2] = a02, a12, a22
    E = mp.eigsy(M, eigvals_only=True)
    return sorted(float(x) for x in E)  # float() rounds correctly to float64


def curated_cases() -> dict:
    cases = dict(gen_cases.ZOO)
    d = 1e-8
    cases["path-near-double"] = (-1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0)
    cases["path-near-triple"] = (1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0)
    cases["frontier-huge-id"] = (1e300, 1e300, 1e300, 0.0, 0.0, 0.0)
    cases["frontier-tiny-id"] = (1e-300, 1e-300, 1e-300, 0.0, 0.0, 0.0)
    return cases


def generate(path: str = GOLDEN_JSON) -> None:
    out = []
    for name, A in curated_cases().items():
        vals = golden_eigvals(A)
        out.append({
            "name": name,
            "matrix": [float(x) for x in A],
            # 17 significant digits round-trip exactly to the same float64.
            "eigvals": [format(v, ".17g") for v in vals],
        })
        print(f"  {name:<18} {vals[0]:.17g}  {vals[1]:.17g}  {vals[2]:.17g}")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"wrote {path} ({len(out)} cases)")


def check(path: str = GOLDEN_JSON) -> int:
    with open(path) as f:
        cases = json.load(f)
    failures = 0
    for case in cases:
        A = tuple(case["matrix"])
        ref = [float(s) for s in case["eigvals"]]
        e, _, _ = mirror.decomposition(A)
        err = max(abs(e[i] - ref[i]) for i in range(3))
        tol = cert_tol(A)
        status = "ok" if err <= tol else "FAIL"
        if err > tol:
            failures += 1
        print(f"  {status:<4} {case['name']:<18} max|err| {err:.3e} "
              f"(gate {tol:.3e})")
    print(f"golden: {failures} failure(s)")
    return 1 if failures else 0


def main() -> int:
    if "--check" in sys.argv:
        return check()
    generate()
    return 0


if __name__ == "__main__":
    sys.exit(main())
