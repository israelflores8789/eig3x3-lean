"""compare.py — the parity suite: eig3x3_cli vs numpy/LAPACK.

Doctrine: Python's job is parity against references, full stop.

  * Eigenvalues are compared directly — well-posed, ascending on both
    sides (LAPACK's convention and ours).
  * Eigenvectors are NEVER compared entrywise across implementations
    (sign- and cluster-ambiguous). The certificate metrics are the
    vector-level evidence, and BOTH sides — this library and LAPACK — must
    meet the same gates.
  * Every CLI payload is checked for bitwise certificate self-consistency
    (recomputed from its own eigenpairs; see lean_cli.py), which is the
    serialization tripwire for the boundary.
  * The golden suite compares against golden.json (50-digit mpmath
    references), reconstructed bit-exactly from dyadic pairs.

Usage: python parity/compare.py            (repo root, after `just cli`)
"""

import argparse
import json
import math
import sys

import numpy as np

from gates import EPS, ORTHO_TOL, cert_tol, max_abs_entry
import gen_cases
import lean_cli

CERT_KEYS = lean_cli.CERT_KEYS


def to_dense(A) -> np.ndarray:
    a00, a11, a22, a01, a02, a12 = A
    return np.array([[a00, a01, a02],
                     [a01, a11, a12],
                     [a02, a12, a22]])


def numpy_reference(A) -> tuple:
    """LAPACK via numpy: eigh (NOT eigvalsh — we need V for LAPACK-side
    certificates). w ascending; columns of V are the eigenvectors."""
    M = to_dense(A)
    w, V = np.linalg.eigh(M)
    res = max(float(np.max(np.abs(M @ V[:, i] - w[i] * V[:, i])))
              for i in range(3))
    ort = float(np.max(np.abs(V.T @ V - np.eye(3))))
    rec = float(np.max(np.abs(V @ np.diag(w) @ V.T - M)))
    certs = {"maxResidual": res, "orthogonality": ort, "reconstruction": rec}
    return (tuple(float(x) for x in w), certs)


def load_golden(path: str) -> tuple:
    """golden.json -> (cases, refs); refs reconstructed bit-exactly."""
    with open(path) as f:
        doc = json.load(f)
    cases, refs = [], []
    for case in doc["cases"]:
        cases.append(tuple(float(x) for x in case["matrix"]))
        refs.append((tuple(math.ldexp(int(sig), int(exp))
                           for sig, exp in case["eigvals"]), None))
    return cases, refs


def ulp_gap(a: float, b: float) -> float:
    if a == b:
        return 0.0
    u = math.ulp(a) if a != 0.0 else (math.ulp(b) if b != 0.0 else 1.0)
    return abs(a - b) / u


def compare_suite(name: str, cases: list, refs: list,
                  results: list) -> dict:
    """refs: [(ref_eigvals, ref_certs_or_None)]; results: lean_cli.Result."""
    worst_eval = 0.0
    worst_cert = {k: 0.0 for k in CERT_KEYS}
    failures = []
    echo_mismatches = 0
    worst_ulp = 0.0
    consistency_skipped = 0
    for A, (e_ref, c_ref), res in zip(cases, refs, results):
        mA = max_abs_entry(A)
        unit = EPS * mA
        err = max(abs(res.eigvals[i] - e_ref[i]) for i in range(3))
        # The zero matrix has zero scale; both sides are exact there anyway.
        worst_eval = max(worst_eval, err / unit if unit > 0.0 else err)
        tol = cert_tol(A)

        # Boundary checks: input fidelity (informational) and output
        # self-consistency (bitwise, gated).
        if res.echo is None:
            consistency_skipped += 1
        else:
            gap = max(ulp_gap(a, b) for a, b in zip(A, res.echo))
            if gap > 0.0:
                echo_mismatches += 1
                worst_ulp = max(worst_ulp, gap)
            if not lean_cli.self_consistent(res):
                failures.append((name, "self-consistency", A))

        # Certificates under the shared gates — both implementations.
        cert_pairs = [(res.certs, "cli")]
        if c_ref is not None:
            cert_pairs.append((c_ref, "numpy"))
        for certs, tag in cert_pairs:
            for key in CERT_KEYS:
                gate = ORTHO_TOL if key == "orthogonality" else tol
                if certs[key] > gate:
                    failures.append((name, key, tag, certs[key], gate))
                if tag == "cli" and gate > 0.0:
                    worst_cert[key] = max(worst_cert[key], certs[key] / gate)

    return {"suite": name, "n": len(cases), "worst_eval_units": worst_eval,
            "worst_cert_frac": worst_cert, "failures": failures,
            "echo_mismatches": echo_mismatches, "worst_ulp": worst_ulp,
            "consistency_skipped": consistency_skipped}


def print_report(rep: dict) -> None:
    wc = rep["worst_cert_frac"]
    print(f"{rep['suite']:<11} n={rep['n']:<6} worst eval err "
          f"{rep['worst_eval_units']:7.2f} eps*maxAbs (budget 64) | "
          f"cert fractions: residual {wc['maxResidual']:.3f}, "
          f"ortho {wc['orthogonality']:.3f}, recon {wc['reconstruction']:.3f}")
    if rep["echo_mismatches"] or rep["worst_ulp"] > 0.0:
        print(f"    input fidelity: {rep['echo_mismatches']} echoed entries "
              f"differ (worst {rep['worst_ulp']:.1f} ulp; ≤2 ulp is the "
              f"documented parse caveat)")
    if rep["consistency_skipped"]:
        print(f"    note: self-consistency skipped for "
              f"{rep['consistency_skipped']} result(s) without a matrix "
              f"echo — apply the echo patch to tests/Cli.lean")
    for f in rep["failures"][:5]:
        print(f"    FAIL {f}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lean-binary", default=lean_cli.DEFAULT_BINARY)
    p.add_argument("--golden", default="golden.json")
    p.add_argument("--n", type=int, default=20000,
                   help="batch size for the random suites")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--suites", nargs="+",
                   default=["zoo", "random", "logscale", "paths",
                            "smallscale", "frontier"])
    args = p.parse_args()

    if not lean_cli.available(args.lean_binary):
        p.error(f"eig3x3_cli not found at {args.lean_binary!r}; "
                "run `just cli` first")

    rng = np.random.default_rng(args.seed)
    print(f"# eig3x3 parity — seed={args.seed} "
          f"(eval errors in eps*maxAbs units; budget is 64)")
    total_fail = 0

    for suite in args.suites:
        cases = gen_cases.cases_for(suite, args.n, rng)
        refs = [numpy_reference(A) for A in cases]
        results = lean_cli.run(cases, args.lean_binary)
        rep = compare_suite(suite, cases, refs, results)
        print_report(rep)
        total_fail += len(rep["failures"])

    try:
        gcases, grefs = load_golden(args.golden)
    except FileNotFoundError:
        print(f"# golden skipped: {args.golden} not found "
              f"(run `just golden` from the repo root)")
        gcases = None
    if gcases:
        results = lean_cli.run(gcases, args.lean_binary)
        rep = compare_suite("golden", gcases, grefs, results)
        print_report(rep)
        total_fail += len(rep["failures"])

    if total_fail:
        print(f"parity: {total_fail} gate violation(s)")
        return 1
    print("parity: all gates held")
    return 0


if __name__ == "__main__":
    sys.exit(main())
