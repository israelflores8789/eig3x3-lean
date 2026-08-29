# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""parity.py — parity test suite

This module takes the results of running eig3x3_cli and performs the parity
comparison against the results of running numpy/LAPACK.

Doctrine:

  * Eigenvalues are compared directly — well-posed, ascending on both
    sides (LAPACK's convention and the Lean implementation).
  * Eigenvectors are NEVER compared entrywise across implementations
    since eigenvectors are sign- and cluster-ambiguous. Instead, the
    certificate metrics are the vector-level evidence, and BOTH sides —
    this library and LAPACK — MUST meet the same test gates.
  * Every CLI payload is checked for bitwise certificate self-consistency
    (recomputed from its own eigenpairs; see lean_cli.py), which acts as
    a serialization guard.
  * The "golden" test suite compares against golden.json (50-digit mpmath
    references), reconstructed bit-exactly from dyadic pairs.

Usage:
    python parity.py                    # after `just build_cli`
"""

import argparse
import json
import math
import os
import sys

import numpy as np
import pytest

import gen_cases
import lean_cli
from gates import EPS, ORTHO_TOL, cert_tol, max_abs_entry

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_GOLDEN = os.path.join(PROJECT_ROOT, "generated", "golden.json")

CERT_KEYS = lean_cli.CERT_KEYS


def to_dense(A) -> np.ndarray:
    a00, a11, a22, a01, a02, a12 = A
    return np.array([[a00, a01, a02], [a01, a11, a12], [a02, a12, a22]])


def numpy_reference(A) -> tuple:
    """LAPACK via numpy: eigh (NOT eigvalsh — we need V for LAPACK-side
    certificates). w ascending; columns of V are the eigenvectors."""
    M = to_dense(A)
    w, V = np.linalg.eigh(M)
    res = max(float(np.max(np.abs(M @ V[:, i] - w[i] * V[:, i]))) for i in range(3))
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
        refs.append(
            (
                tuple(math.ldexp(int(sig), int(exp)) for sig, exp in case["eigvals"]),
                None,
            )
        )
    return cases, refs


def ulp_gap(a: float, b: float) -> float:
    if a == b:
        return 0.0
    u = math.ulp(a) if a != 0.0 else (math.ulp(b) if b != 0.0 else 1.0)
    return abs(a - b) / u


def compare_suite(name: str, cases: list, refs: list, results: list) -> dict:
    """refs: [(ref_eigvals, ref_certs_or_None)]; results: lean_cli.Result."""
    worst_eval = 0.0
    worst_cert = {k: 0.0 for k in CERT_KEYS}
    failures: list[tuple] = []
    echo_mismatches = 0
    worst_ulp = 0.0
    for A, (e_ref, c_ref), res in zip(cases, refs, results, strict=True):
        mA = max_abs_entry(A)
        unit = EPS * mA
        err = max(abs(res.eigvals[i] - e_ref[i]) for i in range(3))
        # The zero matrix has zero scale; both sides are exact there anyway.
        worst_eval = max(worst_eval, err / unit if unit > 0.0 else err)
        tol = cert_tol(A)

        # Boundary checks: input fidelity (informational) and output
        # self-consistency (bitwise, gated).
        gap = max(ulp_gap(a, b) for a, b in zip(A, res.echo, strict=True))
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

    return {
        "suite": name,
        "n": len(cases),
        "worst_eval_units": worst_eval,
        "worst_cert_frac": worst_cert,
        "failures": failures,
        "echo_mismatches": echo_mismatches,
        "worst_ulp": worst_ulp,
    }


def print_report(rep: dict) -> None:
    wc = rep["worst_cert_frac"]
    print(
        f"{rep['suite']:<11} n={rep['n']:<6} worst eval err "
        f"{rep['worst_eval_units']:7.2f} eps*maxAbs (budget 64) | "
        f"cert fractions: residual {wc['maxResidual']:.3f}, "
        f"ortho {wc['orthogonality']:.3f}, recon {wc['reconstruction']:.3f}"
    )
    if rep["echo_mismatches"] or rep["worst_ulp"] > 0.0:
        print(
            f"    input fidelity: {rep['echo_mismatches']} echoed entries "
            f"differ (worst {rep['worst_ulp']:.1f} ulp; ≤2 ulp is the "
            f"documented parse caveat)"
        )
    for f in rep["failures"][:5]:
        print(f"    FAIL {f}")


# ---------------------------------------------------------------------------
# Pytest Test Entry Points (Granular Suite Execution)
# ---------------------------------------------------------------------------


@pytest.mark.zoo
def test_zoo(lean_cli_binary: str) -> None:
    """Parity against NumPy/LAPACK on the curated ZOO matrices."""
    cases = gen_cases.cases_for("zoo", 0, np.random.default_rng(0))
    refs = [numpy_reference(A) for A in cases]
    results = lean_cli.run(cases, lean_cli_binary)
    rep = compare_suite("zoo", cases, refs, results)
    assert not rep["failures"], f"zoo failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


@pytest.mark.random
def test_random(lean_cli_binary: str, batch_n: int, rng: np.random.Generator) -> None:
    """Parity against NumPy/LAPACK on uniform random matrices in [-1, 1]."""
    cases = gen_cases.cases_for("random", batch_n, rng)
    refs = [numpy_reference(A) for A in cases]
    results = lean_cli.run(cases, lean_cli_binary)
    rep = compare_suite("random", cases, refs, results)
    assert not rep["failures"], f"random failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


@pytest.mark.logscale
def test_logscale(lean_cli_binary: str, batch_n: int, rng: np.random.Generator) -> None:
    """Parity across 600 orders of magnitude (1e-300 to 1e300)."""
    cases = gen_cases.cases_for("logscale", batch_n, rng)
    refs = [numpy_reference(A) for A in cases]
    results = lean_cli.run(cases, lean_cli_binary)
    rep = compare_suite("logscale", cases, refs, results)
    assert not rep["failures"], f"logscale failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


@pytest.mark.paths
def test_paths(lean_cli_binary: str) -> None:
    """Parity along Habera-Zilian adversarial perturbation paths."""
    cases = gen_cases.cases_for("paths", 0, np.random.default_rng(0))
    refs = [numpy_reference(A) for A in cases]
    results = lean_cli.run(cases, lean_cli_binary)
    rep = compare_suite("paths", cases, refs, results)
    assert not rep["failures"], f"paths failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


@pytest.mark.smallscale
def test_smallscale(lean_cli_binary: str) -> None:
    """Parity on vanishing double eigenvalues (s, s, 2s) at 1e-300 to 1e-1."""
    cases = gen_cases.cases_for("smallscale", 0, np.random.default_rng(0))
    refs = [numpy_reference(A) for A in cases]
    results = lean_cli.run(cases, lean_cli_binary)
    rep = compare_suite("smallscale", cases, refs, results)
    assert not rep["failures"], f"smallscale failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


@pytest.mark.frontier
def test_frontier(lean_cli_binary: str) -> None:
    """Parity on frontier scale-edge cases (2^±997, extreme dynamic range)."""
    cases = gen_cases.cases_for("frontier", 0, np.random.default_rng(0))
    refs = [numpy_reference(A) for A in cases]
    results = lean_cli.run(cases, lean_cli_binary)
    rep = compare_suite("frontier", cases, refs, results)
    assert not rep["failures"], f"frontier failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


@pytest.mark.golden
def test_golden(lean_cli_binary: str, golden_json: str) -> None:
    """Parity against 50-digit mpmath reference vectors from golden.json."""
    gcases, grefs = load_golden(golden_json)
    results = lean_cli.run(gcases, lean_cli_binary)
    rep = compare_suite("golden", gcases, grefs, results)
    assert not rep["failures"], f"golden failures: {rep['failures']}"
    assert rep["worst_eval_units"] <= 64.0


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--lean-binary", default=lean_cli.DEFAULT_BINARY)
    p.add_argument("--golden", default=DEFAULT_GOLDEN)
    p.add_argument(
        "--n", type=int, default=20000, help="batch size for the random suites"
    )
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--suites",
        nargs="+",
        default=["zoo", "random", "logscale", "paths", "smallscale", "frontier"],
    )
    args = p.parse_args()

    if not lean_cli.available(args.lean_binary):
        p.error(
            f"eig3x3_cli not found at {args.lean_binary!r}; run `just build_cli` first"
        )

    rng = np.random.default_rng(args.seed)
    print(
        f"# eig3x3 parity — seed={args.seed} "
        f"(eval errors in eps*maxAbs units; budget is 64)"
    )
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
        print(
            f"# golden skipped: {args.golden} not found "
            f"(run `just golden` from the repo root)"
        )
        gcases, grefs = None, None
    if gcases and grefs:
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
