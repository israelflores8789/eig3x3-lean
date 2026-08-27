# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Parity against numpy.linalg.eigh.

Policy:
Eigenvalues are compared directly (well-posed, ordered on both sides —
LAPACK returns them ascending, so do we).

Eigenvectors are NEVER compared entrywise. They are sign-ambiguous, and
individually ill-posed inside clusters. The certificate metrics are the
vector-level evidence, and BOTH sides (this library and LAPACK) must meet
the same test gates.

Usage:
    python compare.py --impl mirror        # against the python replica
    python compare.py --impl lean          # against the Lean source
"""

import argparse
import sys

import numpy as np

import gen_cases
import lean_cli
import mirror
from gates import EPS, ORTHO_TOL, cert_tol, max_abs_entry

CERT_KEYS = ("maxResidual", "orthogonality", "reconstruction")


def to_dense(A) -> np.ndarray:
    a00, a11, a22, a01, a02, a12 = A
    return np.array([[a00, a01, a02], [a01, a11, a12], [a02, a12, a22]])


def numpy_decomposition(A) -> tuple:
    """LAPACK via numpy: w ascending, columns of V are the eigenvectors."""
    M = to_dense(A)
    w, V = np.linalg.eigh(M)
    res = max(float(np.max(np.abs(M @ V[:, i] - w[i] * V[:, i]))) for i in range(3))
    ort = float(np.max(np.abs(V.T @ V - np.eye(3))))
    rec = float(np.max(np.abs(V @ np.diag(w) @ V.T - M)))
    certs = {"maxResidual": res, "orthogonality": ort, "reconstruction": rec}
    return (
        tuple(float(x) for x in w),
        (
            tuple(float(x) for x in V[:, 0]),
            tuple(float(x) for x in V[:, 1]),
            tuple(float(x) for x in V[:, 2]),
        ),
        certs,
    )


def mirror_batch(cases: list) -> list:
    return [mirror.decomposition(A) for A in cases]


def lean_batch(binary: str):
    def run_batch(cases: list) -> list:
        return lean_cli.run(cases, binary)

    return run_batch


def compare_suite(name: str, cases: list, impl_batch) -> dict:
    refs = [numpy_decomposition(A) for A in cases]
    outs = impl_batch(cases)
    worst_eval = 0.0
    worst_cert = {k: 0.0 for k in CERT_KEYS}
    failures = []
    for A, (e_ref, _, c_ref), (e, _, c) in zip(cases, refs, outs, strict=True):
        mA = max_abs_entry(A)
        unit = EPS * mA
        err = max(abs(e[i] - e_ref[i]) for i in range(3))
        # The zero matrix has zero scale; both sides are exact there anyway.
        worst_eval = max(worst_eval, err / unit if unit > 0.0 else err)
        for key in CERT_KEYS:
            gate = ORTHO_TOL if key == "orthogonality" else cert_tol(A)
            for certs, tag in ((c, "impl"), (c_ref, "numpy")):
                value = certs[key]
                if value > gate:
                    failures.append((name, key, tag, value, gate))
                if tag == "impl" and gate > 0.0:
                    worst_cert[key] = max(worst_cert[key], value / gate)
    return {
        "suite": name,
        "n": len(cases),
        "worst_eval_units": worst_eval,
        "worst_cert_frac": worst_cert,
        "failures": failures,
    }


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--impl", choices=["mirror", "lean"], default="mirror")
    p.add_argument("--lean-binary", default=lean_cli.DEFAULT_BINARY)
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

    if args.impl == "lean" and not lean_cli.available(args.lean_binary):
        p.error(
            f"Lean binary not found at {args.lean_binary!r}; "
            "build eig3x3_cli first (the CLI task)"
        )
    impl_batch = mirror_batch if args.impl == "mirror" else lean_batch(args.lean_binary)

    rng = np.random.default_rng(args.seed)
    print(
        f"# eig3x3 parity — impl={args.impl} seed={args.seed} "
        f"(errors in eps*maxAbs units; budget is 64)"
    )
    total_fail = 0
    for suite in args.suites:
        cases = gen_cases.cases_for(suite, args.n, rng)
        rep = compare_suite(suite, cases, impl_batch)
        wc = rep["worst_cert_frac"]
        print(
            f"{suite:<11} n={rep['n']:<6} worst eval err "
            f"{rep['worst_eval_units']:7.2f} eps*maxAbs | "
            f"worst cert fractions of gate: residual {wc['maxResidual']:.3f}, "
            f"ortho {wc['orthogonality']:.3f}, recon {wc['reconstruction']:.3f}"
        )
        for f in rep["failures"][:5]:
            print(f"    FAIL {f}")
        total_fail += len(rep["failures"])
    if total_fail:
        print(f"parity: {total_fail} gate violation(s)")
        return 1
    print("parity: all gates held")
    return 0


if __name__ == "__main__":
    sys.exit(main())
