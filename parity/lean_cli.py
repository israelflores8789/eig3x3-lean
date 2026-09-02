# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""lean_cli.py — subprocess bridge to the Lean eig3x3_cli binary

This module owns the JSON schema (the contract Cli.lean implements).
One process invocation per batch, never per matrix.

Schema
------
Input (stdin, JSON):
    {"matrices": [[a00, a11, a22, a01, a02, a12], ...]}
The six entries are Eig3x3.SymmMat3 field order.

Output (stdout, JSON):
    {"results": [{"matrix": [a00, ..., a12],
                  "eigvals": [l1, l2, l3],
                  "eigvecs": [c1x, c1y, c1z,
                              c2x, c2y, c2z,
                              c3x, c3y, c3z],
                  "certificates": {"maxResidual": r,
                                   "orthogonality": o,
                                   "reconstruction": rc}},
                 ...]}
eigvecs are COLUMN-MAJOR: column i is the unit eigenvector for l_i (consistent
with NumPy / SciPy / PyTorch convention). Getting this order wrong can cause a
silent failure. The schema is explicit about this on purpose.

Boundary verification (two distinct bitwise checks):
  * input fidelity: the echoed "matrix" is the CLI's parsed view of the
    input; comparing it against what was sent measures the input parser.
    Informational — a ≤1–2 ulp gap is the documented OfScientific caveat,
    far below the parity gates.
  * output self-consistency: the payload's certificate fields are recomputed
    from the payload's own eigenpairs and echoed matrix (see
    `recompute_certificates`) and must match BITWISE. Certificates use only
    +, −, ×, ÷ — correctly rounded to 1/2 ulp and deterministic under IEEE 754 —
    so any corruption of the payload (transposition, field misorder, a 1-ulp
    slip) breaks the equality. This is an op-order mirror of Certificates.lean.
"""

import json
import os
import subprocess
from typing import NamedTuple

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_BINARY = os.path.join(PROJECT_ROOT, ".lake", "build", "bin", "eig3x3_cli")

CERT_KEYS = ("maxResidual", "orthogonality", "reconstruction")


class Result(NamedTuple):
    """One decomposition payload from the CLI."""

    sent: tuple  # the matrix as handed to the CLI
    echo: tuple  # the CLI's parsed view of it
    eigvals: tuple  # (l1, l2, l3)
    eigvecs: tuple  # (c1, c2, c3), columns
    certs: dict  # CERT_KEYS


def available(binary: str = DEFAULT_BINARY) -> bool:
    return os.path.isfile(binary) and os.access(binary, os.X_OK)


def run(matrices: list, binary: str = DEFAULT_BINARY) -> list:
    """Run the CLI on a batch of 6-tuples; returns a list of Result."""
    if not available(binary):
        raise FileNotFoundError(
            f"eig3x3_cli not found at {binary!r}; run `just build-cli` first"
        )
    payload = {"matrices": [[float(x) for x in A] for A in matrices]}
    proc = subprocess.run(
        [binary], input=json.dumps(payload), text=True, capture_output=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"eig3x3_cli exited {proc.returncode}: {proc.stderr.strip()}"
        )
    data = json.loads(proc.stdout)
    out = []
    for A, r in zip(matrices, data["results"], strict=True):
        if "matrix" not in r:
            raise RuntimeError(
                "eig3x3_cli result is missing the matrix echo field; "
                "try rebuilding with `just build-cli`.\n"
                "If the problem persists, check that Cli.lean dispatches "
                'an echoed "matrix".'
            )
        echo = tuple(float(x) for x in r["matrix"])
        e = tuple(float(x) for x in r["eigvals"])
        v = [float(x) for x in r["eigvecs"]]
        Q = (tuple(v[0:3]), tuple(v[3:6]), tuple(v[6:9]))
        certs = {k: float(val) for k, val in r["certificates"].items()}
        out.append(Result(A, echo, e, Q, certs))
    return out


def decomposition(A, binary: str = DEFAULT_BINARY) -> Result:
    """Single-matrix convenience wrapper."""
    return run([A], binary)[0]


# --- Certificate recomputation: op-order mirror of Certificates.lean -------
# Pure +, -, x over the payload's own numbers, so IEEE 754 determinism makes
# the recomputation bitwise-identical to Lean's `certify` whenever the
# payload is faithful.


def _mul_vec(A, v):
    a00, a11, a22, a01, a02, a12 = A
    return (
        a00 * v[0] + a01 * v[1] + a02 * v[2],
        a01 * v[0] + a11 * v[1] + a12 * v[2],
        a02 * v[0] + a12 * v[1] + a22 * v[2],
    )


def _dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _residual(A, lam, v):
    av = _mul_vec(A, v)
    r = (av[0] - lam * v[0], av[1] - lam * v[1], av[2] - lam * v[2])
    return max(abs(r[0]), max(abs(r[1]), abs(r[2])))


def _orthogonality_error(Q):
    c1, c2, c3 = Q
    m1 = max(abs(_dot(c1, c1) - 1.0), abs(_dot(c2, c2) - 1.0))
    m2 = max(abs(_dot(c3, c3) - 1.0), abs(_dot(c1, c2)))
    m3 = max(abs(_dot(c1, c3)), abs(_dot(c2, c3)))
    return max(m1, max(m2, m3))


def _reconstruct(e, Q):
    def comp(m, v):
        return (
            m * v[0] * v[0],
            m * v[1] * v[1],
            m * v[2] * v[2],
            m * v[0] * v[1],
            m * v[0] * v[2],
            m * v[1] * v[2],
        )

    b1, b2, b3 = comp(e[0], Q[0]), comp(e[1], Q[1]), comp(e[2], Q[2])
    return tuple(b1[i] + b2[i] + b3[i] for i in range(6))


def _reconstruction_error(A, e, Q):
    B = _reconstruct(e, Q)
    m1 = max(abs(B[0] - A[0]), abs(B[1] - A[1]))
    m2 = max(abs(B[2] - A[2]), abs(B[3] - A[3]))
    m3 = max(abs(B[4] - A[4]), abs(B[5] - A[5]))
    return max(m1, max(m2, m3))


def recompute_certificates(A, e, Q) -> dict:
    """Certificates of (e, Q) against A, exactly as Certificates.lean."""
    return {
        "maxResidual": max(_residual(A, e[i], Q[i]) for i in range(3)),
        "orthogonality": _orthogonality_error(Q),
        "reconstruction": _reconstruction_error(A, e, Q),
    }


def self_consistent(res: Result) -> bool | None:
    """Bitwise self-consistency certificate.

    The payload's certificate fields must equal the certificates recomputed
    from its own eigenpairs and echoed matrix. Works because the certificates
    use only +, −, ×, ÷, which are deterministic under IEEE 754 — any
    serialization corruption (transposition, field misorder, a 1-ulp slip)
    breaks the equality.

    The echo field is mandatory; `run()` rejects binaries that lack it.
    """
    return recompute_certificates(res.echo, res.eigvals, res.eigvecs) == res.certs
