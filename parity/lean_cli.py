# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Subprocess bridge to the eig3x3-cli Lean binary.

This module defines the contract the Lean side implements (the spec for the
CLI task). One process invocation per batch, never per matrix.

Schema
------
Input (stdin, JSON):
    {"matrices": [[a00, a11, a22, a01, a02, a12], ...]}
The six entries are Eig3x3.SymmMat3 field order.

Output (stdout, JSON):
    {"results": [{"eigvals": [l1, l2, l3],
                  "eigvecs": [c1x, c1y, c1z,
                              c2x, c2y, c2z,
                              c3x, c3y, c3z],
                  "certificates": {"maxResidual": r,
                                   "orthogonality": o,
                                   "reconstruction": rc}},
                 ...]}
eigvecs are COLUMN-MAJOR: c_i is the unit eigenvector for l_i (the NumPy /
SciPy / PyTorch convention). Getting this order wrong is the classic silent
failure at this boundary — the schema says it twice on purpose.

Serialization note: output floats must round-trip exactly. Lean's default
Float display prints six significant digits, which is lossy — the CLI must
emit full precision (e.g. via the exact dyadic decimal expansion).
"""

import json
import os
import subprocess

DEFAULT_BINARY = os.path.join(".lake", "build", "bin", "eig3x3-cli")


def available(binary: str = DEFAULT_BINARY) -> bool:
    return os.path.isfile(binary) and os.access(binary, os.X_OK)


def run(matrices: list, binary: str = DEFAULT_BINARY) -> list:
    """Run the CLI on a batch; returns [(eigvals, eigvec columns, certs)]."""
    if not available(binary):
        raise FileNotFoundError(
            f"eig3x3-cli not found at {binary!r}; build it with "
            "`lake build eig3x3-cli` (the CLI task)")
    payload = {"matrices": [[float(x) for x in A] for A in matrices]}
    proc = subprocess.run([binary], input=json.dumps(payload),
                          text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(f"eig3x3-cli exited {proc.returncode}: "
                           f"{proc.stderr.strip()}")
    data = json.loads(proc.stdout)
    out = []
    for r in data["results"]:
        e = tuple(float(x) for x in r["eigvals"])
        v = [float(x) for x in r["eigvecs"]]
        Q = (tuple(v[0:3]), tuple(v[3:6]), tuple(v[6:9]))
        certs = {k: float(val) for k, val in r["certificates"].items()}
        out.append((e, Q, certs))
    return out


def decomposition(A, binary: str = DEFAULT_BINARY):
    """Uniform harness interface for a single matrix."""
    return run([A], binary)[0]
