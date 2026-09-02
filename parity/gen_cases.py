# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""gen_cases.py — case generation

The zoo, random matrices across scales, Habera-Zilian's proposed adversarial
paths, and frontier cases.

Everything is seeded; the seed travels in every report header. Matrices are
6-tuples (a00, a11, a22, a01, a02, a12), matching Lean's Eig3x3.SymmMat3.
"""

import numpy as np

# The same six matrices as Tests/Util.lean. Keep in sync.
ZOO = {
    "worked": (2.0, 2.0, 2.0, 1.0, 0.0, 1.0),  # eigenvalues 2±√2, 2
    "diagonal": (1.0, 2.0, 3.0, 0.0, 0.0, 0.0),
    "scaled-id": (3.7, 3.7, 3.7, 0.0, 0.0, 0.0),
    "zero": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    "near-double": (1.0, 1.0, -1.0, 1.0e-8, 0.0, 0.0),
    "near-triple": (1.0, 1.0, 1.0 + 1.0e-8, 0.0, 0.0, 0.0),
    "regression": (-2.0, 2.0, -9.0, 9.0, 8.0, 6.0),  # the r10 matrix
}

# Log-spaced perturbation sizes for the adversarial paths.
DELTAS = [float(d) for d in np.logspace(-16, -4, 25)]


def random_uniform(rng: np.random.Generator, n: int) -> list:
    """Entries uniform in [-1, 1] — the bread-and-butter parity batch."""
    return [tuple(float(x) for x in rng.uniform(-1.0, 1.0, 6)) for _ in range(n)]


def random_logscale(
    rng: np.random.Generator, n: int, emin: float = -300.0, emax: float = 300.0
) -> list:
    """Entries uniform in [-1, 1] times 10^k, k log-uniform in [emin, emax].

    Scale coverage is what exercises the preconditioner: 600 orders of
    magnitude, up to the edge of the documented frontier.
    """
    out = []
    for _ in range(n):
        scale = 10.0 ** float(rng.uniform(emin, emax))
        out.append(tuple(float(x) * scale for x in rng.uniform(-1.0, 1.0, 6)))
    return out


def adversarial_paths() -> list:
    """The paper's benchmark paths (H-Z 2025): eigenvalue clusters driven by
    a log-spaced perturbation delta.

    * diag(-1, 1, 1+d):  near-double at the top of the spectrum
    * diag(1, 1, 1+d):   near-triple
    * [[1, d], [d, 1]] + diag(-1): off-diagonal cluster (the D2-style path
      that exposes the naive discriminant's cancellation)
    """
    cases = []
    for d in DELTAS:
        cases.append((-1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0))
        cases.append((1.0, 1.0, 1.0 + d, 0.0, 0.0, 0.0))
        cases.append((1.0, 1.0, -1.0, d, 0.0, 0.0))
    return cases


def double_at_small_scale() -> list:
    """Exact double eigenvalue (s, s, 2s) at vanishing scale — the boundary
    adversarial family from the validation report."""
    return [
        tuple(float(s) * x for x in (1.0, 1.0, 2.0, 0.0, 0.0, 0.0))
        for s in np.logspace(-300, -1, 25)
    ]


def frontier_cases() -> list:
    """The edge of the working range. The scale factors here (2^±997-ish)
    are representable, so all of these must succeed bit-cleanly."""
    return [
        (1e300, 1e300, 1e300, 0.0, 0.0, 0.0),  # scaled identity, huge
        (1e-300, 1e-300, 1e-300, 0.0, 0.0, 0.0),  # scaled identity, tiny
        (1e300, 1.0, 1e-300, 0.0, 0.0, 0.0),  # extreme dynamic range
    ]


def integer_matrices(
    rng: np.random.Generator, n: int, lo: int = -9, hi: int = 9
) -> list:
    """Small integer-valued matrices — the inputs for exact.py's rational
    identity checks."""
    return [tuple(int(x) for x in rng.integers(lo, hi + 1, 6)) for _ in range(n)]


def cases_for(suite: str, n: int, rng: np.random.Generator) -> list:
    if suite == "zoo":
        return list(ZOO.values())
    if suite == "random":
        return random_uniform(rng, n)
    if suite == "logscale":
        return random_logscale(rng, n)
    if suite == "paths":
        return adversarial_paths()
    if suite == "smallscale":
        return double_at_small_scale()
    if suite == "frontier":
        return frontier_cases()
    raise ValueError(f"unknown suite: {suite}")
