# Copyright (c) 2026 Israel Flores-Arbolay. All rights reserved.
# Released under Apache 2.0 license as described in the file LICENSE.
# Authors: Israel Flores-Arbolay

"""Op-for-op float64 mirror of the Eig3x3 Lean library.

Ports Eigendecomp.lean / Eigenvalues.lean / Eigenvectors.lean /
Certificates.lean line for line, in the same operation order, so on the same
platform its results should be bit-identical to the Lean binary. The harness
uses it in two roles:

  1. stand-in for the CLI until eig3x3-cli exists (run everything today);
  2. cross-check thereafter: mirror and CLI must agree bitwise, and any
     disagreement is a serialization bug, not numerics.

`delta_naive` lives here (not in the library) for the naive-vs-present
benchmark study in bench.py.
"""

import math

TWO_PI_OVER_3 = 2.0943951023931953

Symm = tuple   # (a00, a11, a22, a01, a02, a12)
Vec3 = tuple   # (x, y, z)
Mat3 = tuple   # (c1, c2, c3) — columns

# ---- Vec3 / Mat3 operations (Basic.lean) ----

def cross(u: Vec3, v: Vec3) -> Vec3:
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])

def dot(u: Vec3, v: Vec3) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

def norm_sq(u: Vec3) -> float:
    return u[0] * u[0] + u[1] * u[1] + u[2] * u[2]

def vscale(v: Vec3, s: float) -> Vec3:
    return (s * v[0], s * v[1], s * v[2])

def vsub(u: Vec3, v: Vec3) -> Vec3:
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])

def mul_vec(A: Symm, v: Vec3) -> Vec3:
    a00, a11, a22, a01, a02, a12 = A
    return (a00 * v[0] + a01 * v[1] + a02 * v[2],
            a01 * v[0] + a11 * v[1] + a12 * v[2],
            a02 * v[0] + a12 * v[1] + a22 * v[2])

def max_abs_entry(A: Symm) -> float:
    a00, a11, a22, a01, a02, a12 = A
    return max(max(abs(a00), abs(a01)),
               max(max(abs(a02), abs(a11)), max(abs(a12), abs(a22))))

def mat_det(Q: Mat3) -> float:
    return dot(cross(Q[0], Q[1]), Q[2])

# ---- Eigenvalues (Eigenvalues.lean) ----

def i1(A: Symm) -> float:
    return A[0] + A[1] + A[2]

def j2(A: Symm) -> float:
    a00, a11, a22, a01, a02, a12 = A
    d0, d1, d2 = a00 - a11, a00 - a22, a11 - a22
    return ((d0 * d0 + d1 * d1 + d2 * d2) / 6.0
            + (a01 * a01 + a02 * a02 + a12 * a12))

def j3(A: Symm) -> float:
    a00, a11, a22, a01, a02, a12 = A
    d0, d1, d2 = a00 - a11, a00 - a22, a11 - a22
    t1, t2, t3 = d1 + d2, d0 - d2, -d0 - d1
    offdiag = 2.0 * a01 * a12 * a02
    mixed = (a01 * a01 * t1 + a02 * a02 * t2 + a12 * a12 * t3) / 3.0
    diag = t1 * t2 * t3 / 27.0
    return offdiag + mixed - diag

def delta(A: Symm) -> float:
    """Algorithm 8 discriminant (sum of squares)."""
    a00, a11, a22, p, q, r = A
    d0, d1, d2 = a00 - a11, a00 - a22, a11 - a22
    r1  = p * r * q - q * p * r
    r2  = -p * q * d2 + p * p * r - q * q * r
    r3  = p * r * d1 - p * p * q + q * r * r
    r4  = q * r * d0 + p * r * r - q * q * p
    r5  = p * r * d1 - p * q * p + q * r * r
    r6  = q * r * d0 - p * q * q + p * r * r
    r7  = -q * p * d2 + p * p * r - q * r * q
    r8  = r * d0 * d1 - q * p * d1 + p * p * r - r * r * r
    r9  = r * d0 * d1 - q * p * d0 + q * r * q - r * r * r
    r10 = p * d1 * d2 + q * r * d2 + p * q * q - p * p * p
    r11 = p * d1 * d2 + q * r * d1 + p * r * r - p * p * p
    r12 = -q * d0 * d2 + p * r * d0 + q * r * r - q * q * q
    r13 = q * d0 * d2 + p * r * d2 - p * q * p + q * q * q
    r14 = d0 * d1 * d2 - p * p * d0 + q * q * d1 - r * r * d2
    return (9.0 * r1 * r1 + 6.0 * (r2 * r2 + r3 * r3 + r4 * r4)
            + 8.0 * (r5 * r5 + r6 * r6 + r7 * r7)
            + 2.0 * (r8 * r8 + r9 * r9 + r10 * r10 + r11 * r11
                     + r12 * r12 + r13 * r13)
            + r14 * r14)

def delta_naive(J2: float, J3: float) -> float:
    """Naive discriminant 4J2^3 - 27J3^2, clamped. Benchmarking only."""
    d = 4.0 * J2 * J2 * J2 - 27.0 * J3 * J3
    return 0.0 if d < 0.0 else d

def eigvals(A: Symm, use_naive: bool = False) -> Vec3:
    I1 = i1(A)
    J2 = j2(A)
    if J2 == 0.0:
        # Exact scaled identity (H-Z Eq. 62): eigenvalues are the common
        # diagonal entry, bit-exactly.
        return (A[0], A[0], A[0])
    J3 = j3(A)
    d = delta_naive(J2, J3) if use_naive else delta(A)
    phi = math.atan2(math.sqrt(27.0 * d), 27.0 * J3)
    c = 2.0 * math.sqrt(3.0 * J2)
    lam = lambda k: (I1 + c * math.cos(phi / 3.0 + TWO_PI_OVER_3 * k)) / 3.0
    return (lam(1.0), lam(2.0), lam(3.0))

# ---- Eigenvectors (Eigenvectors.lean) ----

def eigvec_isolated(A: Symm, lam: float) -> Vec3:
    a00, a11, a22, a01, a02, a12 = A
    r0 = (a00 - lam, a01, a02)
    r1 = (a01, a11 - lam, a12)
    r2 = (a02, a12, a22 - lam)
    c0 = cross(r0, r1)   # Eberly's labeling: r0xr1, r0xr2, r1xr2
    c1 = cross(r0, r2)
    c2 = cross(r1, r2)
    d0, d1, d2 = norm_sq(c0), norm_sq(c1), norm_sq(c2)
    if d1 > d0:
        best, dmax = (c2, d2) if d2 > d1 else (c1, d1)
    else:
        best, dmax = (c2, d2) if d2 > d0 else (c0, d0)
    if dmax == 0.0:
        return (1.0, 0.0, 0.0)
    return vscale(best, 1.0 / math.sqrt(dmax))

def orthonormal_complement(w: Vec3) -> tuple:
    if abs(w[1]) < abs(w[0]):
        inv = 1.0 / math.sqrt(w[0] * w[0] + w[2] * w[2])
        u = (-w[2] * inv, 0.0, w[0] * inv)
    else:
        inv = 1.0 / math.sqrt(w[1] * w[1] + w[2] * w[2])
        u = (0.0, w[2] * inv, -w[1] * inv)
    return (u, cross(w, u))

def eigvec_in_plane(A: Symm, v0: Vec3, lam: float) -> Vec3:
    u, v = orthonormal_complement(v0)
    au, av = mul_vec(A, u), mul_vec(A, v)
    m00, m01, m11 = dot(u, au) - lam, dot(u, av), dot(v, av) - lam
    if abs(m00) < abs(m11):
        max_abs = max(abs(m11), abs(m01))
        if max_abs == 0.0:
            return u
        if abs(m11) < abs(m01):
            t = m11 / m01; n = 1.0 / math.sqrt(1.0 + t * t)
            return vsub(vscale(u, t * n), vscale(v, n))
        else:
            t = m01 / m11; n = 1.0 / math.sqrt(1.0 + t * t)
            return vsub(vscale(u, n), vscale(v, t * n))
    else:
        max_abs = max(abs(m00), abs(m01))
        if max_abs == 0.0:
            return u
        if abs(m00) < abs(m01):
            t = m00 / m01; n = 1.0 / math.sqrt(1.0 + t * t)
            return vsub(vscale(u, n), vscale(v, t * n))
        else:
            t = m01 / m00; n = 1.0 / math.sqrt(1.0 + t * t)
            return vsub(vscale(u, t * n), vscale(v, n))

def eigvecs(B: Symm, e: Vec3) -> Mat3:
    l1, l2, l3 = e
    if l2 - l1 < l3 - l2:
        v3 = eigvec_isolated(B, l3)
        v2 = eigvec_in_plane(B, v3, l2)
        return (cross(v2, v3), v2, v3)
    else:
        v1 = eigvec_isolated(B, l1)
        v2 = eigvec_in_plane(B, v1, l2)
        return (v1, v2, cross(v1, v2))

# ---- Pipeline (Eigendecomp.lean) ----

def eigendecomp(A: Symm, use_naive: bool = False) -> tuple:
    mA = max_abs_entry(A)
    if mA == 0.0:
        return ((0.0, 0.0, 0.0),
                ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
    m, _ = math.frexp(mA)        # mA = m * 2^e, m in [0.5, 1)
    s = m / mA                   # exactly 2^(-e): the quotient is a power of 2
    B = tuple(x * s for x in A)
    e = eigvals(B, use_naive)
    Q = eigvecs(B, e)
    inv = mA / m                 # exactly 2^e
    return (tuple(x * inv for x in e), Q)

# ---- Certificates (Certificates.lean) ----

def residual(A: Symm, lam: float, v: Vec3) -> float:
    av = mul_vec(A, v)
    r = (av[0] - lam * v[0], av[1] - lam * v[1], av[2] - lam * v[2])
    return max(abs(r[0]), max(abs(r[1]), abs(r[2])))

def orthogonality_error(Q: Mat3) -> float:
    c1, c2, c3 = Q
    m1 = max(abs(dot(c1, c1) - 1.0), abs(dot(c2, c2) - 1.0))
    m2 = max(abs(dot(c3, c3) - 1.0), abs(dot(c1, c2)))
    m3 = max(abs(dot(c1, c3)), abs(dot(c2, c3)))
    return max(m1, max(m2, m3))

def reconstruct(e: Vec3, Q: Mat3) -> Symm:
    def comp(m: float, v: Vec3) -> Symm:
        return (m * v[0] * v[0], m * v[1] * v[1], m * v[2] * v[2],
                m * v[0] * v[1], m * v[0] * v[2], m * v[1] * v[2])
    b1, b2, b3 = comp(e[0], Q[0]), comp(e[1], Q[1]), comp(e[2], Q[2])
    return tuple(b1[i] + b2[i] + b3[i] for i in range(6))

def reconstruction_error(A: Symm, e: Vec3, Q: Mat3) -> float:
    B = reconstruct(e, Q)
    return max(abs(B[i] - A[i]) for i in range(6))

def certify(A: Symm, e: Vec3, Q: Mat3) -> dict:
    return {
        "maxResidual": max(residual(A, e[i], Q[i]) for i in range(3)),
        "orthogonality": orthogonality_error(Q),
        "reconstruction": reconstruction_error(A, e, Q),
    }

def decomposition(A: Symm, use_naive: bool = False) -> tuple:
    """Uniform harness interface: (eigvals, eigvec columns, certificates)."""
    e, Q = eigendecomp(A, use_naive)
    return (e, Q, certify(A, e, Q))
