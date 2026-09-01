# Eig3x3

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Lean 4](https://img.shields.io/badge/Lean-v4.34.0--rc2-blue.svg)](https://github.com/leanprover/lean4)
[![Reservoir](https://img.shields.io/badge/Reservoir-v1.0.0-green.svg)](https://reservoir.lean-lang.org/@israelflores8789/Eig3x3)
[![Zero Dependencies](https://img.shields.io/badge/Dependencies-0-brightgreen.svg)](#)

A pure, high-performance **Lean 4** library for the closed-form eigendecomposition of $3 \times 3$ real symmetric matrices over IEEE 754 64-bit floating point (`Float`).

**Eig3x3** delivers machine-precision eigenvalues and eigenvectors without external C/FFI bindings and without Mathlib dependencies.

## Table of Contents

- [Overview & Highlights](#overview--highlights)
- [Motivation](#motivation)
- [Installation](#installation)
  - [Supported Toolchains & Downstream Compatibility](#supported-toolchains--downstream-compatibility)
- [Quick Start](#quick-start)
- [Core Architecture & Numerical Methods](#core-architecture--numerical-methods)
- [Accuracy & Numerical Stability](#accuracy--numerical-stability)
- [3D Vector & Matrix Algebra Vocabulary](#3d-vector--matrix-algebra-vocabulary)
- [Runtime Certificates](#runtime-certificates)
- [Verification & Testing](#verification--testing)
  - [Lean Test Suite](#lean-test-suite)
  - [Python Parity & Differential Suite](#python-parity--differential-suite)
- [Performance Benchmarks](#performance-benchmarks)
- [Citation](#citation)
- [Credits](#credits)

## Overview & Highlights

- **Pure Lean 4 / Zero Dependencies**: No C/FFI toolchain dependencies, no Mathlib dependency. Fully self-contained and ready as a drop-in Lake package.
- **Numerically Stabilized Closed-Form**: Combines the invariant-based eigenvalue algorithm of **Habera & Zilian (2025)** with the robust null-space eigenvector construction of **Eberly (2014)**.
- **Deterministic Sub-Microsecond Speed**: Closed-form computation eliminates loop branches and iterative convergence checks, achieving $\sim 0.35\,\mu\text{s}$ per decomposition ($\sim 2.8\times 10^6$ matrices/sec).
- **Bit-Transparent Preconditioning**: Employs `Float.frExp` power-of-two scaling to guard against underflow/overflow over 600 orders of magnitude ($10^{-300}$ to $10^{300}$) with zero mantissa rounding distortion.
- **Runtime Error Certificates**: Built-in verification module (`certify`) computing per-instance residual, orthonormality, and reconstruction errors.
- **Mathlib-Consistent Operator Layer**: Opt-in Unicode notation (`open scoped Eig3x3`) providing `Qᵀ`, `M⁻¹`, `u ⬝ᵥ v`, `A ⬝ₘ B`, `u ⊗ᵥ v`, `u ⨯₃ v`, `s • v`, `u ⊙ v`, `‖v‖`, `|x|`, and fast natural powers `x ^ⁿ n`.

## Motivation

Eigendecomposition of $3 \times 3$ symmetric matrices is a core primitive across scientific computing, geometric processing, robotics, and physics simulations.

Standard production workflows typically bind via FFI to iterative LAPACK routines (e.g., `dsyev`). While robust, iterative routines:
1. Impose external C build dependencies and FFI overhead.
2. Suffer variable execution times dependent on convergence criteria.
3. Complicate deployment across pure-Lean environments, embedded contexts, and web/Wasm targets.

Classical closed-form solutions (Cardano/Viète) avoid iterations but suffer catastrophic floating-point cancellation near degenerate (double or triple) eigenvalues. **Eig3x3** brings recent advances in numerically stable closed-form solvers to the Lean 4 ecosystem, providing a lightweight, robust foundation for scientific computing in Lean.

## Installation

### Using Lean Reservoir / `lakefile.toml` (Recommended)

Add `Eig3x3` to your project's `lakefile.toml`:

```toml
[[require]]
name = "Eig3x3"
scope = "israelflores8789"
version = ">= 1.0.0"
```

*Alternatively, require directly from Git:*

```toml
[[require]]
name = "Eig3x3"
git = "https://github.com/israelflores8789/eig3x3-lean"
rev = "v1.0.0"
```

### Using `lakefile.lean`

```lean
require «Eig3x3» from "israelflores8789" / "eig3x3-lean"
-- or from git:
require «Eig3x3» from git
  "https://github.com/israelflores8789/eig3x3-lean" @ "v1.0.0"
```

Then run `lake update` and `lake build`.

### Supported Toolchains & Downstream Compatibility

> [!NOTE]
> When `Eig3x3` is consumed as a Lake dependency in your project, its internal `lean-toolchain` file is **ignored** by Lake — your root project's toolchain governs the entire workspace.

Because `Eig3x3` has zero external dependencies (no Mathlib, no Batteries) and uses only standard Lean 4 core language primitives (`Float`, structures, and typeclasses), it maintains broad compatibility across Lean 4 versions:

- **Supported Toolchains**: Lean 4 `v4.27.0` or later.
- **CI Toolchain Matrix**: Built and continuously tested across active stable and candidate Lean 4 releases, including:
  - `v4.27.0`
  - `v4.30.0`
  - `v4.34.0-rc2`

## Quick Start

```lean
import Eig3x3

open scoped Eig3x3

def main : IO Unit := do
  -- 1. Define a symmetric 3×3 matrix via its 6 unique entries:
  --    [[2.0, 1.0, 0.0],
  --     [1.0, 2.0, 1.0],
  --     [0.0, 1.0, 2.0]]
  let A : Eig3x3.SymmMat3 := ⟨2.0, 2.0, 2.0, 1.0, 0.0, 1.0⟩

  -- 2. Compute the full eigendecomposition: A = Q Λ Qᵀ
  let decomp := Eig3x3.eigendecomp A
  let e := decomp.eigvals  -- Ordered: l₀ ≤ l₁ ≤ l₂
  let Q := decomp.eigvecs  -- Right-handed orthonormal matrix (det Q = 1)

  IO.println s!"Eigenvalues: [{e.l₀}, {e.l₁}, {e.l₂}]"
  -- Expected: [0.5857864376269049, 2.0, 3.414213562373095]

  -- 3. Verify numerical quality with runtime certificates:
  let certs := Eig3x3.certify A decomp
  IO.println s!"Max Residual ‖Av - λv‖∞: {certs.maxResidual}"
  IO.println s!"Orthogonality Error:    {certs.orthogonality}"
  IO.println s!"Reconstruction Error:   {certs.reconstruction}"
  -- Errors are typically ≈ 1e-16 to 1e-15

  -- 4. Convenient vector & matrix arithmetic:
  let v : Eig3x3.Vec3 := ⟨1.0, 0.0, 0.0⟩
  let transformed := Q * v
  let norm := ‖transformed‖
  IO.println s!"Transformed vector norm: {norm}"
```

## Core Architecture & Numerical Methods

The calculation of eigenvectors is fundamentally a null-space computation of $(A - \lambda I)$ that directly consumes the computed eigenvalues $\lambda$. Consequently, eigenvalue accuracy dictates the stability of the entire pipeline.

```
                  ┌─────────────────────────────────────────┐
                  │          Input SymmMat3 (A)             │
                  └────────────────────┬────────────────────┘
                                       │
                         [Float.frExp Preconditioning]
                                       │
                  ┌────────────────────▼────────────────────┐
                  │    1. Habera–Zilian (2025) Pipeline     │
                  │   - Trace recentering (I₁)              │
                  │   - Deviatoric invariants (J₂, J₃)      │
                  │   - Sum-of-squares discriminant (Δ)     │
                  │   - Quadrant-safe atan2 angle (φ)       │
                  │   - 3-element permutation sort guard    │
                  └────────────────────┬────────────────────┘
                                       │ Ordered Eigenvalues (l₀ ≤ l₁ ≤ l₂)
                  ┌────────────────────▼────────────────────┐
                  │    2. Eberly (2014) Eigenvector Engine  │
                  │   - Spectral gap comparison             │
                  │   - Cross-product null-space (isolated) │
                  │   - 2×2 planar orthogonal complement    │
                  │   - Right-handed completion (c₁ ⨯₃ c₂)  │
                  └────────────────────┬────────────────────┘
                                       │
                            [Power-of-2 Rescaling]
                                       │
                  ┌────────────────────▼────────────────────┐
                  │    Decomposition { eigvals, eigvecs }   │
                  │      + Optional Runtime Certify         │
                  └─────────────────────────────────────────┘
```

### Key Algorithmic Improvements & Fixes

1. **Sum-of-Squares Discriminant ($\Delta$)**: Uses Habera–Zilian Algorithm 8 sum-of-squares formulation ($\Delta = 4J_2^3 - 27J_3^2$) rather than subtractive cubics, preventing catastrophic cancellation when eigenvalues are close. Includes the verified $r_{10}$ correction (`+ q·r·d₂`).
2. **Quadrant-Safe $\text{atan2}$ Angle Formulation**: Replaces the classical $\arccos$ formulation with $\varphi = \text{atan2}(\sqrt{27\Delta}, 27J_3)$, maintaining forward stability across the full domain.
3. **Ordering Contract Guard**: While Habera–Zilian guarantees $\lambda_0 \le \lambda_1 \le \lambda_2$ in exact arithmetic, floating-point evaluation of transcendental $\cos$ at degenerate angles can vary by $\sim 1\text{ ULP}$. A final 3-element compare-exchange sort enforces the strict ordering contract.
4. **Spectral Gap Isolation**: Compares $(\lambda_1 - \lambda_0)$ vs $(\lambda_2 - \lambda_1)$ to construct the isolated eigenvector first from the longest row cross-product, completing the remaining eigenvectors via 2D planar reduction and cross product.

## Accuracy & Numerical Stability

Max residual error ($\|Av - \lambda v\|_\infty / \|A\|_\infty$) across characteristic problem regimes in double precision (`Float` / `float64`):

| Test Regime | Naive Cardano / Viète | NumPy / LAPACK (`dsyev`) | `Eig3x3` (Lean 4) | Status / Notes |
| :--- | :---: | :---: | :---: | :--- |
| **Uniform Random** ($A_{ij} \in [-1, 1]$) | $\sim 10^{-15}$ | $\sim 10^{-16}$ | $\mathbf{\sim 10^{-16}}$ | Machine precision across all solvers |
| **Near-Double Eigenvalues** ($\delta = 10^{-8}$) | $\approx 3.6 \times 10^{-9}$ | $\approx 2.2 \times 10^{-16}$ | $\mathbf{\approx 2.2 \times 10^{-16}}$ | Naive loses ~8 digits; Eig3x3 matches LAPACK |
| **Near-Triple Eigenvalues** ($\delta = 10^{-8}$) | $\approx 1.5 \times 10^{-8}$ | $\approx 2.2 \times 10^{-16}$ | $\mathbf{\approx 2.2 \times 10^{-16}}$ | $J_2 \to 0$ stabilized via diagonal differences |
| **Scaled Identity** ($cI$) | $\sim 10^{-15}$ | $0.0$ | $\mathbf{0.0}$ | Exact zero fast path ($J_2 = 0$) |
| **Dynamic Range** ($10^{-300}$ to $10^{300}$) | Overflow / Underflow | Fails / Subnormal | $\mathbf{\sim 10^{-16}}$ | Bit-transparent `Float.frExp` scaling |

## 3D Vector & Matrix Algebra Vocabulary

`import Eig3x3` includes a complete, standalone 3D linear algebra toolkit with `Float`-type structures `Vec3` and `Mat3` for a 3-element vector and 3-column matrix, respectively, without Mathlib dependencies. Activating `open scoped Eig3x3` enables the following mathematical notations with Mathlib conventions:

| Notation | Operation | Lean Declaration | Precedence / Associativity |
| :--- | :--- | :--- | :--- |
| `u + v` / `A + B` | Addition | `HAdd.hAdd` | `infixl:65` |
| `u - v` / `A - B` | Subtraction | `HSub.hSub` | `infixl:65` |
| `-v` / `-A` | Negation | `Neg.neg` | Prefix |
| `A * B` | Matrix product | `HMul.hMul` | `infixl:70` |
| `A * v` | Matrix-vector product | `HMul.hMul` | `infixl:70` |
| `v * A` | Row-vector matrix product ($v^T A$) | `HMul.hMul` | `infixl:70` |
| `v / s` / `A / s` | Scalar division | `HDiv.hDiv` | `infixl:70` |
| `u ⊗ᵥ v` | Vector outer product ($u v^T$) | `Vec3.outer` | `infixl:70` |
| `u ⬝ᵥ v` | Vector dot product | `Vec3.dot` | `infixl:72` |
| `A ⬝ₘ B` | Matrix Frobenius inner product | `Mat3.dot` | `infixl:72` |
| `s • v` / `s • A` | Scalar multiplication | `HSMul.hSMul` | `infixr:73` |
| `u ⨯₃ v` | Vector cross product ($\mathbb{R}^3$) | `Vec3.cross` | `infixl:74` |
| `x ^ⁿ n` / `A ^ⁿ n`| Fast natural power (repeated mul / squaring) | `PowNat.powNat`| `infixr:80` |
| `u ⊙ v` / `A ⊙ B` | Hadamard (entrywise) product | `Hadamard.hadamard`| `infixl:100` |
| `\|x\|` / `\|v\|` / `\|A\|` | Absolute value / entrywise magnitude | `Abs.abs` | Delimited (`\|v\|`) |
| `‖v‖` / `‖A‖` | Euclidean (vector) / Frobenius (matrix) norm | `Norm.norm` | Delimited (`‖v‖`) |
| `‖v‖²` / `‖A‖²` | Squared norm | `NormSq.normSq` | Delimited (`‖v‖²`) |
| `Aᵀ` | Matrix transpose | `Mat3.transpose` | `postfix:max` |
| `A⁻¹` | Matrix inverse (via cofactors) | `Inv.inv` | `postfix:max` |

## Runtime Certificates

Because mathematical proofs over IEEE 754 `Float` are inherently empirical without interval arithmetic, `Eig3x3.certify` provides runtime numerical validation:

```lean
let certs := Eig3x3.certify A decomp
```

- **`maxResidual`**: $\max_i \|A v_i - \lambda_i v_i\|_\infty$ (checks eigenpair validity).
- **`orthogonality`**: $\max_{i,j} |c_i \cdot c_j - \delta_{ij}|$ (checks orthonormality of $Q$).
- **`reconstruction`**: $\max_{i,j} |(Q \Lambda Q^T)_{ij} - A_{ij}|$ (checks spectral synthesis).

### Scale-Aware Quality Standard

All decompositions satisfy strict machine-epsilon gates ($\varepsilon = 2^{-52} \approx 2.22 \times 10^{-16}$):
- $\text{Residual} \le 64\varepsilon \cdot \max_{ij} |A_{ij}|$
- $\text{Reconstruction} \le 64\varepsilon \cdot \max_{ij} |A_{ij}|$
- $\text{Orthogonality} \le 16\varepsilon$ (dimensionless)

## Verification & Testing

Cloning the repository offers various tests that validate the eigensolver's computation. This repo uses [`just`](https://github.com/casey/just) as the command runner.

### Lean Test Suite

Run the full native Lean verification suite with `just test` or `lake test`:

```bash
just test
# or
lake test
```

The Lean test harness executes:
- **`KnownAnswer`**: Exact-arithmetic validation for vectors, matrices, operators, and special matrices (zero matrix, scaled identity, textbook cases).
- **`Golden`**: Bit-exact verification against 50-digit `mpmath` reference vectors (loaded via dyadic pairs $[sig, exp]$ to prevent decimal parser rounding).
- **`Properties`**: 5,000 deterministic pseudo-random matrices (via splitmix64 PRNG) testing trace invariants, determinant identities, right-handedness ($\det Q = 1$), and scale invariance.
- **`Regression`**: Pinned historical boundary cases including the Habera–Zilian $r_{10}$ discriminant check and clustered perturbation paths.
- **`Certificates`**: Runtime error assertions on the curated case zoo under the $64\varepsilon / 16\varepsilon$ standard.

### Python Parity & Differential Suite

A comprehensive differential test harness validates the Lean CLI binary against NumPy/LAPACK and 50-digit `mpmath` references.

```bash
# Run the full parity suite (default n=1000, seed=42)
just parity

# Run in parallel across all CPU cores (pytest-xdist) (useful for large n)
just parity-parallel

# Run granular sub-suites
just parity zoo          # 7 curated edge cases
just parity random       # Uniform random symmetric matrices A_ij ∈ [-1, 1]
just parity logscale     # Matrices spanning 600 orders of magnitude (10^-300 to 10^300)
just parity paths        # Adversarial perturbation paths (diag(-1,1,1+δ), diag(1,1,1+δ))
just parity smallscale   # Vanishing scale double eigenvalues (s, s, 2s)
just parity frontier     # Extreme dynamic range cases (2^±997)
just parity golden       # 50-digit mpmath golden cases

# Generate CI JUnit XML report
just parity-ci
```

## Performance Benchmarks

Benchmarks measured on $n = 20{,}000$ random symmetric matrices comparing in-process Lean 4 native execution against optimized NumPy/LAPACK (`dgeev`/`dsyev`):

| Implementation / Workload | Latency (µs / matrix) | Throughput (matrices / sec) | Speedup vs. Single LAPACK |
| :--- | :---: | :---: | :---: |
| **`Eig3x3` (Lean in-process)** | **0.35 µs** | **~2,850,000 / s** | **~35× faster** |
| **`Eig3x3` + `certify` (Lean in-process)** | **0.44 µs** | **~2,270,000 / s** | **~28× faster** |
| `numpy.linalg.eigh` (single-matrix loop) | 12.90 µs | ~77,500 / s | 1.0× (baseline) |
| `numpy.linalg.eigh` (vectorized batch) | 1.55 µs | ~645,000 / s | ~8.3× faster |
| `Eig3x3` CLI (end-to-end JSON IPC) | 6.22 µs | ~160,000 / s | ~2.1× faster |

## Citation

If you use `Eig3x3` in your academic research or software project, please cite:

```bibtex
@software{FloresArbolay_Eig3x3_2026,
  author       = {Israel D. Flores-Arbolay},
  title        = {Eig3x3: Pure Lean 4 Closed-Form 3x3 Symmetric Eigensolver},
  year         = {2026},
  publisher    = {GitHub},
  url          = {https://github.com/israelflores8789/eig3x3-lean},
  license      = {Apache-2.0}
}
```

## Credits

**This library is an implementation of the algorithms defined in the following scholarly articles and publishings:**

David H. Eberly, "A Robust Eigensolver for 3x3 Symmetric Matrices," Geometric Tools,
LLC, 2014. <https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf>

Michal Habera and Andreas Zilian, "Numerically stable evaluation of closed-form
expressions for eigenvalues of 3×3 matrices," arXiv:2511.00292 [math.NA], 2025.
<https://doi.org/10.48550/arXiv.2511.00292>

Michal Habera and Andreas Zilian, "Symbolic spectral decomposition of 3x3 matrices,"
arXiv:2111.02117 [math.NA], 2021. <https://doi.org/10.48550/arXiv.2111.02117>

**While not referenced directly, each author maintains C implementations of their work which can be found on their GitHub here:**

David Eberly's C implementation of their 2014 paper at Geometric Tools: [`davideberly/GeometricTools`](https://github.com/davideberly/GeometricTools)

Michal Habera and Andreas Zilian's C implementation of their 2025 paper: [`michalhabera/eig3x3`](https://github.com/michalhabera/eig3x3)

**Please give them your support!**
