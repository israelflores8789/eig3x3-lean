# Contributing to Eig3x3

Thank you for your interest in contributing to **Eig3x3**!

Eig3x3 is a pure Lean 4, zero-dependency, closed-form eigensolver for $3 \times 3$ real symmetric matrices operating over IEEE 754 64-bit floating-point numbers (`Float`). It combines the invariant-based eigenvalue algorithm of Habera & Zilian (2025) with the null-space eigenvector construction of Eberly (2014).

## Getting Started

### Prerequisites

| Tool | Purpose | Recommended Installation |
| :--- | :--- | :--- |
| **`elan`** | Lean version manager | `curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf \| sh` |
| **`uv`** | Python package manager (>= 0.4.0) | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **`just`** | Command runner | `cargo install just` or `brew install just` or package manager |
| **`typos`** | Source code spell checker | `cargo install typos-cli` or `brew install typos-cli` |

The Lean version is pinned in `lean-toolchain` (e.g., `leanprover/lean4:v4.34.0-rc2`); `elan` installs it automatically the first time you run `lake` inside the repository.

### Setup

```bash
git clone https://github.com/israelflores8789/eig3x3-lean.git && cd eig3x3
uv sync --dev               # Python dev dependencies into .venv (parity harness)
uv run pre-commit install   # formatting, YAML, and spell-check hooks
just build                  # core library, CLI harness, and benchmarks
```

Verify the build with a CLI smoke test:

```bash
echo '{"matrices":[[2.0,2.0,2.0,1.0,0.0,1.0]]}' | .lake/build/bin/eig3x3_cli
```

**Expected output:** eigenvalues `[0.58578643762690497, 2.0, 3.4142135623730949]` with certificates (`maxResidual`, `orthogonality`, `reconstruction`) around $10^{-16}$. If the output diverges or throws a JSON parse error, stop and investigate before proceeding.

That's the whole setup. Run `just --list` — or read the `justfile`, which is the source of truth — for every other recipe (individual build targets, parity sub-suites, docs, benchmarks). This guide documents workflows and policy, not individual commands.

### Editor Setup

If you use VS Code, open the workspace root and accept the recommended extensions when prompted (`.vscode/extensions.json`: Lean 4, Python, Ruff, Pyrefly, typos, TOML, and `justfile` syntax). Workspace settings in `.vscode/settings.json` are preconfigured: 2-space Lean indentation with format on save, and Ruff handling Python formatting and import sorting. For other editors, a Lean 4 language server and your usual Python tooling are sufficient.

### Spell Checking

We maintain zero spelling errors across code, comments, and documentation. Run `just spell-diff` to check and `just spell-fix` to apply unambiguous fixes.

> [!IMPORTANT]
> If a flagged word is a legitimate domain term, proper noun (e.g., author names like *Habera*, *Zilian*, *Eberly*), or mathematical identifier, **do not** introduce misspelled workarounds. Propose adding the word to `_typos.toml` under `[default.extend-words]`.

## Core Principles & Architectural Mandates

1. **Zero External Dependencies**: The core library (`src/`) must **never** depend on Mathlib, Batteries, or external C/FFI bindings. It must compile purely against the standard Lean 4 toolchain as a self-contained Lake package.
2. **Deterministic Closed-Form Speed**: Algorithms must be non-iterative, branch-minimized, and execute in sub-microsecond time ($\sim 0.35\,\mu\text{s}$).
3. **Rigorous Numerical Stability**: Floating-point operations must adhere strictly to IEEE 754 double precision (`Float`), using power-of-two preconditioning (`Float.frExp`) and sum-of-squares formulations to eliminate catastrophic cancellation.
4. **Empirical Certification**: Because pure floating-point arithmetic is not formally verified via Lean proofs, runtime certificates (`Eig3x3.certify`) validate accuracy per decomposition. Floating-point proofs are outside the scope of this project.

## Numerical Standards & Invariants

Eig3x3 targets 64-bit double precision floating point (IEEE 754 `binary64`, machine epsilon $\varepsilon = 2^{-52} \approx 2.220446049250313 \times 10^{-16}$).

### Scale-Aware Machine-Epsilon Gates

Both the native Lean test suite (`tests/Tests/Util.lean`) and the Python parity harness (`parity/gates.py`) enforce identical tolerance gates:

| Metric | Tolerance Gate | Formula / Standard |
| :--- | :---: | :--- |
| **Residual Error** ($\|Av_i - \lambda_i v_i\|_\infty$) | $\le 64\varepsilon \cdot \max_{ij} \|A_{ij}\|$ | `certTol(A) = 64 * EPS * maxAbsEntry(A)` |
| **Reconstruction Error** ($\|Q \Lambda Q^T - A\|_\infty$) | $\le 64\varepsilon \cdot \max_{ij} \|A_{ij}\|$ | `certTol(A) = 64 * EPS * maxAbsEntry(A)` |
| **Eigenvalue Parity vs Reference** ($|\lambda_i - \lambda_i^{\text{ref}}|$) | $\le 64\varepsilon \cdot \max_{ij} \|A_{ij}\|$ | `evalTol(A) = 64 * EPS * maxAbsEntry(A)` |
| **Orthogonality Error** ($\|Q^T Q - I\|_\infty$) | $\le 16\varepsilon$ | `orthoTol = 16 * EPS` (dimensionless) |
| **Right-Handed Determinant** ($\det Q$) | $| \det Q - 1.0 | \le 16\varepsilon$ | Right-handed orthonormal frame |
| **Ordering Contract** | $\lambda_0 \le \lambda_1 \le \lambda_2$ | Monotonically non-decreasing |

### Strict Invariant Rules

- **NEVER loosen a tolerance gate or modify expected test outputs to make a failing test pass.** If a gate fails, the mathematical algorithm or implementation has regressed.
- **Preserve Numerical Literals**: Values in test cases (e.g., $1.0 + 1.0\times 10^{-8}$, $2^{\pm 997}$, ULP perturbations) are crafted specifically to exercise corner cases (such as near-double roots or exponent extremes). Do not round or alter them.
- **Never edit generated files directly**: The file `generated/golden.json` is generated via 50-digit `mpmath` arithmetic. Update it only using `just golden`.
- **Mirroring Mandate**: `parity/gates.py` and `tests/Tests/Util.lean` must always remain strictly synchronized.

## Verification & Testing

Run the entire local CI gate before opening a pull request:

```bash
just ci   # Lean tests, Ruff, Pyrefly, pytest parity suite (with JUnit XML), typos
```

The sections below describe what the gate covers. For individual recipes (`just test`, `just parity`, `just lint`, `just typecheck`, `just spell-diff`, …), see `just --list`.

### Lean Test Suite (`just test`)

The test runner (`tests/Main.lean`) executes five dedicated modules:

1. **`KnownAnswer`**: Exact arithmetic checks on textbook matrices, identity, zero matrix, and 3D algebra operators.
2. **`Golden`**: Bit-exact comparison against 50-digit `mpmath` reference vectors loaded via exact dyadic pairs `[sig, exp]` (bypassing float parser rounding).
3. **`Properties`**: 5,000 randomized matrices generated with SplitMix64 verifying trace invariants ($\text{tr}(A) = \sum \lambda_i$), determinant preservation ($\det A = \prod \lambda_i$), scale invariance, and right-handedness ($\det Q = 1$).
4. **`Regression`**: Pinned historical boundary cases, including the Habera–Zilian $r_{10}$ discriminant check and clustered perturbation paths.
5. **`Certificates`**: Verification that all matrices in the case zoo satisfy the $64\varepsilon / 16\varepsilon$ runtime gates.

The Lean suite certifies its own output; it does not depend on Python to be correct.

### Python Parity Suite (`just parity`)

The Python harness is a differential reference, not a second implementation: it feeds matrices through the compiled CLI and compares against NumPy/LAPACK (`dgeev`/`dsyev`) and 50-digit `mpmath` golden values. The default run uses n=1000 random samples; `just parity-parallel` distributes it across cores with pytest-xdist. Named sub-suites (`zoo`, `random`, `logscale`, `paths`, `smallscale`, `frontier`, `golden`) target specific regimes — e.g., `just parity frontier` for the $2^{\pm 997}$ dynamic-range cases. See the `justfile` for the full list.

### Errata (`just errata`)

Executes exhibits discovered during the initial algorithm implementation (e.g., the Habera–Zilian §7 Algorithm 8 $r_{10}$ paper typo).

> [!NOTE]
> `parity/errata.py` is a permanent record of mathematical errata. Do not modify existing exhibits.

### Benchmarks (`just bench`)

Compares in-process Lean 4 execution against single-matrix and vectorized `numpy.linalg.eigh`. Timings cover the decomposition only; the `certify` step is excluded because NumPy/LAPACK has no equivalent cost. Residual validation belongs to the test suites, not the benchmark.

## CI

The CI workflow (`.github/workflows/ci.yml`) runs on pull requests targeting `main` or `dev` and on pushes to `dev`:

- Matrix: `ubuntu-latest` across `leanprover/lean4:v4.27.0` (minimum supported), `v4.30.0`, and `v4.34.0-rc2` (pinned development toolchain).
- Steps: `typos` spell check, `lake build` of all targets, then `just ci`. The parity JUnit report is uploaded as a workflow artifact (`parity-junit-report-*`).

`just ci` locally reproduces the same gate the CI enforces.

## Code Style & Conventions

### Lean 4

- **Language Settings**: We enforce `autoImplicit = false`, `builtinLint = true`, `linter.unusedVariables = true`, `linter.deprecated = true`, and `linter.missingDocs = true` in `lakefile.toml`.
- **Naming Conventions**: Types, structures, and typeclasses use `UpperCamelCase` (e.g., `SymmMat3`, `Vec3`, `Decomposition`); functions, definitions, and theorems use `lowerCamelCase` (e.g., `eigendecomp`, `maxAbsEntry`, `residual`). Place Unicode math operators under `open scoped Eig3x3` to avoid polluting downstream namespaces.
- **Explicit Types**: Always provide explicit return types and argument types for public declarations.
- **Float Operations**: Prefer idiomatic Lean float operations (e.g., `x.abs`, `x.sqrt`, `x.frExp`) and use `PowNat` / `^ⁿ` for integer exponentiation.

### Documentation & Docstrings

Because `linter.missingDocs = true` is enabled, all public types, structures, fields, and definitions in `src/` **must** have docstrings (`/-- ... -/`). Module documentation (`/-! ... -/`) should describe high-level mathematical concepts and design decisions. Cite sources with bracket syntax linked to `references.bib` (e.g., `[HaberaZilian2025]`, `[Eberly2014]`); when introducing new mathematical sources, add entries to both `references.bib` and `CITATION.cff`.

### Python

- **Formatting & Linting**: Enforced via Ruff (`pyproject.toml`, line length 88).
- **Type Annotations**: All functions in `parity/` must have complete type signatures validated by Pyrefly (`just typecheck`).
- **No Side Effects**: Differential test scripts must avoid global state mutations and support `pytest-xdist` parallel workers.

## Documentation (`doc-gen4`)

API documentation is generated with `doc-gen4` from the isolated subproject in `docbuild/`:

- **Zero-Dependency Guarantee**: The root `lakefile.toml` must remain dependency-free. `doc-gen4` is required only inside `docbuild/lakefile.toml`.
- **Toolchain Pinning**: The `rev` in `docbuild/lakefile.toml` must match the root `lean-toolchain`. When upgrading the Lean toolchain, bump both files together.

Use `just docs-dev` to build locally and `just docs-serve` to preview on http://localhost:8000.

> [!CAUTION]
> Do not use `just docs` for local development; it is reserved for the automated GitHub Pages deployment workflow (`.github/workflows/docs.yml`).

## Areas for Contribution (Roadmap)

We welcome contributions in several key areas. If you plan to work on a significant feature, please open an issue first for architectural discussion.

1. **Iterative Fallback Solver**: An optional, robust iterative eigensolver (e.g., Jacobi rotations or QR with shifts) for pathological inputs where closed-form transcendental evaluations lose accuracy.
2. **Complex Hermitian $3 \times 3$ Eigensolver**: Real eigenvalues and unitary eigenvector matrices ($Q^* Q = I$, $\det Q = 1$) for $A = A^*$.
3. **Additional $3 \times 3$ Decompositions**: SVD via eigendecomposition of $A^T A$; generalized symmetric eigensolvers ($A v = \lambda B v$); polar decomposition ($A = U P$); fast Cholesky and $LDL^T$ factorizations.
4. **Formal Verification & Bounds Proofs**: A separate, opt-in companion package (e.g., `Eig3x3.Verify` with Mathlib) formalizing correctness proofs, condition numbers, and error bounds — the core library stays proof-free.
5. **Platform Portability & Hardware Acceleration**: Characterization and optimization across ARM64 (NEON), x86-64 (AVX2/FMA), WebAssembly, and embedded targets; SIMD exploration once native vector primitives stabilize in Lean 4.
6. **Differential Parity Corpus Expansion**: Additional adversarial generators (e.g., Wilkinson-type clustered spectra, random orthogonal similarity orbits $Q D Q^T$).

## Submitting Changes

### 1. Issue First

For non-trivial enhancements, bug fixes, or new mathematical routines, open an issue first using one of the GitHub templates: **Bug Report**, **Numerical Discrepancy** (precision/tolerance issues), or **Feature Request**.

### 2. Fork and Branch

Fork the repository, clone your fork, and create a descriptive branch off `dev` (or `main` for critical patches):

```bash
git checkout -b feat/iterative-fallback
# or: git checkout -b fix/atan2-quadrant-guard
```

### 3. Validate Locally

```bash
just spell-diff && just ci
```

### 4. Commit Guidelines

We follow [Conventional Commits](https://www.conventionalcommits.org/) — `feat:`, `fix:`, `perf:`, `test:`, `docs:`, `refactor:`. Keep messages concise, descriptive, and focused on the *why* of the change (e.g., `fix: correct angle normalization in planar reduction`).

### 5. Submitting a Pull Request

1. Push your branch and open a Pull Request targeting the `dev` branch.
2. Fill out the PR template completely: summary and motivation, algorithmic/mathematical references, verification steps taken, and confirmation that all numerical gates and CI checks pass.
3. Engage with code review feedback promptly. Once approved and all status checks pass, your contribution will be merged!

---

*Thank you for helping make scientific computing in Lean 4 faster, more robust, and more accessible!*
