# Contributing to Eig3x3

Thank you for your interest in contributing to **Eig3x3**!

Eig3x3 is a pure Lean 4, zero-dependency, closed-form eigensolver for $3 \times 3$ real symmetric matrices operating over IEEE 754 64-bit floating-point numbers (`Float`). It combines the invariant-based eigenvalue algorithm of Habera & Zilian (2025) with the null-space eigenvector construction of Eberly (2014).

This document outlines development setup, coding standards, numerical precision requirements, testing workflows, and the contribution lifecycle.

## Table of Contents

- [Core Principles & Architectural Mandates](#core-principles--architectural-mandates)
- [Development Environment Setup](#development-environment-setup)
  - [Prerequisites](#prerequisites)
  - [VS Code Setup & Recommended Extensions](#vs-code-setup--recommended-extensions)
  - [Python Virtual Environment with `uv`](#python-virtual-environment-with-uv)
  - [Spell Checking with `typos`](#spell-checking-with-typos)
  - [Pre-commit Hooks](#pre-commit-hooks)
- [Building the Project](#building-the-project)
  - [Build Targets](#build-targets)
  - [CLI Smoke Test](#cli-smoke-test)
- [IEEE 754 Numerical Standards & Invariants](#ieee-754-numerical-standards--invariants)
  - [Scale-Aware Machine-Epsilon Gates](#scale-aware-machine-epsilon-gates)
  - [Strict Invariant Rules](#strict-invariant-rules)
  - [Precision & Literal Policy](#precision--literal-policy)
- [Verification & Testing](#verification--testing)
  - [Lean Test Suite (`lake test`)](#lean-test-suite-lake-test)
  - [Python Parity & Differential Suite (`pytest`)](#python-parity--differential-suite-pytest)
  - [Errata & Historical Regressions](#errata--historical-regressions)
  - [Benchmarks](#benchmarks)
- [CI/CD Lifecycle & Test Reports](#cicd-lifecycle--test-reports)
- [Code Style & Conventions](#code-style--conventions)
  - [Lean 4 Conventions](#lean-4-conventions)
  - [Documentation & Docstrings](#documentation--docstrings)
  - [Python Conventions](#python-conventions)
- [Documentation System (`doc-gen4`)](#documentation-system-doc-gen4)
- [Areas for Contribution (Roadmap)](#areas-for-contribution-roadmap)
- [Contribution Procedure & Git Workflow](#contribution-procedure--git-workflow)
  - [1. Issue First](#1-issue-first)
  - [2. Fork and Branch](#2-fork-and-branch)
  - [3. Local Validation](#3-local-validation)
  - [4. Commit Guidelines](#4-commit-guidelines)
  - [5. Submitting a Pull Request](#5-submitting-a-pull-request)

---

## Core Principles & Architectural Mandates

1. **Zero External Dependencies**: The core library (`src/`) must **never** depend on Mathlib, Batteries, or external C/FFI bindings. It must compile purely against the standard Lean 4 toolchain as a self-contained Lake package.
2. **Deterministic Closed-Form Speed**: Algorithms must be non-iterative, branch-minimized, and execute in sub-microsecond time ($\sim 0.35\,\mu\text{s}$).
3. **Rigorous Numerical Stability**: Floating-point operations must adhere strictly to IEEE 754 double precision (`Float`), using power-of-two preconditioning (`Float.frExp`) and sum-of-squares formulations to eliminate catastrophic cancellation.
4. **Empirical Certification**: Because pure floating-point arithmetic is not formally verified via Lean proofs, runtime certificates (`Eig3x3.certify`) validate accuracy per decomposition.

---

## Development Environment Setup

### Prerequisites

Install the following tools on your system:

| Tool | Purpose | Recommended Installation |
| :--- | :--- | :--- |
| **`elan`** | Lean version manager | `curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf \| sh` |
| **`uv`** | Fast Python package manager (>= 0.4.0) | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **`just`** | Command runner | `cargo install just` or `brew install just` or package manager |
| **`typos`** | Fast source code spell checker | `cargo install typos-cli` or `brew install typos-cli` |

The repository pins the active Lean version in `lean-toolchain` (e.g., `leanprover/lean4:v4.34.0-rc2`). Running `lake` or `elan` commands inside the repository will automatically resolve the pinned toolchain.

### VS Code Setup & Recommended Extensions

If you use Visual Studio Code, open the workspace root. When prompted, install the recommended extensions defined in `.vscode/extensions.json`:

- **Lean 4** (`leanprover.lean4`): Lean language support and integrated LSP server.
- **Python** (`ms-python.python`): Python language support.
- **Debugpy** (`ms-python.debugpy`): Python debugging.
- **Pyrefly** (`meta.pyrefly`): High-performance type checker for the parity suite.
- **Ruff** (`charliermarsh.ruff`): Fast Python linter and formatter.
- **Just Syntax** (`nefrob.vscode-just-syntax`): Syntax highlighting for `justfile`.
- **Even Better TOML** (`tamasfe.even-better-toml`): TOML support with schema validation for `lakefile.toml` and `pyproject.toml`.
- **Typos** (`tekumara.typos-vscode`): In-editor spell checking.

#### LSP & Formatting Settings

Workspace settings in `.vscode/settings.json` are preconfigured:
- Lean: indentation is 2 spaces, format on save enabled (`leanprover.lean4`).
- Python: formatting and import sorting handled by Ruff on save; `python.languageServer` set to `None` in favor of Pyrefly.

### Python Virtual Environment with `uv`

The Python differential harness (`parity/`) validates Lean output against NumPy, LAPACK, and 50-digit `mpmath` reference values.

Initialize and synchronize the virtual environment:

```bash
# Install development dependencies into local .venv
uv sync --dev
```

You can run commands directly through `uv run` or activate the virtual environment manually:

```bash
# Option A: Activate in your shell
source .venv/bin/activate

# Option B: Run via uv or just (recommended)
uv run pytest parity/parity.py
# or
just parity
```

### Spell Checking with `typos`

We maintain zero spelling errors across all code, comments, and documentation.

```bash
# Check for spelling errors without modifying files
just spell-diff

# Automatically apply fixes for unambiguous typos
just spell-fix
```

> [!IMPORTANT]
> If a flagged word is a legitimate domain term, proper noun (e.g., author names like *Habera*, *Zilian*, *Eberly*), or mathematical identifier, **do not** introduce misspelled workarounds. Propose adding the word to `_typos.toml` under `[default.extend-words]`.

### Pre-commit Hooks

Install pre-commit hooks to enforce formatting, YAML validity, and spell checks on every commit:

```bash
uv run pre-commit install
```

---

## Building the Project

We use `just` recipes to standardize builds across Lean and Python. Run `just --list` to inspect all available tasks.

### Build Targets

```bash
# Build the core Lean library (Eig3x3)
just build
# or: lake build

# Build the JSON CLI harness for Python parity testing
just build-cli
# or: lake build eig3x3_cli

# Build the in-process microbenchmark executable
just build-bench
# or: lake build eig3x3_bench

# Build all Lean binaries and libraries
just build-all
```

### CLI Smoke Test

After compiling the CLI binary, verify communication with this quick smoke test:

```bash
echo '{"matrices":[[2.0,2.0,2.0,1.0,0.0,1.0]]}' | .lake/build/bin/eig3x3_cli
```

**Expected output:**
- Eigenvalues: `[0.58578643762690497, 2.0, 3.4142135623730949]`
- Certificates (`maxResidual`, `orthogonality`, `reconstruction`): $\approx 10^{-16}$.

If the output diverges or throws a JSON parse error, stop and investigate before proceeding.

---

## IEEE 754 Numerical Standards & Invariants

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

---

## Verification & Testing

Every proposed change must pass the full verification matrix.

```
Local Verification Pipeline:
  just test      ──> Lean unit & property tests (lake test)
  just lint      ──> Ruff linter & formatter checks
  just typecheck ──> Pyrefly static type checker
  just parity-ci ──> Pytest parity suite with JUnit XML report
  just spell-diff──> Typos check
```

### Lean Test Suite (`lake test`)

Run the complete native Lean verification suite:

```bash
just test
# or: lake test
```

The test runner (`tests/Main.lean`) executes five dedicated test modules:
1. **`KnownAnswer`**: Exact arithmetic checks on textbook matrices, identity, zero matrix, and 3D algebra operators.
2. **`Golden`**: Bit-exact comparison against 50-digit `mpmath` reference vectors loaded via exact dyadic pairs `[sig, exp]` (bypassing float parser rounding).
3. **`Properties`**: 5,000 randomized matrices generated with SplitMix64 verifying trace invariants ($\text{tr}(A) = \sum \lambda_i$), determinant preservation ($\det A = \prod \lambda_i$), scale invariance, and right-handedness ($\det Q = 1$).
4. **`Regression`**: Pinned historical boundary cases, including the Habera–Zilian $r_{10}$ discriminant check and clustered perturbation paths.
5. **`Certificates`**: Verification that all matrices in the case zoo satisfy the $64\varepsilon / 16\varepsilon$ runtime gates.

### Python Parity & Differential Suite (`pytest`)

The Python test suite compares the compiled Lean CLI binary against NumPy/LAPACK (`dgeev`/`dsyev`):

```bash
# Run full parity suite with default n=1000 random samples
just parity

# Run in parallel across all CPU cores (pytest-xdist)
just parity-parallel

# Run granular sub-suites
just parity zoo          # 7 curated edge cases
just parity random       # Uniform random symmetric matrices A_ij ∈ [-1, 1]
just parity logscale     # Matrices spanning 600 orders of magnitude (10^-300 to 10^300)
just parity paths        # Adversarial perturbation paths (diag(-1,1,1+δ), diag(1,1,1+δ))
just parity smallscale   # Vanishing scale double eigenvalues (s, s, 2s)
just parity frontier     # Extreme dynamic range cases (2^±997)
just parity golden       # 50-digit mpmath golden cases
```

### Errata & Historical Regressions

Run `just errata` to execute exhibits discovered during the initial algorithm implementation (e.g., the Habera–Zilian §7 Algorithm 8 $r_{10}$ paper typo).

> [!NOTE]
> `parity/errata.py` is a permanent record of mathematical errata. Do not modify existing exhibits.

### Benchmarks

Evaluate performance against NumPy/LAPACK:

```bash
just bench
```

This compares in-process Lean 4 execution against single-matrix and vectorized `numpy.linalg.eigh`.

---

## CI/CD Lifecycle & Test Reports

Our Continuous Integration workflow (`.github/workflows/ci.yml`) runs on every pull request targeting `main` or `dev`, and on pushes to `dev`.

### CI Workflow Stages

1. **Spell Check**: Executes `crate-ci/typos` across all repository files.
2. **Lean & Parity Matrix**: Runs on `ubuntu-latest` across multiple Lean 4 toolchains:
   - `leanprover/lean4:v4.27.0` (minimum supported version)
   - `leanprover/lean4:v4.30.0`
   - `leanprover/lean4:v4.34.0-rc2` (pinned development toolchain)
3. **Execution**:
   - Compiles all Lean targets (`lake build`).
   - Runs `just ci`, executing:
     - `lake test`
     - `uv run ruff check parity` & `uv run ruff format --check parity`
     - `uv run pyrefly check parity`
     - `pytest parity/parity.py --junitxml=generated/junit.xml`
4. **Artifact Upload**: Publishes `generated/junit.xml` as a workflow artifact (`parity-junit-report-*`) for test tracking.

Run the entire CI gate locally before opening a PR:

```bash
just ci
```

---

## Code Style & Conventions

### Lean 4 Conventions

- **Language Settings**: We enforce `autoImplicit = false`, `builtinLint = true`, `linter.unusedVariables = true`, `linter.deprecated = true`, and `linter.missingDocs = true` in `lakefile.toml`.
- **Naming Conventions**:
  - Types, structures, and typeclasses: `UpperCamelCase` (e.g., `SymmMat3`, `Vec3`, `Decomposition`).
  - Functions, definitions, and theorems: `lowerCamelCase` (e.g., `eigendecomp`, `maxAbsEntry`, `residual`).
  - Scoped notation: Place Unicode math operators under `open scoped Eig3x3` to avoid polluting downstream namespaces.
- **Explicit Types**: Always provide explicit return types and argument types for public declarations.
- **Float Operations**: Prefer idiomatic Lean float operations (e.g., `x.abs`, `x.sqrt`, `x.frExp`) and use `PowNat` / `^ⁿ` for integer exponentiation.

### Documentation & Docstrings

Because `linter.missingDocs = true` is enabled, all public types, structures, fields, and definitions in `src/` **must** have docstrings (`/-- ... -/`).

- Module documentation (`/-! ... -/`) should describe high-level mathematical concepts and design decisions.
- Use BibTeX citations linking to `references.bib` using the bracket syntax (e.g., `[HaberaZilian2025]` or `[Eberly2014]`).
- When introducing new mathematical sources, add the corresponding entries to `references.bib` and `CITATION.cff`.

### Python Conventions

- **Formatting & Linting**: Enforced via Ruff (`pyproject.toml`, line length 88).
- **Type Annotations**: All Python functions in `parity/` must have complete type signatures validated by Pyrefly (`just typecheck`).
- **No Side Effects**: Differential test scripts must avoid global state mutations and support `pytest-xdist` parallel workers.

---

## Documentation System (`doc-gen4`)

API documentation is generated using `doc-gen4` via the isolated subproject in `docbuild/`.

- **Zero-Dependency Guarantee**: The root `lakefile.toml` must remain dependency-free. `doc-gen4` is strictly required only inside `docbuild/lakefile.toml`.
- **Toolchain Pinning**: The `rev` in `docbuild/lakefile.toml` must match the root `lean-toolchain`. When upgrading the Lean toolchain, bump both files together.

```bash
# Build documentation locally
just docs-dev

# Serve documentation on http://localhost:8000
just docs-serve

# Stop documentation server
just docs-serve-stop
```

> [!CAUTION]
> Do not use `just docs` for local development; it is reserved for the automated GitHub Pages deployment workflow (`.github/workflows/docs.yml`).

---

## Areas for Contribution (Roadmap)

We welcome contributions in several key areas. If you plan to work on a significant feature, please open an issue first for architectural discussion.

1. **Iterative Fallback Solver**:
   - Implement an optional, robust iterative eigensolver (e.g., Jacobi rotations or QR with shifts) as a fallback mechanism for pathological inputs or extreme precision edge cases where closed-form transcendental evaluations lose accuracy.
2. **Complex Hermitian $3 \times 3$ Eigensolver**:
   - Extend or package a dedicated eigensolver for $3 \times 3$ complex Hermitian matrices ($A = A^*$), computing real eigenvalues and unitary eigenvector matrices ($Q^* Q = I$, $\det Q = 1$).
3. **Additional $3 \times 3$ Matrix Decompositions**:
   - $3 \times 3$ Singular Value Decomposition (SVD) based on symmetric eigendecomposition of $A^T A$.
   - Generalized symmetric eigensolvers ($A v = \lambda B v$ for symmetric positive-definite $B$).
   - Polar decomposition ($A = U P$).
   - Fast $3 \times 3$ Cholesky and $LDL^T$ factorizations.
4. **Formal Verification & Bounds Proofs**:
   - Create a separate, opt-in verification companion package (e.g., `Eig3x3.Verify` with Mathlib) formalizing mathematical correctness proofs, condition numbers, and error bounds for the closed-form algorithms.
5. **Platform Portability & Hardware Acceleration**:
   - Performance characterization and optimizations across ARM64 (NEON), x86-64 (AVX2/FMA), WebAssembly (Wasm), and embedded environments.
   - Explorations of SIMD vectorization once native vector primitives stabilize in Lean 4.
6. **Differential Parity Corpus Expansion**:
   - Additional adversarial matrix generators (e.g., Wilkinson-type clustered spectra, random orthogonal similarity orbits $Q D Q^T$).

---

## Contribution Procedure & Git Workflow

### 1. Issue First

For non-trivial enhancements, bug fixes, or new mathematical routines, please open an issue first using one of our GitHub issue templates:
- **Bug Report**
- **Numerical Discrepancy** (for precision/tolerance issues)
- **Feature Request**

### 2. Fork and Branch

1. Fork the repository on GitHub.
2. Clone your fork and create a descriptive branch off `dev` (or `main` for critical patches):
   ```bash
   git checkout -b feat/iterative-fallback
   # or: git checkout -b fix/atan2-quadrant-guard
   ```

### 3. Local Validation

Before committing, ensure your local changes pass all checks:

```bash
# 1. Check spelling
just spell-diff

# 2. Run the full local CI suite
just ci
```

### 4. Commit Guidelines

We follow [Conventional Commits](https://www.conventionalcommits.org/):

- `feat: add 3x3 SVD decomposition`
- `fix: correct angle normalization in planar reduction`
- `perf: optimize matrix-vector multiplication`
- `test: add Wilkinson cluster test cases`
- `docs: update Habera-Zilian 2025 citations`
- `refactor: extract symmetric matrix cofactors`

Keep commit messages concise, descriptive, and focused on the *why* of the change.

### 5. Submitting a Pull Request

1. Push your branch to your fork:
   ```bash
   git push origin feat/iterative-fallback
   ```
2. Open a Pull Request targeting the `dev` branch.
3. Fill out the Pull Request template completely, detailing:
   - Summary of changes and motivation.
   - Algorithmic and mathematical references.
   - Verification steps taken and local test outputs.
   - Confirmation that all numerical gates and CI checks pass.
4. Engage with code review feedback promptly. Once approved and all status checks pass, your contribution will be merged!

---

*Thank you for helping make scientific computing in Lean 4 faster, more robust, and more accessible!*
