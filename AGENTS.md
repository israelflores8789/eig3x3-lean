# Eig3x3 Lean - Agent Context

## Purpose
Pure Lean4 `Float` type computable library for the eigendecomposition of 3×3 real symmetric matrices using eigenvalue decomposition via Habera-Zilian (2025), with the §7 stabilized discriminant, and eigenvector decomposition via Eberly (2014) null-space construction. No FFI, no Mathlib, no dependencies. Designed to be a drop-in Lake dependency for other projects.

## Repository Structure
```
eig3x3-lean/
├── src/                       # ALL Lean source
│   ├── Eig3x3.lean            # root module, re-exports the library
│   ├── Eig3x3/
│   │   ├── Basic.lean         # Type definitions and primitives
│   │   ├── Eigenvalues.lean   # Habera–Zilian port with improvements
│   │   ├── Eigenvectors.lean  # Eberly port with improvements
│   │   ├── Eigendecomp.lean   # Perform eigendecomposition
│   │   └── Certificates.lean  # Runtime calculation assurances
├── tests/
│   ├── Main.lean              # lean_exe: test driver (`lake test`)
│   ├── Cli.lean               # lean_exe: eig3x3_cli, parity/bench JSON-driven CLI
│   ├── JsonMiniReader.lean    # shared minimal JSON reader (Cli + Golden)
│   ├── Tests/
│   │   ├── Certificates.lean  # Per-case certificate assertions under machine-epsilon gates
│   │   ├── Golden.lean        # Assertions against "golden" reference vectors
│   │   ├── KnownAnswer.lean   # Exact known-answer assertions
│   │   ├── Regression.lean    # Pinned historical failures (r₁₀ discriminant, near-double path)
│   │   └── Util.lean          # Shared assertions, machine-epsilon gates, case zoo
├── parity/                    # pytest parity suite (Lean exe vs references)
│   ├── gates.py               # shared numerical standard; mirrors Tests/Util.lean
│   ├── gen_cases.py           # case generation
│   ├── lean_cli.py            # subprocess bridge to eig3x3_cli; owns JSON schema
│   ├── parity.py              # pytest test harness & numpy/LAPACK parity test
│   ├── errata.py              # exhibits of errata found during implementation
│   ├── golden.py              # 50-digit mpmath golden vectors generator; writes golden.json
│   ├── conftest.py            # pytest configuration, auto-scaffolding fixtures & CLI options
│   └── bench.py               # performance benchmark of Lean against numpy/LAPACK
├── scripts/                   # dev/CI helper scripts
├── golden.json                # generated golden vectors (generated, versioned, never hand-edited)
├── lakefile.toml              # Lean4 project configuration
├── lean-toolchain             # Lean4 toolchain version
├── lake-manifest.json         # Lean4 project dependencies and metadata
├── pyproject.toml             # Python env (non-packaged) + pytest/ruff config
├── uv.lock
└── justfile                   # command runner
```

### Lean source package dependency flow (`Eig3x3/`)
```
Basic ──┬──> Eigenvalues ───┐
        ├──> Eigenvectors ──┴──> Eigendecomp
        └──> Certificates
```

### Lean test suite dependency flow (`tests/`)
```
Eig3x3 ──┬──> Util ──┬──> KnownAnswer ───┐
         │           ├──> Golden ────────┤
         │           ├──> Properties ────┤
         │           ├──> Regression ────┤
         │           └──> Certificates ──┴──> Main  (lean_exe) <──┬── JsonMiniReader
         ├──────────────────────────────────> Cli   (lean_exe) <──┘
         └──────────────────────────────────> Bench (lean_exe)
```

### Python parity test dependency flow (`parity/`)
```
Parity Test:
                                         generated/golden.json ──┐
  gen_cases.py ──> lean_cli.py ──┬──> eig3x3_cli (Lean binary) ──┼──> parity.py
                                 │                    gates.py ──┘
                                 └──> golden.py ────────────────────> conftest.py (pytest)

Generate "Golden" Vectors:
  golden.py (generate vectors) ──> generated/golden.json ──> Tests/Golden.lean

Benchmark:
  bench.py  ──> eig3x3_cli (Lean binary)
```

## Commands
All commands run via `just` from the repo root; `just --list` shows everything.

| Command | What it runs | Pass condition |
|---|---|---|
| `just test` | Lean test suite (`lake test`): known-answer, regression, certificates, golden | exit 0, all assertions `ok` |
| `just errata` | errata exhibits during algorithm implementation | 0 failures |
| `just parity` | executes python parity test suite: eigenvalue error in ε·maxAbs units + certificates under the shared 64ε/16ε gates | all gates held |
| `just bench` | performance benchmark of Lean binary against numpy/LAPACK | informational only |
| `just ci` | all of the above plus ruff and pyrefly | all green |

- USE `just build-all` to build the Lean source, bench, and CLI binaries.
- **CLI Smoke Test**: After running `just build-cli`, RUN `echo '{"matrices":[[2.0,2.0,2.0,1.0,0.0,1.0]]}' | .lake/build/bin/eig3x3_cli` and EXPECT eigenvalues 0.58578643762690497, 2.0, 3.4142135623730949 with certificates ≈ 1e-16. IF the result does NOT match expectations, STOP and REPORT the failure.
- USE `just spell-diff` to check for spelling errors FIRST, VERIFY the errors presented make semantic sense AND are not correctly spelled terms (e.g. author's names like "Michal Habera" etc), THEN USE `just spell-fix` to correct spelling errors package-wide IF safe to do so. IF an error presented is identified as a correctly spelled term, STOP and PROPOSE to the user a change to the `_typos.toml` config file (create ONLY if missing and required). Do NOT use `just spell`.

### Rules
- NEVER loosen a tolerance or edit an expected value to make a failure pass. STOP and REPORT the failure.
- Numerical literals are INTENTIONALLY close (e.g. 1 ULP apart). Do NOT modify these figures; their purpose is to test mathematic float-point execution.
- the directory `generated/` contains generated artifacts (e.g. `golden.json`). NEVER edit the contents inside `generated/`.
  - USE `just golden` to generate or update `generated/golden.json`.
- `parity/gates.py` mirrors `Tests/Util.lean`; IF requested by the user, modify them together or STOP and REPORT a conflict.
- Do NOT modify `parity/errata.py`.
- USE `just <cmd>` for all tests. Do NOT perform manual command line equivalents UNLESS you believe the `just` command to be stale. IF you believe a dedicated `just` command is stale, REPORT to the user your proposed update.

## Known Issues
- The Habera–Zilian reference paper contains a typo in: §7, Algorithm 8, `r₁₀`. The corrected form is implemented in `src/Eig3x3/Eigenvalues.lean`. Do NOT "fix" it back.

## Boundaries
- Do NOT add Lake dependencies unless explicitly asked.
- ALWAYS prefer the latest fully **stable** version of the Lean 4 toolchain. Do NOT change `lean-toolchain` without notifying the user.

## Documentation (doc-gen4)
- API docs are generated with **doc-gen4** from the nested `docbuild/` Lake package. The root `lakefile.toml` must remain dependency-free. NEVER add `doc-gen4` to it. The `docbuild/` package is inert for downstream consumers.
- The doc-gen4 `rev` in `docbuild/lakefile.toml` is pinned to match the root `lean-toolchain` stable release (`v4.x`). Bump the two *together*, never independently.
- Do NOT build documentation UNLESS asked.
- Do NOT build documentation as a matter of routine when running tests or verifying source package changes.

### Commands
- `just docs-dev` — build docs into `docbuild/.lake/build/doc/`
- `just docs-update` — the ONLY way doc-gen4 gets updated. Do NOT update every build unless requested.
- Do NOT use `just docs` — reserved for a separate docs CI/CD GitHub workflow (`docs.yml`).

### Citations
- `references.bib` is the docs bibliography. Cite in docstrings as `[Key]` or `[text][Key]` to allow doc-gen4 to render citations as links.
- `Eigenvalues.lean` must cite [HaberaZilian2025]; `Eigenvectors.lean` must cite [Eberly2014]. New algorithmic sources get a bib entry first, then citations.
- Keep `references.bib` in sync with `CITATION.cff`'s `references:` list.

## References
- Habera–Zilian (2025) — eigenvalue decomposition
  - title: "Numerically stable evaluation of closed-form expressions for eigenvalues of 3×3 matrices"
  - doi: 10.48550/arXiv.2511.00292
- Eberly (2014) — eigenvector decomposition
  - title: "A Robust Eigensolver for 3×3 Symmetric Matrices"
  - article: <https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf>
