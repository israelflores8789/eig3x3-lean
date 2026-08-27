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
│   ├── Cli.lean               # lean_exe: parity/bench CLI (JSON in, JSON out)
│   ├── Tests/
│   │   ├── Certificates.lean  # Per-case certificate assertions under machine-epsilon gates
│   │   ├── Golden.lean        # Embedded 50-digit mpmath references
│   │   ├── KnownAnswer.lean   # Exact known-answer assertions
│   │   ├── Regression.lean    # Pinned historical failures (r₁₀ discriminant, near-double path)
│   │   └── Util.lean          # Shared assertions, machine-epsilon gates, case zoo
├── parity/                    # pytest parity suite (Lean exe vs references)
│   ├── gates.py               # shared numerical standard; mirrors Tests/Util.lean
│   ├── gen_cases.py           # case generation
│   ├── mirror.py              # op-for-op bit-exact float64 cross-check port of the Lean Eig3x3 package
│   ├── lean_cli.py            # subprocess bridge to eig3x3-cli; owns JSON schema
│   ├── compare.py             # numpy.linalg.eigvalsh parity plus certificates
│   ├── exact.py               # exact rat-math checks
│   ├── golden.py              # 50-digit mpmath golden vectors; writes golden.json
│   ├── properties.py          # ordering, trace, det, handedness, certificates, bit-exact 2^k scale invariance
│   └── bench.py               # timing and the naive-vs-present discriminant study
├── scripts/                   # dev/CI helper scripts
├── golden.json                # generated golden vectors (generated, versioned, never hand-edited)
├── lakefile.toml              # Lean4 project configuration
├── lean-toolchain             # Lean4 toolchain version
├── lake-manifest.json         # Lean4 project dependencies and metadata
├── pyproject.toml             # Python env (non-packaged) + pytest/ruff config
├── uv.lock
└── justfile                   # command runner
```

### Lean package dependency flow (`Eig3x3/`)
```
Basic ──┬──> Eigenvalues ───┐
        ├──> Eigenvectors ──┴──> Eigendecomp
        └──> Certificates
```

### Python parity layout and dependency flow (`parity/`)
```
layer 0 (no internal deps):  gates    gen_cases    mirror    lean_cli
layer 1:                     exact        → gen_cases
                             golden       → gates, gen_cases, mirror
                             properties   → gates, gen_cases, mirror
layer 2:                     compare      → gates, gen_cases, mirror, lean_cli
                             bench        → compare, gen_cases, mirror
external boundary:  lean_cli ──JSON on stdin/stdout──> eig3x3-cli ──imports──> Eig3x3
artifact flow:      golden.py ──writes──> golden.json ──embedded──> Tests/Golden.lean
```

Data flow: `gen_cases` produces matrices; an implementation (`mirror` now, `lean_cli`→`eig3x3-cli` later) produces decompositions plus certificates; `gates` judges everything in eps·maxAbs units; the four check modules are pass/fail suites; `bench.py` reports and never gates.

Flat by design: nine single-purpose modules, each runnable as a script with a nonzero exit code on failure. Do not restructure into packages or split out a tests/ directory — the layout is documented here and is deliberate.

## Commands

All commands run via `just` from the repo root; `just --list` shows everything.

| Command | What it runs | Pass condition |
|---|---|---|
| `just test` | Lean test suite (`lake test`): known-answer, regression, certificates, golden | exit 0, all assertions `ok` |
| `just exact` | Exact rat-math identity checks | 0 failures |
| `just props` | Property invariants | 0 failures |
| `just parity` | numpy parity vs mirror: eigenvalue error in ε·maxAbs units + certificates under the shared 64ε/16ε gates | all gates held |
| `just golden-check` | validates `eigendecomp` against golden.json | 0 failures |
| `just bench` | timing + naive-vs-present discriminant study | informational only |
| `just ci` | all of the above plus ruff and pyrefly | all green |

- USE `uv run pytest parity -k <expr>` to run a subset of parity tests
- **CLI Smoke Test**: After running `just cli --impl lean`, RUN `echo '{"matrices":[[2.0,2.0,2.0,1.0,0.0,1.0]]}' | .lake/build/bin/eig3x3-cli` and EXPECT eigenvalues 0.58578643762690497, 2.0, 3.4142135623730949 with certificates ≈ 1e-16. IF the result does NOT match expectations, STOP and REPORT the failure.

### Rules
- NEVER loosen a tolerance or edit an expected value to make a failure pass. STOP and REPORT the failure.
- Numerical literals are INTENTIONALLY close (e.g. 1 ULP apart). Do NOT modify these figures; their purpose is to test mathematic float-point execution.
- `golden.json` and the literals in `Tests/Golden.lean` are generated artifacts. NEVER hand-edit. Regenerate via `just golden` to keep `golden.py`'s curated set and `Tests/Golden.lean` in sync. IF `golden.json` is missing, RUN `just golden` to generate the file.
- The parity harness runs against the Python mirror by default. USE the `--impl lean` flag on `just cli` to perform the parity comparison against Lean. This is REQUIRED before release.

## Known Issues
- The Habera–Zilian reference paper contains a typo in: §7, Algorithm 8, `r₁₀`. The corrected form is implemented in `src/Eig3x3/Eigenvalues.lean`. Do NOT "fix" it back.

## Boundaries
- Do NOT add Lake dependencies unless explicitly asked.
- ALWAYS prefer the latest fully **stable** version of the Lean 4 toolchain. Do NOT change `lean-toolchain` without notifying the user. 

## References
- Habera–Zilian (2025) — eigenvalue decomposition
  - title: "Numerically stable evaluation of closed-form expressions for eigenvalues of 3×3 matrices"
  - doi: 10.48550/arXiv.2511.00292
- Eberly (2014) — eigenvector decomposition
  - title: "A Robust Eigensolver for 3×3 Symmetric Matrices"
  - article: <https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf>
