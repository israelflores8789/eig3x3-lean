# Eig3x3 Lean - Agent Context

## Purpose
Pure Lean4 `Float` type computable library for the eigendecomposition of 3×3 real symmetric matrices using eigenvalue decomposition via Habera-Zilian (2025), with the §7 stabilized discriminant, and eigenvector decomposition via Eberly (2014) null-space construction. No FFI, no Mathlib, no dependencies. Designed to be a drop-in Lake dependency for other projects.

## Repository Structure
```
eig3x3-lean/
├── src/                       # ALL Lean source
│   ├── Eig3x3.lean            # root module, re-exports the library
│   ├── Eig3x3/
│   │   ├── Eigenvalues.lean   # Habera–Zilian port with improvements
│   │   ├── Eigenvectors.lean  # Eberly port with improvements
│   │   └── ...
│   └── Main.lean              # lean_exe: parity/bench CLI driver
├── parity/                    # pytest parity suite (Lean exe vs references)
├── benchmarks/                # performance harness
├── scripts/                   # dev/CI helper scripts
├── testdata/                  # generated golden vectors — never hand-edit
├── lakefile.toml              # Lean4 project configuration
├── lean-toolchain             # Lean4 toolchain version
├── lake-manifest.json         # Lean4 project dependencies and metadata
├── pyproject.toml             # Python env (non-packaged) + pytest/ruff config
├── uv.lock
└── justfile                   # command runner
```

## Known Issues
- The Habera–Zilian reference paper contains a typo in: §7, Algorithm 8, `r₁₀`. The corrected form
  is implemented in `src/Eig3x3/Eigenvalues.lean`. Do NOT "fix" it back.

## Boundaries
- Do NOT and Lake dependencies unless explicitly asked.
- ALWAYS prefer the latest fully **stable** version of the Lean 4 toolchain. Do NOT change `lean-toolchain` without notifying the user. 

## References
- Habera–Zilian (2025) — eigenvalue decomposition
  - title: "Numerically stable evaluation of closed-form expressions for eigenvalues of 3×3 matrices"
  - doi: 10.48550/arXiv.2511.00292
  - repo: <https://github.com/michalhabera/eig3x3>
- Eberly (2014) — eigenvector decomposition
  - title: "A Robust Eigensolver for 3×3 Symmetric Matrices"
  - book: <https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf>
  - repo: <https://github.com/davideberly/GeometricTools>
