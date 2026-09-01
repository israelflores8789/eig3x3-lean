## Description

<!-- Provide a concise summary of the changes proposed in this Pull Request. -->

## Motivation & Context

<!-- Why is this change required? What problem does it solve? If it fixes an open issue, please link it here using `Fixes #...` or `Closes #...`. -->

## Type of Change

- [ ] 🐛 Bug fix (non-breaking change which fixes an issue)
- [ ] ✨ New feature (non-breaking change which adds functionality)
- [ ] ⚡ Performance improvement (optimization without altering semantics)
- [ ] 📐 Mathematical / numerical adjustment (must strictly maintain 64ε / 16ε gates)
- [ ] 📝 Documentation update (docstrings, references, or guides)
- [ ] 🧪 Tests / CI update (parity tests, regression cases, or workflow configuration)
- [ ] 🔧 Refactoring / code maintenance

## Mathematical & Algorithmic Details

<!-- If applicable, describe the mathematical derivation, citing relevant literature (e.g., Habera & Zilian 2025, Eberly 2014) or formulas. -->

## Verification & Testing

Please confirm the verification steps performed:

- [ ] `just test` (Lean test suite: KnownAnswer, Golden, Properties, Regression, Certificates)
- [ ] `just lint` (Ruff formatting & lint checks)
- [ ] `just typecheck` (Pyrefly static type checks)
- [ ] `just parity-ci` (Python differential parity suite under 64ε / 16ε gates)
- [ ] `just spell-diff` (No spelling violations)
- [ ] `just ci` passes cleanly end-to-end locally

### Test Output / Evidence

<!-- Paste relevant command outputs, benchmark results, or certificate outputs if applicable. -->

```text

```

## Checklist

- [ ] My code adheres to the project's style conventions and Lean community standards.
- [ ] I have maintained the zero-dependency rule for the root `Eig3x3` package (no Mathlib/Batteries in `lakefile.toml`).
- [ ] I have added docstrings (`/-- ... -/`) for all new public declarations (`linter.missingDocs = true`).
- [ ] I have not loosened any numerical tolerance gates or edited expected test values to force a pass.
- [ ] If new citations are introduced, they have been added to `references.bib` and `CITATION.cff`.
