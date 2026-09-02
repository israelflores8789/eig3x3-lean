# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-09-01

### Added

- Closed-form eigendecomposition of 3x3 real symmetric matrices in pure Lean 4
  over `Float`, with no Mathlib dependency
- Eigenvalues via the Habera-Zilian (2025) closed-form method
- Eigenvector construction following Eberly's (2014) r0 × r1 cross-product
  convention
- Core types: `SymmMat3`, `Vec3`, `Mat3`, and `Decomposition` (with
  `eigvals` and `eigvecs` fields)
- Numerical certificate checks at maximum justified machine precision:
  residual, orthonormality, reconstruction, determinant, and eigenvalue
  ordering
- Python parity test suite validating results against NumPy/LAPACK
- Supported on Lean toolchains v4.27.0, v4.30.0, and v4.34.0-rc2

[1.0.0]: https://github.com/israelflores8789/eig3x3-lean/releases/tag/v1.0.0
