module

public import all Eig3x3.Basic
import Eig3x3.Eigenvalues
import all Eig3x3.Eigenvectors

/-!
# Eig3x3.Eigendecomp — the pipeline driver

`eigendecomp` is the library's primary entry point and the only module that
knows both stages exist. It owns the pipeline policy: zero fast path,
max-abs preconditioning (Eberly's overflow guard), ordered eigenvalues on
the scaled matrix (Habera–Zilian), isolated-eigenvalue-first eigenvector
assembly (Eberly), and eigenvalue rescaling.
-/

namespace Eig3x3

/-- Full eigendecomposition of a real symmetric 3×3 matrix:
    `A = QΛQᵀ = Σᵢ λᵢ cᵢcᵢᵀ` with `l₁ ≤ l₂ ≤ l₃` and `Q` right-handed
    orthonormal.

    Pipeline: zero fast path → max-abs preconditioning (Eberly) → ordered
    eigenvalues on the scaled matrix (Habera–Zilian) → eigenvector assembly
    consuming those eigenvalues (Eberly) → eigenvalue rescaling. -/
public def eigendecomp (A : SymmMat3) : Decomposition :=
  let maxAbs := A.maxAbsEntry
  if maxAbs == 0.0 then
    { eigvals := ⟨0.0, 0.0, 0.0⟩, eigvecs := Mat3.id }
  else
    let s := 1.0 / maxAbs
    let B := A.scale s
    let e := eigvals B
    let Q := eigvecs B e
    -- maxAbs > 0 here, so rescaling preserves the l₁ ≤ l₂ ≤ l₃ ordering.
    { eigvals := e.scale maxAbs, eigvecs := Q }

end Eig3x3
