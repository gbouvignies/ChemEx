"""Direct test: Are the Liouvillian matrices passed to calculate_propagators Hermitian?

This intercepts actual Liouvillian matrices during a real ChemEx run.
"""

import sys
from pathlib import Path

import numpy as np

# Patch calculate_propagators to capture matrices
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

captured_matrices = []


def is_hermitian(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if a matrix is Hermitian."""
    if matrix.size == 0:
        return True
    return np.allclose(matrix, matrix.conj().T, atol=tol, rtol=tol)


def is_anti_hermitian(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if a matrix is anti-Hermitian."""
    if matrix.size == 0:
        return True
    return np.allclose(matrix, -matrix.conj().T, atol=tol, rtol=tol)


# Monkey-patch to capture matrices
from chemex.nmr import spectrometer

original_calculate_propagators = spectrometer.calculate_propagators


def patched_calculate_propagators(liouv, delays, *, dephasing=False):
    """Capture the Liouvillian matrix and call original."""
    global captured_matrices
    captured_matrices.append(liouv.copy())
    return original_calculate_propagators(liouv, delays, dephasing=dephasing)


spectrometer.calculate_propagators = patched_calculate_propagators

# Now run a real ChemEx simulation
print("=" * 70)
print("Testing Liouvillian Hermiticity During Real ChemEx Run")
print("=" * 70)
print("\nRunning CPMG simulation to capture Liouvillian matrices...")

import glob
import subprocess
import tempfile

example_dir = Path(__file__).parent.parent / "examples/Experiments/CPMG_15N_IP_0013"

if example_dir.exists():
    # Get all experiment TOML files
    experiment_files = list((example_dir / "Experiments").glob("*.toml"))

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = ["chemex", "simulate"]
        for exp_file in experiment_files:
            cmd.extend(["-e", str(exp_file)])
        cmd.extend(["-p", str(example_dir / "Parameters/parameters.toml")])
        cmd.extend(["-o", tmpdir])

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
            cwd=example_dir.parent.parent,
        )

        if result.returncode != 0:
            print("Error running simulation:")
            print(result.stderr)
        else:
            print(f"✓ Simulation completed successfully")
            print(f"✓ Captured {len(captured_matrices)} Liouvillian matrices")

            # Analyze captured matrices
            print("\n" + "=" * 70)
            print("Analysis of Captured Liouvillian Matrices")
            print("=" * 70)

            hermitian_count = 0
            anti_hermitian_count = 0
            neither_count = 0

            unique_shapes = set()
            sample_matrices = []

            for i, matrix in enumerate(captured_matrices[:100]):  # Check first 100
                unique_shapes.add(matrix.shape)

                is_herm = is_hermitian(matrix)
                is_anti = is_anti_hermitian(matrix)

                if is_herm:
                    hermitian_count += 1
                    if len(sample_matrices) < 3:
                        sample_matrices.append(("Hermitian", matrix, i))
                elif is_anti:
                    anti_hermitian_count += 1
                    if len(sample_matrices) < 3:
                        sample_matrices.append(("Anti-Hermitian", matrix, i))
                else:
                    neither_count += 1
                    if len(sample_matrices) < 3:
                        sample_matrices.append(("Neither", matrix, i))

            total = min(100, len(captured_matrices))
            print(f"\nMatrices analyzed: {total}")
            print(f"Unique shapes: {unique_shapes}")
            print(f"\nResults:")
            print(f"  Hermitian:     {hermitian_count:4d} ({hermitian_count/total*100:.1f}%)")
            print(
                f"  Anti-Hermitian: {anti_hermitian_count:4d} ({anti_hermitian_count/total*100:.1f}%)"
            )
            print(f"  Neither:        {neither_count:4d} ({neither_count/total*100:.1f}%)")

            # Check eigenvalues for sample matrices
            print("\n" + "─" * 70)
            print("Sample Matrix Analysis")
            print("─" * 70)

            for matrix_type, matrix, idx in sample_matrices[:3]:
                print(f"\nMatrix #{idx} ({matrix_type}):")
                print(f"  Shape: {matrix.shape}")
                print(f"  Dtype: {matrix.dtype}")

                eigenvalues = np.linalg.eigvals(matrix)
                are_real = np.allclose(eigenvalues.imag, 0, atol=1e-10)

                print(f"  Eigenvalues real: {are_real}")
                if not are_real:
                    max_imag = np.max(np.abs(eigenvalues.imag))
                    print(f"  Max imaginary part: {max_imag:.2e}")

                print(f"  Sample eigenvalues:")
                for j, eig in enumerate(eigenvalues[:5]):
                    if abs(eig.imag) < 1e-10:
                        print(f"    λ{j} = {eig.real:.6f}")
                    else:
                        print(f"    λ{j} = {eig.real:.6f} + {eig.imag:.6f}i")

            # Conclusion
            print("\n" + "=" * 70)
            print("CONCLUSION")
            print("=" * 70)

            if hermitian_count == total:
                print("✓ ALL matrices are Hermitian")
                print("  → Using np.linalg.eigh is SAFE and CORRECT")
            elif anti_hermitian_count == total:
                print("⚠ ALL matrices are anti-Hermitian")
                print("  → Eigenvalues are purely imaginary")
                print("  → np.linalg.eigh will give INCORRECT results")
                print("  → MUST use np.linalg.eig")
            elif hermitian_count > 0 and neither_count > 0:
                print("⚠ MIXED: Some Hermitian, some not")
                print(
                    f"  → {hermitian_count}/{total} are Hermitian, {neither_count}/{total} are not"
                )
                print("  → Current implementation (check + fallback) is CORRECT")
            else:
                print("✗ NO matrices are Hermitian or anti-Hermitian")
                print("  → np.linalg.eigh will give INCORRECT results")
                print("  → MUST revert to np.linalg.eig only")

else:
    print(f"Example directory not found: {example_dir}")
    print("Cannot run test without example data")
