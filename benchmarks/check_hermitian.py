"""Check if ChemEx Liouvillian matrices are Hermitian.

This determines whether we can use np.linalg.eigh (faster) or must use
np.linalg.eig for eigenvalue decomposition.
"""

import sys
from pathlib import Path

import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.nmr.basis import Basis


def is_hermitian(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if a matrix is Hermitian (equal to its conjugate transpose).

    Args:
        matrix: Matrix to check
        tol: Tolerance for comparison

    Returns:
        True if matrix is Hermitian within tolerance
    """
    return np.allclose(matrix, matrix.conj().T, atol=tol, rtol=tol)


def check_basis_matrices() -> None:
    """Check if basis matrices used to construct Liouvillians are Hermitian."""
    print("=" * 70)
    print("Checking if Basis Matrices are Hermitian")
    print("=" * 70)

    # Test different basis types (from real examples in codebase)
    test_bases = [
        ("ixyz", "nh"),        # 15N CEST
        ("ixyzsz", "nh"),      # 15N CPMG TR
        ("ixyzsz", "hn"),      # 1H CPMG
        ("iz", "nh"),          # Relaxation
    ]

    for basis_type, spin_system in test_bases:
        print(f"\n{'-' * 70}")
        print(f"Basis: type={basis_type}, spin_system={spin_system}")
        print(f"{'-' * 70}")

        try:
            basis = Basis(type=basis_type, spin_system=spin_system)

            hermitian_count = 0
            non_hermitian_count = 0

            print(f"Total matrices: {len(basis.matrices)}")

            for name, matrix in basis.matrices.items():
                if is_hermitian(matrix):
                    hermitian_count += 1
                else:
                    non_hermitian_count += 1
                    print(f"  NON-HERMITIAN: {name} (shape: {matrix.shape})")

            print(f"\nSummary:")
            print(f"  Hermitian: {hermitian_count}")
            print(f"  Non-Hermitian: {non_hermitian_count}")

            if non_hermitian_count == 0:
                print(f"  ✓ All matrices are Hermitian")
            else:
                print(f"  ✗ Some matrices are NOT Hermitian")

        except Exception as e:
            print(f"  Error processing {spin_name}: {e}")


def check_constructed_liouvillian() -> None:
    """Check a fully constructed Liouvillian matrix."""
    print("\n" + "=" * 70)
    print("Checking Constructed Liouvillian")
    print("=" * 70)

    # Create a simple test Liouvillian
    try:
        basis = Basis(type="ixyz", spin_system="nh")

        # Build a test Liouvillian by summing basis matrices with realistic parameters
        par_values = {
            "r2_i_a": 10.0,
            "r1_i_a": 2.0,
            "dw_i_a": 2.0,
            "kex_ab": 500.0,
        }

        # Sum basis matrices (simplified version of _build_base_liouvillian)
        liouvillian = sum(
            (
                basis.matrices[name] * par_values.get(name, 0.0)
                for name in basis.matrices
            ),
            start=np.zeros_like(list(basis.matrices.values())[0]),
        )

        print(f"\nLiouvillian shape: {liouvillian.shape}")
        print(f"Liouvillian dtype: {liouvillian.dtype}")

        if is_hermitian(liouvillian):
            print("✓ Liouvillian IS Hermitian")
            print("  → Can use np.linalg.eigh for 2-3x speedup!")
        else:
            print("✗ Liouvillian is NOT Hermitian")
            print("  → Must use np.linalg.eig (no speedup from eigh)")

            # Check if it's anti-Hermitian (common for Liouvillians)
            if is_hermitian(1j * liouvillian):
                print("  → But it IS anti-Hermitian (L† = -L)")

    except Exception as e:
        print(f"Error: {e}")


def main() -> None:
    """Run all checks."""
    print("\n" + "=" * 70)
    print("ChemEx Liouvillian Hermiticity Check")
    print("=" * 70)
    print("\nThis determines if we can use np.linalg.eigh (fast) or")
    print("must use np.linalg.eig (slower) for eigenvalue decomposition.")

    check_basis_matrices()
    check_constructed_liouvillian()

    print("\n" + "=" * 70)
    print("Analysis Complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
