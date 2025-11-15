"""Comprehensive test to verify if ChemEx Liouvillians are truly Hermitian.

This is critical because using eigh on non-Hermitian matrices will give
incorrect results.
"""

import sys
from pathlib import Path

import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.configuration.conditions import Conditions
from chemex.models.factory import create_model
from chemex.models.model import model
from chemex.nmr.basis import Basis
from chemex.nmr.liouvillian import LiouvillianIS


def is_hermitian(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if a matrix is Hermitian."""
    return np.allclose(matrix, matrix.conj().T, atol=tol, rtol=tol)


def is_anti_hermitian(matrix: np.ndarray, tol: float = 1e-10) -> bool:
    """Check if a matrix is anti-Hermitian (A = -A†)."""
    return np.allclose(matrix, -matrix.conj().T, atol=tol, rtol=tol)


def check_real_liouvillians() -> None:
    """Test with actual Liouvillian objects as used in ChemEx."""
    print("=" * 70)
    print("Testing REAL Liouvillians from ChemEx")
    print("=" * 70)

    # Set up a 2-state model (most common)
    model.name = "2st"
    model._model = create_model("2st")

    # Test different experiment types
    test_cases = [
        {
            "name": "CEST 15N (ixyz, nh)",
            "basis": Basis(type="ixyz", spin_system="nh"),
            "spin_system": "nh",
        },
        {
            "name": "CPMG 15N (ixyzsz, nh)",
            "basis": Basis(type="ixyzsz", spin_system="nh"),
            "spin_system": "nh",
        },
        {
            "name": "CPMG 1H (ixyzsz, hn)",
            "basis": Basis(type="ixyzsz", spin_system="hn"),
            "spin_system": "hn",
        },
    ]

    for test_case in test_cases:
        print(f"\n{'─' * 70}")
        print(f"Test: {test_case['name']}")
        print(f"{'─' * 70}")

        try:
            basis = test_case["basis"]

            # Create conditions (realistic NMR conditions)
            conditions = Conditions(
                h_larmor_frq=600.0,  # 600 MHz spectrometer
                temperature=298.0,
                p_total=0.0,
                l_total=0.0,
            )

            # Create Liouvillian with realistic parameters
            from chemex.parameters.spin_system import SpinSystem

            spin_system = SpinSystem.from_value(test_case["spin_system"])

            liouv_obj = LiouvillianIS(
                spin_system=spin_system,
                basis=basis,
                conditions=conditions,
            )

            # Set realistic parameters for 2-state exchange
            liouv_obj.par_values = {
                "r2_i_a": 10.0,  # R2 relaxation state A
                "r1_i_a": 2.0,  # R1 relaxation state A
                "r2_i_b": 15.0,  # R2 relaxation state B
                "r1_i_b": 2.5,  # R1 relaxation state B
                "cs_i_a": 120.0,  # Chemical shift state A (ppm)
                "cs_i_b": 118.0,  # Chemical shift state B (ppm)
                "kex_ab": 500.0,  # Exchange rate A→B (s⁻¹)
                "kex_ba": 500.0,  # Exchange rate B→A (s⁻¹)
                "pa": 0.95,  # Population state A
                "pb": 0.05,  # Population state B
            }

            # Set offsets
            liouv_obj.carrier_i = 119.0  # Carrier frequency
            liouv_obj.offset_i = 0.0  # On-resonance

            # Build the Liouvillian
            liouv_obj._build_base_liouvillian()

            # Get the matrix
            liouv_matrix = liouv_obj._l_base

            print(f"Matrix shape: {liouv_matrix.shape}")
            print(f"Matrix dtype: {liouv_matrix.dtype}")

            # Check properties
            is_herm = is_hermitian(liouv_matrix)
            is_anti_herm = is_anti_hermitian(liouv_matrix)

            print(f"\nHermitian: {is_herm}")
            print(f"Anti-Hermitian: {is_anti_herm}")

            # Check eigenvalues
            eigenvalues_eig = np.linalg.eigvals(liouv_matrix)
            are_real = np.allclose(eigenvalues_eig.imag, 0, atol=1e-10)
            print(f"Eigenvalues are real: {are_real}")

            if not are_real:
                print(
                    f"  Max imaginary part: {np.max(np.abs(eigenvalues_eig.imag)):.2e}"
                )

            # Compare eig vs eigh if Hermitian
            if is_herm:
                print("\n✓ Matrix IS Hermitian - eigh is safe to use")
                eigenvalues_eigh = np.linalg.eigvalsh(liouv_matrix)
                print(f"  eigh eigenvalues (sample): {eigenvalues_eigh[:3]}")
            else:
                print("\n✗ Matrix is NOT Hermitian - must use eig!")
                if is_anti_herm:
                    print("  → But it IS anti-Hermitian")
                    print("  → Eigenvalues should be purely imaginary")

        except Exception as e:
            print(f"Error: {e}")
            import traceback

            traceback.print_exc()


def check_with_offset() -> None:
    """Test with RF offset (may break Hermiticity)."""
    print("\n" + "=" * 70)
    print("Testing with RF OFFSET (CEST-like)")
    print("=" * 70)

    model.name = "2st"
    model._model = create_model("2st")

    basis = Basis(type="ixyz", spin_system="nh")
    conditions = Conditions(h_larmor_frq=600.0, temperature=298.0)

    from chemex.parameters.spin_system import SpinSystem

    spin_system = SpinSystem.from_value("nh")

    liouv_obj = LiouvillianIS(
        spin_system=spin_system,
        basis=basis,
        conditions=conditions,
    )

    liouv_obj.par_values = {
        "r2_i_a": 10.0,
        "r1_i_a": 2.0,
        "r2_i_b": 15.0,
        "r1_i_b": 2.5,
        "cs_i_a": 120.0,
        "cs_i_b": 118.0,
        "kex_ab": 500.0,
        "kex_ba": 500.0,
        "pa": 0.95,
        "pb": 0.05,
    }

    # Test with different offsets
    offsets = [0.0, 100.0, 500.0, 1000.0]  # Hz

    for offset in offsets:
        liouv_obj.carrier_i = 119.0
        liouv_obj.offset_i = offset
        liouv_obj._build_base_liouvillian()

        liouv_matrix = liouv_obj._l_base

        is_herm = is_hermitian(liouv_matrix)
        eigenvalues = np.linalg.eigvals(liouv_matrix)
        are_real = np.allclose(eigenvalues.imag, 0, atol=1e-10)

        print(f"\nOffset = {offset} Hz:")
        print(f"  Hermitian: {is_herm}")
        print(f"  Eigenvalues real: {are_real}")

        if not is_herm or not are_real:
            print(f"  ⚠️  WARNING: Not Hermitian or complex eigenvalues!")


def main() -> None:
    """Run comprehensive Hermiticity checks."""
    print("\n" + "=" * 70)
    print("COMPREHENSIVE LIOUVILLIAN HERMITICITY VERIFICATION")
    print("=" * 70)
    print("\nThis checks if ChemEx Liouvillians are truly Hermitian")
    print("to validate the use of np.linalg.eigh optimization.")

    check_real_liouvillians()
    check_with_offset()

    print("\n" + "=" * 70)
    print("Analysis Complete")
    print("=" * 70)
    print("\nConclusion:")
    print("  If ALL tests show Hermitian=True → eigh is safe")
    print("  If ANY test shows Hermitian=False → MUST use eig!")
    print("=" * 70)


if __name__ == "__main__":
    main()
