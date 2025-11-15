"""Validate vectorized Liouvillian construction.

Compare old (sum-based) vs new (einsum-based) implementation for:
1. Numerical accuracy
2. Performance in realistic scenarios
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def old_build_base_liouvillian(liouvillian):
    """Original implementation using sum with generator."""
    return sum(
        (
            liouvillian.basis.matrices[name] * liouvillian.par_values.get(name, 0.0)
            for name in liouvillian.basis.matrices
        ),
        start=np.array(0.0),
    )


def test_numerical_accuracy():
    """Test that vectorized implementation gives identical results."""
    from chemex.nmr.basis import Basis
    from chemex.nmr.liouvillian import LiouvillianIS
    from chemex.parameters.spin_system import SpinSystem

    print("=" * 70)
    print("NUMERICAL ACCURACY TEST")
    print("=" * 70)

    spin_system = SpinSystem.from_name("10N-HN")
    basis = Basis(type="ixyz", spin_system="nh")

    class MockConditions:
        h_larmor_frq = 600.0
        temperature = 298.0
        p_total = 1e-3
        l_total = 1e-3
        d2o = 0.1

    conditions = MockConditions()
    liouvillian = LiouvillianIS(spin_system, basis, conditions)

    # Test with various parameter sets
    test_cases = [
        {
            "r2_i_a": 10.0,
            "r2_i_b": 15.0,
            "r1_i_a": 1.5,
            "r1_i_b": 2.0,
            "kab": 200.0,
            "kba": 50.0,
            "cs_i_a": 0.0,
            "cs_i_b": 3.0,
        },
        {
            "r2_i_a": 5.0,
            "r2_i_b": 25.0,
            "r1_i_a": 2.5,
            "r1_i_b": 3.0,
            "kab": 500.0,
            "kba": 125.0,
            "cs_i_a": -1.0,
            "cs_i_b": 5.0,
        },
        {
            "r2_i_a": 20.0,
            "r2_i_b": 8.0,
            "r1_i_a": 0.8,
            "r1_i_b": 1.2,
            "kab": 100.0,
            "kba": 25.0,
            "cs_i_a": 2.0,
            "cs_i_b": -2.0,
        },
    ]

    all_pass = True
    for i, params in enumerate(test_cases):
        liouvillian.par_values = params
        old_result = old_build_base_liouvillian(liouvillian)
        liouvillian._build_base_liouvillian()
        new_result = liouvillian._l_base

        diff = np.max(np.abs(old_result - new_result))
        relative_diff = diff / (np.max(np.abs(old_result)) + 1e-20)

        print(f"\nTest case {i+1}:")
        print(f"  Max absolute difference: {diff:.2e}")
        print(f"  Max relative difference: {relative_diff:.2e}")

        if diff < 1e-14:
            print("  ✓ PASS (identical within machine precision)")
        else:
            print("  ✗ FAIL")
            all_pass = False

    return all_pass


def test_realistic_workflow():
    """Simulate a realistic CPMG fitting workflow."""
    from chemex.nmr.basis import Basis
    from chemex.nmr.liouvillian import LiouvillianIS
    from chemex.parameters.spin_system import SpinSystem

    print("\n" + "=" * 70)
    print("REALISTIC WORKFLOW SIMULATION")
    print("=" * 70)

    # Simulate 96 profiles × 100 iterations = 9,600 Liouvillian updates
    num_profiles = 96
    num_iterations = 100
    total_updates = num_profiles * num_iterations

    print(f"Simulating {num_profiles} profiles × {num_iterations} iterations")
    print(f"Total Liouvillian updates: {total_updates}")

    # Create Liouvillians for multiple profiles
    liouvillians = []
    for i in range(num_profiles):
        spin_system = SpinSystem.from_name(f"{10+i}N-HN")
        basis = Basis(type="ixyz", spin_system="nh")

        class MockConditions:
            h_larmor_frq = 600.0
            temperature = 298.0
            p_total = 1e-3
            l_total = 1e-3
            d2o = 0.1

        conditions = MockConditions()
        liouvillian = LiouvillianIS(spin_system, basis, conditions)
        liouvillians.append(liouvillian)

    # Simulate parameter updates during fitting
    rng = np.random.default_rng(42)
    base_params = {
        "r2_i_a": 10.0,
        "r2_i_b": 15.0,
        "r1_i_a": 1.5,
        "r1_i_b": 2.0,
        "kab": 200.0,
        "kba": 50.0,
        "cs_i_a": 0.0,
        "cs_i_b": 3.0,
    }

    # Time the workflow
    start = time.perf_counter()
    for _iteration in range(num_iterations):
        # Simulate parameter changes (like in fitting)
        params = {
            k: v * (1 + 0.01 * rng.standard_normal()) for k, v in base_params.items()
        }
        # Update all profiles with new parameters
        for liouvillian in liouvillians:
            liouvillian.update(params)
    elapsed = time.perf_counter() - start

    print(f"\nTotal time for {total_updates} updates: {elapsed:.4f}s")
    print(f"Time per update: {elapsed/total_updates*1000:.4f}ms")
    print(f"Updates per second: {total_updates/elapsed:.0f}")

    # Compare to baseline (0.024ms from benchmark)
    baseline_time_per_update = 0.024  # ms
    expected_baseline_total = baseline_time_per_update * total_updates / 1000  # seconds
    speedup = expected_baseline_total / elapsed

    print(f"\nExpected baseline time: {expected_baseline_total:.4f}s")
    print(f"Measured speedup: {speedup:.2f}x")

    if speedup > 2.0:
        print("✓ Significant speedup achieved!")
    elif speedup > 1.5:
        print("✓ Moderate speedup achieved")
    else:
        print("⚠ Limited speedup - may need investigation")


if __name__ == "__main__":
    accuracy_ok = test_numerical_accuracy()
    if accuracy_ok:
        print("\n✓ All numerical accuracy tests passed!")
    else:
        print("\n✗ Numerical accuracy tests failed!")
        sys.exit(1)

    test_realistic_workflow()

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)
