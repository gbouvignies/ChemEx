"""Benchmark Liouvillian construction performance.

Tests current vs vectorized implementation of _build_base_liouvillian.
"""

import sys
import time
from pathlib import Path

import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def benchmark_current_implementation():
    """Profile the current _build_base_liouvillian implementation."""
    from chemex.nmr.basis import Basis
    from chemex.nmr.liouvillian import LiouvillianIS
    from chemex.parameters.spin_system import SpinSystem

    # Create a realistic Liouvillian
    spin_system = SpinSystem.from_name("10N-HN")
    basis = Basis(type="ixyz", spin_system="nh")

    # Simulate conditions
    class MockConditions:
        h_larmor_frq = 600.0
        temperature = 298.0
        p_total = 1e-3
        l_total = 1e-3
        d2o = 0.1

    conditions = MockConditions()
    liouvillian = LiouvillianIS(spin_system, basis, conditions)

    # Typical parameter values for fitting
    par_values = {
        "r2_i_a": 10.0,
        "r2_i_b": 15.0,
        "r1_i_a": 1.5,
        "r1_i_b": 2.0,
        "kab": 200.0,
        "kba": 50.0,
        "cs_i_a": 0.0,
        "cs_i_b": 3.0,
    }

    print("=" * 70)
    print("Benchmarking Liouvillian Construction")
    print("=" * 70)
    print(f"Basis matrices: {len(liouvillian.basis.matrices)}")
    print(f"Matrix size: {liouvillian.size}x{liouvillian.size}")
    print(f"Parameters: {len(par_values)}")

    # Benchmark current implementation
    iterations = 10000
    start = time.perf_counter()
    for _ in range(iterations):
        liouvillian.update(par_values)
    elapsed = time.perf_counter() - start

    print(f"\nCurrent implementation:")
    print(f"  {iterations} iterations: {elapsed:.4f}s")
    print(f"  Time per call: {elapsed/iterations*1000:.4f}ms")
    print(f"  Calls per second: {iterations/elapsed:.0f}")

    return liouvillian, par_values, elapsed / iterations


def benchmark_vectorized_implementation():
    """Test a vectorized version of _build_base_liouvillian."""
    from chemex.nmr.basis import Basis
    from chemex.nmr.liouvillian import LiouvillianIS
    from chemex.parameters.spin_system import SpinSystem

    # Create same setup
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

    par_values = {
        "r2_i_a": 10.0,
        "r2_i_b": 15.0,
        "r1_i_a": 1.5,
        "r1_i_b": 2.0,
        "kab": 200.0,
        "kba": 50.0,
        "cs_i_a": 0.0,
        "cs_i_b": 3.0,
    }

    # Pre-compute stacked matrices (done once during init)
    matrix_names = list(liouvillian.basis.matrices.keys())
    stacked_matrices = np.stack(
        [liouvillian.basis.matrices[name] for name in matrix_names], axis=0
    )
    print(f"\nVectorized setup:")
    print(f"  Stacked matrices shape: {stacked_matrices.shape}")
    print(f"  Number of matrices: {len(matrix_names)}")

    # Vectorized update function
    def vectorized_build_base(par_values_dict):
        param_array = np.array(
            [par_values_dict.get(name, 0.0) for name in matrix_names]
        )
        return np.einsum("i,ijk->jk", param_array, stacked_matrices)

    # Verify correctness first
    liouvillian.update(par_values)
    original_result = liouvillian._l_base.copy()
    vectorized_result = vectorized_build_base(par_values)

    diff = np.max(np.abs(original_result - vectorized_result))
    print(f"\n  Numerical difference: {diff:.2e}")
    if diff < 1e-10:
        print("  ✓ Results match!")
    else:
        print("  ✗ Results differ!")
        return None, None

    # Benchmark vectorized version
    iterations = 10000
    start = time.perf_counter()
    for _ in range(iterations):
        _ = vectorized_build_base(par_values)
    elapsed = time.perf_counter() - start

    print(f"\nVectorized implementation:")
    print(f"  {iterations} iterations: {elapsed:.4f}s")
    print(f"  Time per call: {elapsed/iterations*1000:.4f}ms")
    print(f"  Calls per second: {iterations/elapsed:.0f}")

    return vectorized_build_base, elapsed / iterations


def benchmark_tensordot_implementation():
    """Test tensordot version (alternative to einsum)."""
    from chemex.nmr.basis import Basis
    from chemex.nmr.liouvillian import LiouvillianIS
    from chemex.parameters.spin_system import SpinSystem

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

    par_values = {
        "r2_i_a": 10.0,
        "r2_i_b": 15.0,
        "r1_i_a": 1.5,
        "r1_i_b": 2.0,
        "kab": 200.0,
        "kba": 50.0,
        "cs_i_a": 0.0,
        "cs_i_b": 3.0,
    }

    matrix_names = list(liouvillian.basis.matrices.keys())
    stacked_matrices = np.stack(
        [liouvillian.basis.matrices[name] for name in matrix_names], axis=0
    )

    def tensordot_build_base(par_values_dict):
        param_array = np.array(
            [par_values_dict.get(name, 0.0) for name in matrix_names]
        )
        return np.tensordot(param_array, stacked_matrices, axes=([0], [0]))

    iterations = 10000
    start = time.perf_counter()
    for _ in range(iterations):
        _ = tensordot_build_base(par_values)
    elapsed = time.perf_counter() - start

    print(f"\nTensordot implementation:")
    print(f"  {iterations} iterations: {elapsed:.4f}s")
    print(f"  Time per call: {elapsed/iterations*1000:.4f}ms")
    print(f"  Calls per second: {iterations/elapsed:.0f}")

    return tensordot_build_base, elapsed / iterations


if __name__ == "__main__":
    liouvillian, par_values, current_time = benchmark_current_implementation()
    vectorized_func, vectorized_time = benchmark_vectorized_implementation()
    tensordot_func, tensordot_time = benchmark_tensordot_implementation()

    if vectorized_time and tensordot_time:
        print("\n" + "=" * 70)
        print("PERFORMANCE COMPARISON")
        print("=" * 70)
        print(f"Current (sum with generator): {current_time*1000:.4f}ms")
        print(f"Vectorized (einsum):          {vectorized_time*1000:.4f}ms")
        print(f"Vectorized (tensordot):       {tensordot_time*1000:.4f}ms")

        speedup_einsum = current_time / vectorized_time
        speedup_tensordot = current_time / tensordot_time

        print(f"\nSpeedup (einsum):    {speedup_einsum:.2f}x")
        print(f"Speedup (tensordot): {speedup_tensordot:.2f}x")

        if speedup_einsum > speedup_tensordot:
            print(f"\n✓ einsum is faster")
        else:
            print(f"\n✓ tensordot is faster")
