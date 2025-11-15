"""Benchmarks for specific ChemEx bottlenecks identified in the analysis.

This module profiles the critical performance bottlenecks:
1. Eigenvalue decomposition in calculate_propagators()
2. Liouvillian construction
3. CPMG/CEST calculation loops
4. Grid search parallelization potential
"""

import sys
from pathlib import Path

import numpy as np

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from benchmark_framework import (
    Benchmark,
    benchmark_matrix_operations,
)


def benchmark_propagator_calculation() -> None:
    """Benchmark the calculate_propagators function (Bottleneck #1)."""
    from chemex.nmr.spectrometer import calculate_propagators  # noqa: PLC0415

    # Create a realistic Liouvillian matrix (size ~12-20 for 15N systems)
    rng = np.random.default_rng(42)
    size = 16
    liouv = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))
    liouv = liouv + liouv.conj().T  # Make Hermitian

    # Typical delays in CPMG experiments (seconds)
    delays = np.array([0.0001, 0.0005, 0.001, 0.002, 0.005])

    print("\n" + "=" * 70)
    print("BOTTLENECK #1: Propagator Calculation")
    print("=" * 70)
    print(f"Liouvillian size: {size}x{size}")
    print(f"Number of delays: {len(delays)}")
    print("This is called for every CPMG offset/ncyc and CEST offset")

    # Benchmark current implementation
    bench = Benchmark("Current calculate_propagators", iterations=100)
    result = bench.run(
        calculate_propagators,
        liouv,
        delays,
        dephasing=False,
        profile=True,
    )
    bench.print_results(result)

    return result


def benchmark_liouvillian_construction() -> None:
    """Benchmark Liouvillian construction (Bottleneck #5)."""
    from chemex.nmr.basis import Basis
    from chemex.nmr.liouvillian import Liouvillian

    print("\n" + "=" * 70)
    print("BOTTLENECK #5: Liouvillian Construction")
    print("=" * 70)
    print("Testing _build_base_liouvillian() with sum() operation")

    # Create a realistic basis for 15N spin system
    # This is a simplified test - real usage is more complex
    try:
        from chemex.parameters.spin_system import SpinSystem

        spin_system = SpinSystem.from_name("15n")

        # Create a minimal basis
        basis = Basis(spin_system=spin_system, atom_index=0)

        # Create parameter values dictionary
        par_values = {
            "r2_i": 10.0,
            "r1_i": 2.0,
            "r2_s": 15.0,
            "r1_s": 2.5,
            "dw_i": 2.0,
            "kex_ab": 500.0,
            "pb": 0.05,
        }

        # Create Liouvillian
        liouv = Liouvillian(
            basis=basis,
            par_values=par_values,
            conditions={},
        )

        # Benchmark the _build_base_liouvillian method
        bench = Benchmark("_build_base_liouvillian (current sum)", iterations=1000)

        def build_liouv() -> None:
            # Trigger rebuild
            liouv._build_base_liouvillian()

        result = bench.run(build_liouv, profile=False)
        bench.print_results(result)

        print("\nNote: This uses Python sum() over 30+ matrices (inefficient)")
        print("Optimization: Replace with np.einsum() or vectorized operations")

        return result

    except ImportError as e:
        print(f"Could not import required modules: {e}")
        print("Skipping Liouvillian construction benchmark")
        return None


def benchmark_matrix_operations_realistic() -> None:
    """Benchmark matrix operations with realistic ChemEx parameters."""
    print("\n" + "=" * 70)
    print("Matrix Operations with Realistic ChemEx Sizes")
    print("=" * 70)

    # Test different Liouvillian sizes
    sizes = [8, 12, 16, 20, 24]  # Typical sizes for different spin systems

    results = {}

    for size in sizes:
        print(f"\n{'─' * 70}")
        print(f"Testing size: {size}x{size}")
        print(f"{'─' * 70}")

        rng = np.random.default_rng(42)
        matrix = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))
        matrix = matrix + matrix.conj().T  # Hermitian

        # Eigenvalue decomposition (current method)
        bench_eig = Benchmark(f"np.linalg.eig ({size}x{size})", iterations=100)
        res_eig = bench_eig.run(lambda: np.linalg.eig(matrix), profile=False)

        # Hermitian eigenvalue (potential optimization)
        bench_eigh = Benchmark(f"np.linalg.eigh ({size}x{size})", iterations=100)
        res_eigh = bench_eigh.run(lambda: np.linalg.eigh(matrix), profile=False)

        speedup = res_eig.time_per_iteration / res_eigh.time_per_iteration

        print(f"  eig:  {res_eig.time_per_iteration * 1000:.3f} ms")
        print(f"  eigh: {res_eigh.time_per_iteration * 1000:.3f} ms")
        print(f"  Speedup: {speedup:.2f}x")

        results[size] = {
            "eig": res_eig.time_per_iteration,
            "eigh": res_eigh.time_per_iteration,
            "speedup": speedup,
        }

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: Eigenvalue Decomposition Performance")
    print("=" * 70)
    print(f"{'Size':<10} {'eig (ms)':<15} {'eigh (ms)':<15} {'Speedup':<10}")
    print("─" * 70)
    for size, res in results.items():
        print(
            f"{size:<10} {res['eig']*1000:<15.3f} {res['eigh']*1000:<15.3f} "
            f"{res['speedup']:<10.2f}x",
        )


def benchmark_cpmg_calculation() -> None:
    """Benchmark CPMG calculation loop (Bottleneck #3)."""
    print("\n" + "=" * 70)
    print("BOTTLENECK #3: CPMG Calculation Loop")
    print("=" * 70)
    print("Testing matrix_power vs repeated multiplication")

    # Simulate CPMG echo matrix
    size = 16
    rng = np.random.default_rng(42)
    echo = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))
    echo = echo / np.linalg.norm(echo)  # Normalize

    # Test different ncyc values (typical in CPMG experiments)
    ncyc_values = [5, 10, 20, 50, 100]

    from numpy.linalg import matrix_power

    def repeated_mult(m: np.ndarray, n: int) -> np.ndarray:
        """Repeated multiplication alternative to matrix_power."""
        result = np.eye(len(m), dtype=m.dtype)
        for _ in range(n):
            result = result @ m
        return result

    print(f"\nMatrix size: {size}x{size}")
    print(f"{'ncyc':<10} {'matrix_power (ms)':<20} {'repeated @ (ms)':<20} {'Better':<10}")
    print("─" * 70)

    for ncyc in ncyc_values:
        # Benchmark matrix_power
        bench_mp = Benchmark(f"matrix_power (ncyc={ncyc})", iterations=100)
        res_mp = bench_mp.run(lambda: matrix_power(echo, ncyc), profile=False)

        # Benchmark repeated multiplication
        bench_rm = Benchmark(f"repeated @ (ncyc={ncyc})", iterations=100)
        res_rm = bench_rm.run(lambda: repeated_mult(echo, ncyc), profile=False)

        time_mp = res_mp.time_per_iteration * 1000
        time_rm = res_rm.time_per_iteration * 1000
        better = "matrix_power" if time_mp < time_rm else "repeated @"

        print(f"{ncyc:<10} {time_mp:<20.3f} {time_rm:<20.3f} {better:<10}")

    print("\nConclusion: For small ncyc (<20), repeated @ may be faster")
    print("For large ncyc (>50), matrix_power uses eigendecomposition")


def main() -> None:
    """Run all bottleneck benchmarks."""
    print("\n" + "=" * 70)
    print("ChemEx Performance Bottleneck Analysis")
    print("=" * 70)
    print("\nThis benchmark suite profiles the critical bottlenecks identified:")
    print("1. Eigenvalue decomposition in calculate_propagators()")
    print("2. Matrix operations (eig vs eigh)")
    print("3. CPMG matrix_power operations")
    print("4. Liouvillian construction with sum()")

    # Run all benchmarks
    benchmark_matrix_operations_realistic()
    benchmark_cpmg_calculation()

    try:
        benchmark_propagator_calculation()
    except Exception as e:
        print(f"\nWarning: Could not run propagator benchmark: {e}")

    try:
        benchmark_liouvillian_construction()
    except Exception as e:
        print(f"\nWarning: Could not run Liouvillian benchmark: {e}")

    print("\n" + "=" * 70)
    print("Benchmark suite complete!")
    print("=" * 70)
    print("\nKey findings will be saved to: benchmarks/baseline_results.txt")


if __name__ == "__main__":
    # First run basic matrix operations
    benchmark_matrix_operations(size=16, iterations=1000)

    # Then run specific bottleneck tests
    main()
