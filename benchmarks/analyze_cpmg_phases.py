"""Analyze CPMG phase cycling pattern for optimization opportunities."""

import sys
from functools import reduce
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.models.loader import register_kinetic_settings
from chemex.experiments.loader import register_experiments

register_kinetic_settings()
register_experiments()


def analyze_phases():
    """Analyze the CPMG phase pattern for optimization opportunities."""
    print("=" * 70)
    print("CPMG PHASE PATTERN ANALYSIS")
    print("=" * 70)

    # The phase pattern from cpmg_15n_ip_0013.py
    cp_phases = np.array(
        [
            [0, 0, 1, 3, 0, 0, 3, 1, 0, 0, 3, 1, 0, 0, 1, 3],
            [1, 3, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 1, 3, 2, 2],
        ],
    )

    print(f"Base phase pattern (repeats every 16 pulses):")
    print(f"  Row 0: {cp_phases[0]}")
    print(f"  Row 1: {cp_phases[1]}")

    # Analyze how phases are used for different ncyc values
    print("\nPhase patterns for different ncyc values:")
    for ncyc in [1, 2, 4, 8, 10, 16, 20, 40]:
        indexes = np.flip(np.arange(2 * int(ncyc)))
        phases = np.take(cp_phases, indexes, mode="wrap", axis=1)

        n_pulses = phases.shape[1]
        n_complete_cycles = n_pulses // 16
        n_remainder = n_pulses % 16

        print(f"\n  ncyc={ncyc:3d}: {n_pulses:3d} pulses, "
              f"{n_complete_cycles} complete 16-pulse cycles + {n_remainder} remainder")

        # Check if the pattern is exactly a multiple of 16
        if n_remainder == 0 and n_complete_cycles > 1:
            print(f"    → Can use matrix_power for {n_complete_cycles} complete cycles")

    # Test optimization: matrix_power vs reduce for large ncyc
    print("\n" + "=" * 70)
    print("OPTIMIZATION TEST: matrix_power for repeating pattern")
    print("=" * 70)

    # Create dummy 6x6 matrices for testing
    np.random.seed(42)
    echo_matrices = np.random.randn(4, 6, 6) * 0.1  # 4 phases

    import time

    ncyc = 40  # 80 pulses = 5 complete 16-pulse cycles
    indexes = np.flip(np.arange(2 * int(ncyc)))
    phases = np.take(cp_phases, indexes, mode="wrap", axis=1)

    # Method 1: Direct reduce
    iterations = 1000
    start = time.perf_counter()
    for _ in range(iterations):
        result1 = reduce(np.matmul, echo_matrices[phases.T])
    t1 = (time.perf_counter() - start) / iterations * 1000
    print(f"\nncyc={ncyc} (80 pulses):")
    print(f"  Method 1 - Direct reduce: {t1:.4f} ms")

    # Method 2: Compute one 16-pulse cycle, then matrix_power
    # First, get the phases for exactly 8 cycles (16 pulses)
    indexes_cycle = np.flip(np.arange(16))
    phases_cycle = np.take(cp_phases, indexes_cycle, mode="wrap", axis=1)

    start = time.perf_counter()
    for _ in range(iterations):
        # Compute the 16-pulse cycle product once
        cycle_product = reduce(np.matmul, echo_matrices[phases_cycle.T])
        # Raise to power for complete cycles
        result2 = np.linalg.matrix_power(cycle_product, 5)
    t2 = (time.perf_counter() - start) / iterations * 1000
    print(f"  Method 2 - matrix_power (5 cycles of 16): {t2:.4f} ms")

    # Check if results match
    print(f"  Results match: {np.allclose(result1, result2)}")
    print(f"  Speedup: {t1/t2:.2f}x")

    # Test for ncyc=20 (40 pulses = 2.5 cycles)
    ncyc = 20
    indexes = np.flip(np.arange(2 * int(ncyc)))
    phases = np.take(cp_phases, indexes, mode="wrap", axis=1)
    n_pulses = phases.shape[1]
    n_complete_cycles = n_pulses // 16
    n_remainder = n_pulses % 16

    start = time.perf_counter()
    for _ in range(iterations):
        result1 = reduce(np.matmul, echo_matrices[phases.T])
    t1 = (time.perf_counter() - start) / iterations * 1000
    print(f"\nncyc={ncyc} ({n_pulses} pulses, {n_complete_cycles} cycles + {n_remainder} remainder):")
    print(f"  Method 1 - Direct reduce: {t1:.4f} ms")

    start = time.perf_counter()
    for _ in range(iterations):
        # Compute the 16-pulse cycle product
        cycle_product = reduce(np.matmul, echo_matrices[phases_cycle.T])
        # Raise to power for complete cycles
        result_cycles = np.linalg.matrix_power(cycle_product, n_complete_cycles)
        # Multiply by remainder
        if n_remainder > 0:
            indexes_rem = np.flip(np.arange(n_remainder))
            phases_rem = np.take(cp_phases, indexes_rem, mode="wrap", axis=1)
            remainder_product = reduce(np.matmul, echo_matrices[phases_rem.T])
            result2 = result_cycles @ remainder_product
        else:
            result2 = result_cycles
    t2 = (time.perf_counter() - start) / iterations * 1000
    print(f"  Method 2 - matrix_power + remainder: {t2:.4f} ms")

    print(f"  Results match: {np.allclose(result1, result2)}")
    print(f"  Speedup: {t1/t2:.2f}x")

    # Analyze overall impact
    print("\n" + "=" * 70)
    print("ESTIMATED IMPACT ON FULL CALCULATION")
    print("=" * 70)

    # Typical ncyc distribution
    ncyc_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40]

    total_time_old = 0
    total_time_new = 0

    for ncyc in ncyc_values:
        indexes = np.flip(np.arange(2 * int(ncyc)))
        phases = np.take(cp_phases, indexes, mode="wrap", axis=1)
        n_pulses = phases.shape[1]
        n_complete_cycles = n_pulses // 16
        n_remainder = n_pulses % 16

        # Time old method
        iterations = 100
        start = time.perf_counter()
        for _ in range(iterations):
            result1 = reduce(np.matmul, echo_matrices[phases.T])
        t_old = (time.perf_counter() - start) / iterations * 1000

        # Time new method
        start = time.perf_counter()
        for _ in range(iterations):
            cycle_product = reduce(np.matmul, echo_matrices[phases_cycle.T])
            if n_complete_cycles > 0:
                result_cycles = np.linalg.matrix_power(cycle_product, n_complete_cycles)
            else:
                result_cycles = np.eye(6)
            if n_remainder > 0:
                indexes_rem = np.flip(np.arange(n_remainder))
                phases_rem = np.take(cp_phases, indexes_rem, mode="wrap", axis=1)
                remainder_product = reduce(np.matmul, echo_matrices[phases_rem.T])
                result2 = result_cycles @ remainder_product
            else:
                result2 = result_cycles
        t_new = (time.perf_counter() - start) / iterations * 1000

        total_time_old += t_old
        total_time_new += t_new

        speedup = t_old / t_new if t_new > 0 else 0
        print(f"ncyc={ncyc:2d}: {t_old:.4f} ms → {t_new:.4f} ms (speedup: {speedup:.2f}x)")

    print(f"\nTotal across all ncyc values:")
    print(f"  Old method: {total_time_old:.4f} ms")
    print(f"  New method: {total_time_new:.4f} ms")
    print(f"  Total speedup: {total_time_old/total_time_new:.2f}x")

    # This is only one iteration; in a fit with 2000 iterations and 48 profiles:
    print(f"\nFor 48 profiles × 2000 iterations:")
    time_saved = (total_time_old - total_time_new) * 48 * 2000 / 1000  # seconds
    print(f"  Potential time saved: {time_saved:.1f} seconds")


if __name__ == "__main__":
    analyze_phases()
