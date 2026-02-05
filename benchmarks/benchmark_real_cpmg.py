#!/usr/bin/env python3
"""Benchmark using real CPMG example from examples/Experiments/CPMG_15N_IP_0013.

This profiles the actual calculation path, not synthetic data.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Must import model registration first
from chemex.models.loader import register_kinetic_settings
from chemex.experiments.loader import register_experiments

register_kinetic_settings()
register_experiments()


def load_real_experiment():
    """Load the real CPMG experiment from examples."""
    from chemex.experiments.builder import build_experiments
    from chemex.configuration.methods import Selection

    example_dir = Path(__file__).parent.parent / "examples" / "Experiments" / "CPMG_15N_IP_0013"

    # Default selection includes everything
    selection = Selection(include="*", exclude=None)

    experiments = build_experiments([
        example_dir / "Experiments" / "600mhz.toml",
        example_dir / "Experiments" / "800mhz.toml",
    ], selection)

    return experiments


def benchmark_residual_calculation():
    """Benchmark the residual calculation (what happens during fitting)."""
    from chemex.parameters import database

    print("=" * 70)
    print("REAL CPMG BENCHMARK")
    print("=" * 70)

    experiments = load_real_experiment()

    # Get all profiles
    all_profiles = []
    for exp in experiments:
        all_profiles.extend(exp.profiles)

    print(f"Total profiles loaded: {len(all_profiles)}")
    print(f"Total data points: {sum(len(p.data.exp) for p in all_profiles)}")

    # Get parameters (built from experiment param_ids)
    params = database.build_lmfit_params(experiments.param_ids)
    print(f"Total parameters: {len(params)}")

    # Warm up
    print("\nWarming up...")
    for _ in range(3):
        experiments.residuals(params)

    # Benchmark single residual calculation (what minimizer calls)
    print("\n--- Benchmark: Single Residual Calculation ---")
    iterations = 100
    start = time.perf_counter()
    for _ in range(iterations):
        experiments.residuals(params)
    elapsed = time.perf_counter() - start

    time_per_call = elapsed / iterations * 1000  # ms
    print(f"Time per residual call: {time_per_call:.4f} ms")
    print(f"Calls per second: {iterations/elapsed:.0f}")

    # Typical fit has ~2000 iterations
    expected_fit_time = time_per_call * 2000 / 1000  # seconds
    print(f"\nExpected fit time (2000 iterations): {expected_fit_time:.1f} seconds")

    # Break down the time
    print("\n--- Breakdown: Profile Components ---")

    # Time spectrometer update
    iterations = 1000
    profile = all_profiles[0]
    start = time.perf_counter()
    for _ in range(iterations):
        profile.update_spectrometer(params)
    elapsed = time.perf_counter() - start
    print(f"update_spectrometer: {elapsed/iterations*1000:.4f} ms")

    # Time calculation
    start = time.perf_counter()
    for _ in range(iterations):
        profile.pulse_sequence.calculate(profile.spectrometer, profile.data)
    elapsed = time.perf_counter() - start
    print(f"pulse_sequence.calculate: {elapsed/iterations*1000:.4f} ms")

    # Time _build_base_liouvillian specifically
    liouvillian = profile.spectrometer.liouvillian
    par_values = profile._get_parameter_values(params)

    start = time.perf_counter()
    for _ in range(iterations):
        liouvillian.par_values = par_values
        liouvillian._build_base_liouvillian()
    elapsed = time.perf_counter() - start
    print(f"_build_base_liouvillian: {elapsed/iterations*1000:.4f} ms")

    # Compare old vs new implementation
    print("\n--- Compare Old vs New _build_base_liouvillian ---")

    # Old implementation (inline)
    import numpy as np
    def old_build(liouvillian):
        return sum(
            (
                liouvillian.basis.matrices[name] * liouvillian.par_values.get(name, 0.0)
                for name in liouvillian.basis.matrices
            ),
            start=np.array(0.0),
        )

    iterations = 10000
    liouvillian.par_values = par_values

    start = time.perf_counter()
    for _ in range(iterations):
        _ = old_build(liouvillian)
    old_time = time.perf_counter() - start
    print(f"Old (sum generator): {old_time/iterations*1000:.4f} ms")

    start = time.perf_counter()
    for _ in range(iterations):
        liouvillian._build_base_liouvillian()
    new_time = time.perf_counter() - start
    print(f"New (einsum):        {new_time/iterations*1000:.4f} ms")

    speedup = old_time / new_time
    print(f"Speedup: {speedup:.2f}x")

    # Calculate total impact
    num_profiles = len(all_profiles)
    num_iterations = 2000
    total_calls = num_profiles * num_iterations

    time_saved_per_call = (old_time - new_time) / iterations  # seconds
    total_time_saved = time_saved_per_call * total_calls

    print(f"\n--- Impact on Full Fit ---")
    print(f"Profiles: {num_profiles}")
    print(f"Iterations: {num_iterations}")
    print(f"Total _build_base_liouvillian calls: {total_calls}")
    print(f"Time saved per call: {time_saved_per_call*1000:.4f} ms")
    print(f"Total time saved: {total_time_saved:.2f} seconds")
    print(f"As percentage of 60s fit: {total_time_saved/60*100:.1f}%")


if __name__ == "__main__":
    benchmark_residual_calculation()
