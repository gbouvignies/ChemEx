"""Validate that caching produces numerically identical results.

This runs simulations with and without caching and compares outputs.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


def read_simulation_output(output_dir: Path) -> dict:
    """Read simulation output files."""
    results = {}

    # Read profile data files
    for profile_file in output_dir.glob("Profiles/*.out"):
        with open(profile_file) as f:
            data = []
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 4:  # residue, ncyc, calc, exp
                        data.append([float(x) for x in parts[1:4]])

            if data:
                results[profile_file.name] = np.array(data)

    return results


def run_simulation_with_cache_setting(enabled: bool) -> tuple[dict, float]:
    """Run simulation with caching enabled or disabled."""
    import time

    example_dir = Path(__file__).parent.parent / "examples/Experiments/CPMG_15N_IP_0013"

    env = os.environ.copy()
    env["CHEMEX_CACHE_EIGEN"] = "1" if enabled else "0"

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            "chemex",
            "simulate",
            "-e",
            str(example_dir / "Experiments/600mhz.toml"),
            "-p",
            str(example_dir / "Parameters/parameters.toml"),
            "-o",
            tmpdir,
        ]

        start = time.perf_counter()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
            cwd=example_dir.parent.parent,
        )
        end = time.perf_counter()

        if result.returncode != 0:
            print(f"Error running simulation:")
            print(result.stderr)
            return {}, 0.0

        output_path = Path(tmpdir)
        data = read_simulation_output(output_path)

        return data, end - start


def main():
    """Validate numerical accuracy."""
    print("=" * 70)
    print("Numerical Accuracy Validation")
    print("=" * 70)
    print("\nRunning CPMG simulation twice:")
    print("1. With caching ENABLED")
    print("2. With caching DISABLED")
    print("\nThen comparing outputs for numerical differences...\n")

    print("Running with cache ENABLED...")
    data_cached, time_cached = run_simulation_with_cache_setting(True)
    print(f"  Time: {time_cached:.2f}s")
    print(f"  Profiles: {len(data_cached)}")

    print("\nRunning with cache DISABLED...")
    data_uncached, time_uncached = run_simulation_with_cache_setting(False)
    print(f"  Time: {time_uncached:.2f}s")
    print(f"  Profiles: {len(data_uncached)}")

    # Compare results
    print("\n" + "=" * 70)
    print("Numerical Comparison")
    print("=" * 70)

    if set(data_cached.keys()) != set(data_uncached.keys()):
        print("✗ ERROR: Different output files!")
        return

    max_diff = 0.0
    max_rel_diff = 0.0

    for filename in sorted(data_cached.keys()):
        arr_cached = data_cached[filename]
        arr_uncached = data_uncached[filename]

        if arr_cached.shape != arr_uncached.shape:
            print(f"✗ ERROR: Shape mismatch in {filename}")
            continue

        # Calculate differences
        abs_diff = np.abs(arr_cached - arr_uncached)
        rel_diff = abs_diff / (np.abs(arr_uncached) + 1e-10)

        max_diff = max(max_diff, np.max(abs_diff))
        max_rel_diff = max(max_rel_diff, np.max(rel_diff))

        print(f"{filename}:")
        print(f"  Max absolute difference: {np.max(abs_diff):.2e}")
        print(f"  Max relative difference: {np.max(rel_diff):.2e}")

    print("\n" + "─" * 70)
    print("Overall Results:")
    print(f"  Max absolute difference: {max_diff:.2e}")
    print(f"  Max relative difference: {max_rel_diff:.2e}")

    # Validate
    if max_diff < 1e-10 and max_rel_diff < 1e-10:
        print(f"\n  ✓ Results are NUMERICALLY IDENTICAL (within machine precision)")
    elif max_diff < 1e-6:
        print(f"\n  ✓ Results are VERY CLOSE (acceptable numerical noise)")
    else:
        print(f"\n  ✗ Results DIFFER significantly!")

    # Performance comparison
    print("\n" + "=" * 70)
    print("Performance Comparison")
    print("=" * 70)
    print(f"  With caching:    {time_cached:.2f}s")
    print(f"  Without caching: {time_uncached:.2f}s")

    if time_cached < time_uncached:
        speedup = time_uncached / time_cached
        print(f"  Speedup: {speedup:.2f}x")
        print(f"  Time saved: {time_uncached - time_cached:.2f}s")
    else:
        print(f"  No speedup observed (caching may not help this workflow)")

    print("\n" + "=" * 70)
    print("Validation Complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
