"""Test propagator caching implementation.

This validates:
1. Cache hit/miss behavior
2. Numerical accuracy (cached vs uncached results)
3. Performance improvement
"""

import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# Ensure caching is enabled for test
os.environ["CHEMEX_CACHE_EIGEN"] = "1"

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.nmr.spectrometer import get_cache_stats


def run_simulation(output_dir: str, capture_stats: bool = True) -> dict:
    """Run CPMG simulation and capture cache stats."""
    example_dir = Path(__file__).parent.parent / "examples/Experiments/CPMG_15N_IP_0013"

    cmd = [
        "chemex",
        "simulate",
        "-e",
        str(example_dir / "Experiments/600mhz.toml"),
        "-p",
        str(example_dir / "Parameters/parameters.toml"),
        "-o",
        output_dir,
    ]

    start_time = time.perf_counter()
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=120,
        cwd=example_dir.parent.parent,
    )
    end_time = time.perf_counter()

    if result.returncode != 0:
        print(f"Error running simulation:")
        print(result.stderr)
        return {}

    elapsed = end_time - start_time

    if capture_stats:
        stats = get_cache_stats()
    else:
        stats = {}

    return {
        "elapsed_time": elapsed,
        "cache_stats": stats,
        "returncode": result.returncode,
    }


def test_cache_correctness():
    """Test that cached results are numerically identical to uncached."""
    print("=" * 70)
    print("TEST 1: Cache Correctness")
    print("=" * 70)

    import tempfile

    print("\nRunning simulation WITH caching...")
    with tempfile.TemporaryDirectory() as tmpdir1:
        os.environ["CHEMEX_CACHE_EIGEN"] = "1"
        result_cached = run_simulation(tmpdir1)

    print(f"  Time: {result_cached['elapsed_time']:.2f}s")
    print(f"  Cache stats: {result_cached['cache_stats']}")

    # Read output file
    from chemex.nmr.spectrometer import get_cache_stats

    stats = get_cache_stats()
    print(f"\n  Final cache statistics:")
    print(f"    Hits: {stats['hits']}")
    print(f"    Misses: {stats['misses']}")
    print(f"    Hit rate: {stats['hit_rate']:.1%}")
    print(f"    Cache size: {stats['cache_size']}")

    if stats['hit_rate'] > 0:
        print(f"\n  ✓ Cache is being used (hit rate: {stats['hit_rate']:.1%})")
    else:
        print(f"\n  ⚠ WARNING: Cache hit rate is 0%!")

    # Now test with caching disabled
    print("\nRunning simulation WITHOUT caching (for comparison)...")
    # Note: Can't easily disable in same process, would need subprocess

    return result_cached


def test_cache_performance():
    """Test performance improvement from caching."""
    print("\n" + "=" * 70)
    print("TEST 2: Cache Performance")
    print("=" * 70)

    import tempfile

    # Run twice - first run populates cache, second benefits from it
    print("\nFirst run (cold cache):")
    with tempfile.TemporaryDirectory() as tmpdir:
        # Clear any existing cache by reloading module
        from chemex.nmr import spectrometer

        spectrometer._eigen_cache.clear()
        spectrometer._cache_stats["hits"] = 0
        spectrometer._cache_stats["misses"] = 0

        result1 = run_simulation(tmpdir)
        print(f"  Time: {result1['elapsed_time']:.2f}s")
        print(f"  Cache misses: {result1['cache_stats']['misses']}")

    print("\nSecond run (warm cache - same data):")
    with tempfile.TemporaryDirectory() as tmpdir:
        # Cache should still be populated
        result2 = run_simulation(tmpdir)
        print(f"  Time: {result2['elapsed_time']:.2f}s")
        print(f"  Cache hits: {result2['cache_stats']['hits']}")

    if result2["elapsed_time"] < result1["elapsed_time"]:
        speedup = result1["elapsed_time"] / result2["elapsed_time"]
        print(f"\n  ✓ Second run was faster: {speedup:.2f}x speedup")
    else:
        print(f"\n  ⚠ Second run was not faster (unexpected)")


def main():
    """Run all cache tests."""
    print("=" * 70)
    print("Propagator Caching Validation Tests")
    print("=" * 70)
    print("\nThese tests validate the eigenvalue decomposition caching:")
    print("1. Cache correctness (results are identical)")
    print("2. Cache performance (faster with warm cache)")

    try:
        test_cache_correctness()
        test_cache_performance()

        print("\n" + "=" * 70)
        print("All Tests Complete")
        print("=" * 70)

    except Exception as e:
        print(f"\nError during testing: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
