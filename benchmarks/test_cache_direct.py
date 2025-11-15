"""Direct test of propagator caching (in-process).

This tests the cache by directly calling calculate_propagators multiple times.
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

# Ensure caching is enabled
os.environ["CHEMEX_CACHE_EIGEN"] = "1"

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.nmr.spectrometer import calculate_propagators, get_cache_stats


def test_cache_with_multiple_delays():
    """Test that cache works when same Liouvillian used with different delays."""
    print("=" * 70)
    print("Direct Cache Test: Multiple Delays, Same Liouvillian")
    print("=" * 70)

    # Create a realistic Liouvillian (16x16 complex matrix)
    rng = np.random.default_rng(42)
    size = 16
    liouv = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))

    # Create different delay arrays
    delays1 = np.array([0.001, 0.002, 0.005, 0.010])  # 4 delays
    delays2 = np.array([0.003, 0.007, 0.015])  # 3 different delays
    delays3 = np.array([0.001, 0.010, 0.020])  # 3 delays, some overlap

    print(f"\nLiouvillian shape: {liouv.shape}")
    print(f"Test 1 delays: {delays1}")
    print(f"Test 2 delays: {delays2}")
    print(f"Test 3 delays: {delays3}")

    # First call - should miss cache
    print("\nCall 1 (should be cache MISS):")
    stats_before = get_cache_stats()
    result1 = calculate_propagators(liouv, delays1)
    stats_after = get_cache_stats()

    print(f"  Result shape: {result1.shape}")
    print(f"  Cache misses: {stats_after['misses'] - stats_before['misses']}")
    print(f"  Cache hits: {stats_after['hits'] - stats_before['hits']}")

    # Second call with SAME Liouvillian but DIFFERENT delays - should hit cache
    print("\nCall 2 (same Liouvillian, different delays - should HIT cache):")
    stats_before = get_cache_stats()
    result2 = calculate_propagators(liouv, delays2)
    stats_after = get_cache_stats()

    print(f"  Result shape: {result2.shape}")
    print(f"  Cache misses: {stats_after['misses'] - stats_before['misses']}")
    print(f"  Cache hits: {stats_after['hits'] - stats_before['hits']}")

    # Third call - should also hit cache
    print("\nCall 3 (same Liouvillian, different delays - should HIT cache):")
    stats_before = get_cache_stats()
    result3 = calculate_propagators(liouv, delays3)
    stats_after = get_cache_stats()

    print(f"  Result shape: {result3.shape}")
    print(f"  Cache misses: {stats_after['misses'] - stats_before['misses']}")
    print(f"  Cache hits: {stats_after['hits'] - stats_before['hits']}")

    # Print overall stats
    final_stats = get_cache_stats()
    print(f"\n{'─' * 70}")
    print("Final Cache Statistics:")
    print(f"  Total calls to eigen decomposition: {final_stats['total']}")
    print(f"  Cache hits: {final_stats['hits']}")
    print(f"  Cache misses: {final_stats['misses']}")
    print(f"  Hit rate: {final_stats['hit_rate']:.1%}")
    print(f"  Cache size: {final_stats['cache_size']}")

    if final_stats["hit_rate"] > 0.5:
        print(f"\n  ✓ Cache is working correctly!")
    else:
        print(f"\n  ✗ Cache hit rate too low")


def test_cache_performance():
    """Benchmark cached vs uncached performance."""
    print("\n" + "=" * 70)
    print("Performance Test: Cached vs Uncached")
    print("=" * 70)

    rng = np.random.default_rng(42)
    size = 16
    liouv = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))

    # Test with many different delay arrays (simulating CPMG experiment)
    delay_arrays = [
        np.array([0.001, 0.002, 0.005, 0.010])
        for _ in range(50)  # 50 calls with same delays
    ]

    # Clear cache
    from chemex.nmr.spectrometer import _eigen_cache

    _eigen_cache.clear()

    # Measure time for 50 calls (first will miss, rest will hit)
    start = time.perf_counter()
    for delays in delay_arrays:
        _ = calculate_propagators(liouv, delays)
    end = time.perf_counter()

    cached_time = end - start
    stats = get_cache_stats()

    print(f"\n50 calls with caching:")
    print(f"  Total time: {cached_time:.4f}s")
    print(f"  Time per call: {cached_time / 50 * 1000:.2f}ms")
    print(f"  Cache hits: {stats['hits']}")
    print(f"  Cache misses: {stats['misses']}")
    print(f"  Hit rate: {stats['hit_rate']:.1%}")

    # Now disable cache and measure again
    os.environ["CHEMEX_CACHE_EIGEN"] = "0"
    # Need to reload module to pick up env var
    import importlib

    import chemex.nmr.spectrometer

    importlib.reload(chemex.nmr.spectrometer)
    from chemex.nmr.spectrometer import calculate_propagators as calc_uncached

    start = time.perf_counter()
    for delays in delay_arrays:
        _ = calc_uncached(liouv, delays)
    end = time.perf_counter()

    uncached_time = end - start

    print(f"\n50 calls WITHOUT caching:")
    print(f"  Total time: {uncached_time:.4f}s")
    print(f"  Time per call: {uncached_time / 50 * 1000:.2f}ms")

    speedup = uncached_time / cached_time
    print(f"\n{'─' * 70}")
    print(f"Speedup from caching: {speedup:.2f}x")
    print(f"Time saved: {uncached_time - cached_time:.4f}s ({(uncached_time - cached_time) / uncached_time * 100:.1f}%)")

    if speedup > 1.5:
        print(f"  ✓ Significant speedup achieved!")
    else:
        print(f"  ⚠ Speedup less than expected")


def main():
    """Run all direct cache tests."""
    test_cache_with_multiple_delays()
    test_cache_performance()

    print("\n" + "=" * 70)
    print("Direct Cache Tests Complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
