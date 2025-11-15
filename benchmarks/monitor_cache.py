#!/usr/bin/env python3
"""Monitor cache statistics during ChemEx execution.

Run this as a wrapper around your ChemEx command to see cache behavior.
"""

import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from chemex.nmr.spectrometer import get_cache_stats

def print_cache_stats():
    """Print current cache statistics."""
    stats = get_cache_stats()
    print("\n" + "=" * 60)
    print("CACHE STATISTICS")
    print("=" * 60)
    print(f"Cache enabled: {stats['enabled']}")
    print(f"Total calls: {stats['total']}")
    print(f"Cache hits: {stats['hits']}")
    print(f"Cache misses: {stats['misses']}")
    print(f"Hit rate: {stats['hit_rate']:.1%}")
    print(f"Cache size: {stats['cache_size']}")
    if "fast_path_calls" in stats:
        print(f"Fast path calls (single delay): {stats['fast_path_calls']}")
    print("=" * 60)


if __name__ == "__main__":
    # Import chemex modules to ensure cache is initialized
    import chemex.nmr.spectrometer as spec_module

    # Monkey-patch the calculate_propagators to track calls
    original_func = spec_module.calculate_propagators
    call_count = 0
    single_delay_count = 0
    multi_delay_count = 0
    unique_liouvillians = set()

    import numpy as np

    def tracked_calculate_propagators(liouv, delays, *, dephasing=False):
        global call_count, single_delay_count, multi_delay_count, unique_liouvillians

        call_count += 1
        delays_array = np.asarray(delays).flatten()

        if delays_array.size == 1:
            single_delay_count += 1
        else:
            multi_delay_count += 1
            # Track unique Liouvillians
            unique_liouvillians.add(liouv.tobytes())

        return original_func(liouv, delays, dephasing=dephasing)

    spec_module.calculate_propagators = tracked_calculate_propagators

    # Now run chemex
    print("Starting ChemEx with cache monitoring...")
    print("Cache is ENABLED" if spec_module._CACHE_ENABLED else "Cache is DISABLED")

    # Import and run main
    from chemex.chemex import main

    try:
        main()
    finally:
        # Print statistics
        print_cache_stats()
        print(f"\n--- Call Analysis ---")
        print(f"Total calculate_propagators calls: {call_count}")
        print(f"Single delay calls (fast path): {single_delay_count}")
        print(f"Multi-delay calls (uses eigen cache): {multi_delay_count}")
        print(f"Unique Liouvillians seen: {len(unique_liouvillians)}")
        if multi_delay_count > 0:
            reuse_rate = (multi_delay_count - len(unique_liouvillians)) / multi_delay_count
            print(f"Liouvillian reuse rate: {reuse_rate:.1%}")
