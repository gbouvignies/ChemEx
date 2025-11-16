"""Profile ChemEx CPMG fit to identify actual bottlenecks.

This instruments the code to measure where time is actually spent.
"""

import cProfile
import pstats
import sys
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def profile_cpmg_fit():
    """Profile a CPMG fit to find actual bottlenecks."""
    from chemex.chemex import main

    print("=" * 70)
    print("PROFILING CHEMEX CPMG FIT")
    print("=" * 70)

    # Run with profiling
    profiler = cProfile.Profile()
    profiler.enable()

    # Temporarily modify sys.argv
    old_argv = sys.argv
    sys.argv = [
        "chemex",
        "fit",
        "-e",
        "examples/Experiments/CPMG_15N_IP_0013/Experiments/600mhz.toml",
        "-e",
        "examples/Experiments/CPMG_15N_IP_0013/Experiments/800mhz.toml",
        "-p",
        "examples/Experiments/CPMG_15N_IP_0013/Parameters/parameters.toml",
        "-o",
        "Output_profile",
    ]

    try:
        main()
    finally:
        sys.argv = old_argv

    profiler.disable()

    # Analyze results
    print("\n" + "=" * 70)
    print("TOP 30 FUNCTIONS BY CUMULATIVE TIME")
    print("=" * 70)

    s = StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.strip_dirs()
    ps.sort_stats("cumulative")
    ps.print_stats(30)
    print(s.getvalue())

    print("\n" + "=" * 70)
    print("TOP 30 FUNCTIONS BY TOTAL TIME (time spent in function itself)")
    print("=" * 70)

    s = StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.strip_dirs()
    ps.sort_stats("tottime")
    ps.print_stats(30)
    print(s.getvalue())

    # Look for specific functions
    print("\n" + "=" * 70)
    print("KEY CHEMEX FUNCTIONS")
    print("=" * 70)

    s = StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.strip_dirs()
    ps.sort_stats("tottime")
    ps.print_stats("liouvillian|spectrometer|propagator|residual|calculate|minimize")
    print(s.getvalue())


if __name__ == "__main__":
    profile_cpmg_fit()
