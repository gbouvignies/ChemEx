#!/usr/bin/env python3
"""Main benchmark runner for ChemEx performance analysis.

Usage:
    python run_benchmarks.py [options]

Options:
    --quick         Run only quick matrix operation benchmarks
    --bottlenecks   Run specific bottleneck benchmarks
    --e2e           Run end-to-end workflow benchmarks
    --all           Run all benchmarks (default)
    --save FILE     Save results to FILE
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

# Add to path
sys.path.insert(0, str(Path(__file__).parent))


def run_quick_benchmarks() -> None:
    """Run quick matrix operation benchmarks."""
    from benchmark_framework import benchmark_matrix_operations

    print("\n" + "=" * 70)
    print("QUICK BENCHMARKS: Matrix Operations")
    print("=" * 70)

    benchmark_matrix_operations(size=16, iterations=1000)


def run_bottleneck_benchmarks() -> None:
    """Run specific bottleneck profiling."""
    from benchmark_bottlenecks import main as run_bottlenecks

    print("\n" + "=" * 70)
    print("BOTTLENECK BENCHMARKS")
    print("=" * 70)

    run_bottlenecks()


def run_e2e_benchmarks() -> None:
    """Run end-to-end workflow benchmarks."""
    from benchmark_endtoend import main as run_e2e

    print("\n" + "=" * 70)
    print("END-TO-END WORKFLOW BENCHMARKS")
    print("=" * 70)

    run_e2e()


def main() -> None:
    """Main benchmark runner."""
    parser = argparse.ArgumentParser(
        description="ChemEx Performance Benchmark Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run all benchmarks
    python run_benchmarks.py --all

    # Run only quick matrix benchmarks
    python run_benchmarks.py --quick

    # Run bottleneck analysis
    python run_benchmarks.py --bottlenecks

    # Run end-to-end workflows
    python run_benchmarks.py --e2e

    # Save results to file
    python run_benchmarks.py --all --save baseline_results.txt
        """,
    )

    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run only quick matrix operation benchmarks",
    )
    parser.add_argument(
        "--bottlenecks",
        action="store_true",
        help="Run specific bottleneck benchmarks",
    )
    parser.add_argument(
        "--e2e",
        action="store_true",
        help="Run end-to-end workflow benchmarks",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all benchmarks (default)",
    )
    parser.add_argument(
        "--save",
        type=str,
        metavar="FILE",
        help="Save results to FILE",
    )

    args = parser.parse_args()

    # If no specific benchmark selected, run all
    if not (args.quick or args.bottlenecks or args.e2e):
        args.all = True

    # Print header
    print("=" * 70)
    print("ChemEx Performance Benchmark Suite")
    print("=" * 70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Redirect output if saving
    original_stdout = None
    if args.save:
        import sys

        original_stdout = sys.stdout
        save_file = Path(args.save)
        sys.stdout = open(save_file, "w")
        print(f"ChemEx Benchmark Results - {datetime.now()}")
        print("=" * 70)

    try:
        # Run selected benchmarks
        if args.quick or args.all:
            run_quick_benchmarks()

        if args.bottlenecks or args.all:
            run_bottleneck_benchmarks()

        if args.e2e or args.all:
            run_e2e_benchmarks()

        # Final summary
        print("\n" + "=" * 70)
        print("BENCHMARK SUITE COMPLETE")
        print("=" * 70)
        print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        if args.save:
            print(f"\nResults saved to: {args.save}")

    finally:
        # Restore stdout
        if original_stdout:
            sys.stdout.close()
            sys.stdout = original_stdout
            print(f"\nResults saved to: {args.save}")


if __name__ == "__main__":
    main()
