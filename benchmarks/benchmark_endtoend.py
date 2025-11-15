"""End-to-end benchmarks using real ChemEx workflows.

This module benchmarks complete ChemEx workflows using example data:
- CPMG fitting
- CEST fitting
- Grid search
- Statistics calculations
"""

import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# Add benchmarks to path
sys.path.insert(0, str(Path(__file__).parent))

from benchmark_framework import Benchmark, ProfilerContext


class ChemExE2EBenchmark:
    """End-to-end benchmark for ChemEx workflows."""

    def __init__(self, examples_dir: Path) -> None:
        """Initialize with path to examples directory.

        Args:
            examples_dir: Path to ChemEx examples directory
        """
        self.examples_dir = examples_dir
        self.results: dict[str, float] = {}

    def benchmark_cpmg_fit(self) -> dict[str, float]:
        """Benchmark a CPMG fitting workflow.

        Returns:
            Dictionary with timing results
        """
        example = self.examples_dir / "Experiments" / "CPMG_15N_IP_0013"

        if not example.exists():
            print(f"Warning: Example not found at {example}")
            return {}

        print("\n" + "=" * 70)
        print("END-TO-END BENCHMARK: CPMG 15N IP Fitting")
        print("=" * 70)
        print(f"Example: {example}")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            # Build command
            cmd = [
                "chemex",
                "fit",
                "-e",
                str(example / "Experiments" / "*.toml"),
                "-p",
                str(example / "Parameters" / "parameters.toml"),
                "-o",
                str(output_dir),
            ]

            print(f"\nCommand: {' '.join(cmd)}")

            # Run with profiling
            with ProfilerContext("CPMG_15N_IP_Fit") as profiler:
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=300,  # 5 minute timeout
                        cwd=example.parent.parent,
                        shell=False,
                    )

                    if result.returncode != 0:
                        print(f"\nError running chemex:")
                        print(result.stderr)
                        return {}

                except subprocess.TimeoutExpired:
                    print("\nBenchmark timed out after 5 minutes")
                    return {}
                except Exception as e:
                    print(f"\nError running benchmark: {e}")
                    return {}

            profiler.print_stats(limit=30)

            timing = {
                "total_time": profiler.end_time - profiler.start_time,
                "experiment": "CPMG_15N_IP",
            }

            self.results["cpmg_fit"] = timing["total_time"]
            return timing

    def benchmark_cest_fit(self) -> dict[str, float]:
        """Benchmark a CEST fitting workflow.

        Returns:
            Dictionary with timing results
        """
        example = self.examples_dir / "Experiments" / "CEST_15N_TR"

        if not example.exists():
            print(f"Warning: Example not found at {example}")
            return {}

        print("\n" + "=" * 70)
        print("END-TO-END BENCHMARK: CEST 15N TR Fitting")
        print("=" * 70)
        print(f"Example: {example}")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            # Build command
            cmd = [
                "chemex",
                "fit",
                "-e",
                str(example / "Experiments" / "*.toml"),
                "-p",
                str(example / "Parameters" / "parameters.toml"),
                "-m",
                str(example / "Methods" / "method.toml"),
                "-d",
                "2st.mf",
                "-o",
                str(output_dir),
            ]

            print(f"\nCommand: {' '.join(cmd)}")

            # Run with profiling
            with ProfilerContext("CEST_15N_TR_Fit") as profiler:
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=300,  # 5 minute timeout
                        cwd=example.parent.parent,
                        shell=False,
                    )

                    if result.returncode != 0:
                        print(f"\nError running chemex:")
                        print(result.stderr)
                        return {}

                except subprocess.TimeoutExpired:
                    print("\nBenchmark timed out after 5 minutes")
                    return {}
                except Exception as e:
                    print(f"\nError running benchmark: {e}")
                    return {}

            profiler.print_stats(limit=30)

            timing = {
                "total_time": profiler.end_time - profiler.start_time,
                "experiment": "CEST_15N_TR",
            }

            self.results["cest_fit"] = timing["total_time"]
            return timing

    def benchmark_simulate(self) -> dict[str, float]:
        """Benchmark a simulation workflow (no fitting).

        Returns:
            Dictionary with timing results
        """
        example = self.examples_dir / "Experiments" / "CPMG_15N_IP_0013"

        if not example.exists():
            print(f"Warning: Example not found at {example}")
            return {}

        print("\n" + "=" * 70)
        print("END-TO-END BENCHMARK: CPMG Simulation")
        print("=" * 70)
        print(f"Example: {example}")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "output"

            # Build command
            cmd = [
                "chemex",
                "simulate",
                "-e",
                str(example / "Experiments" / "*.toml"),
                "-p",
                str(example / "Parameters" / "parameters.toml"),
                "-o",
                str(output_dir),
            ]

            print(f"\nCommand: {' '.join(cmd)}")

            start = time.perf_counter()
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60,  # 1 minute timeout for simulation
                    cwd=example.parent.parent,
                    shell=False,
                )

                if result.returncode != 0:
                    print(f"\nError running chemex simulate:")
                    print(result.stderr)
                    return {}

            except subprocess.TimeoutExpired:
                print("\nBenchmark timed out after 1 minute")
                return {}
            except Exception as e:
                print(f"\nError running benchmark: {e}")
                return {}

            end = time.perf_counter()

            timing = {
                "total_time": end - start,
                "experiment": "CPMG_Simulation",
            }

            print(f"\nSimulation completed in {timing['total_time']:.2f} seconds")

            self.results["simulate"] = timing["total_time"]
            return timing

    def print_summary(self) -> None:
        """Print summary of all benchmark results."""
        if not self.results:
            print("\nNo benchmark results to display")
            return

        print("\n" + "=" * 70)
        print("END-TO-END BENCHMARK SUMMARY")
        print("=" * 70)

        for name, time_sec in self.results.items():
            print(f"{name:.<40} {time_sec:>10.2f} seconds")

        print("=" * 70)


def main() -> None:
    """Run end-to-end benchmarks."""
    # Find examples directory
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent
    examples_dir = repo_root / "examples"

    if not examples_dir.exists():
        print(f"Error: Examples directory not found at {examples_dir}")
        print("Cannot run end-to-end benchmarks without example data")
        return

    print("=" * 70)
    print("ChemEx End-to-End Performance Benchmarks")
    print("=" * 70)
    print(f"\nExamples directory: {examples_dir}")
    print("\nThese benchmarks run complete ChemEx workflows with real data")
    print("to measure actual performance in realistic scenarios.")

    benchmark = ChemExE2EBenchmark(examples_dir)

    # Run benchmarks
    print("\n" + "─" * 70)
    print("Running benchmarks (this may take several minutes)...")
    print("─" * 70)

    # Simulation first (fastest)
    benchmark.benchmark_simulate()

    # Then fitting workflows
    benchmark.benchmark_cpmg_fit()
    benchmark.benchmark_cest_fit()

    # Print summary
    benchmark.print_summary()

    print("\n" + "=" * 70)
    print("End-to-end benchmarks complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
