"""Benchmarking framework for ChemEx performance analysis.

This module provides utilities for profiling ChemEx performance across
different optimization scenarios and identifying bottlenecks.
"""

import cProfile
import pstats
import time
from collections.abc import Callable
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class BenchmarkResult:
    """Results from a benchmark run."""

    name: str
    total_time: float
    iterations: int
    time_per_iteration: float
    memory_peak_mb: float | None = None
    profile_stats: pstats.Stats | None = None
    extra_metrics: dict[str, Any] | None = None


class Benchmark:
    """Base class for performance benchmarks."""

    def __init__(self, name: str, iterations: int = 1) -> None:
        """Initialize benchmark.

        Args:
            name: Descriptive name for this benchmark
            iterations: Number of times to run the benchmark
        """
        self.name = name
        self.iterations = iterations
        self.results: list[BenchmarkResult] = []

    def run(
        self,
        func: Callable,
        *args: Any,
        profile: bool = True,
        **kwargs: Any,
    ) -> BenchmarkResult:
        """Run a benchmark on a function.

        Args:
            func: Function to benchmark
            *args: Positional arguments for the function
            profile: Whether to run cProfile profiling
            **kwargs: Keyword arguments for the function

        Returns:
            BenchmarkResult with timing and profiling data
        """
        times = []

        if profile and self.iterations == 1:
            # Run with profiling
            profiler = cProfile.Profile()
            profiler.enable()
            start = time.perf_counter()
            result = func(*args, **kwargs)
            end = time.perf_counter()
            profiler.disable()

            # Get stats
            stream = StringIO()
            stats = pstats.Stats(profiler, stream=stream)
            stats.strip_dirs()
            stats.sort_stats("cumulative")

            times.append(end - start)
        else:
            # Run without detailed profiling for multiple iterations
            stats = None
            for _ in range(self.iterations):
                start = time.perf_counter()
                result = func(*args, **kwargs)
                end = time.perf_counter()
                times.append(end - start)

        # Calculate statistics
        total_time = sum(times)
        avg_time = total_time / len(times)

        benchmark_result = BenchmarkResult(
            name=self.name,
            total_time=total_time,
            iterations=len(times),
            time_per_iteration=avg_time,
            profile_stats=stats,
        )

        self.results.append(benchmark_result)
        return benchmark_result

    def print_results(self, result: BenchmarkResult | None = None) -> None:
        """Print benchmark results in a readable format.

        Args:
            result: Specific result to print, or None for all results
        """
        results_to_print = [result] if result else self.results

        for res in results_to_print:
            print(f"\n{'=' * 70}")
            print(f"Benchmark: {res.name}")
            print(f"{'=' * 70}")
            print(f"Total time: {res.total_time:.4f} seconds")
            print(f"Iterations: {res.iterations}")
            print(f"Time per iteration: {res.time_per_iteration:.4f} seconds")

            if res.memory_peak_mb:
                print(f"Peak memory: {res.memory_peak_mb:.2f} MB")

            if res.extra_metrics:
                print("\nExtra metrics:")
                for key, value in res.extra_metrics.items():
                    print(f"  {key}: {value}")

            if res.profile_stats:
                print("\nTop 20 time-consuming functions:")
                res.profile_stats.print_stats(20)


class ProfilerContext:
    """Context manager for profiling code blocks."""

    def __init__(self, name: str) -> None:
        """Initialize profiler context.

        Args:
            name: Name for this profiling session
        """
        self.name = name
        self.profiler = cProfile.Profile()
        self.start_time = 0.0
        self.end_time = 0.0

    def __enter__(self) -> "ProfilerContext":
        """Start profiling."""
        self.start_time = time.perf_counter()
        self.profiler.enable()
        return self

    def __exit__(self, *args: Any) -> None:
        """Stop profiling."""
        self.profiler.disable()
        self.end_time = time.perf_counter()

    def print_stats(self, sort_by: str = "cumulative", limit: int = 20) -> None:
        """Print profiling statistics.

        Args:
            sort_by: Sort key for stats (cumulative, time, calls, etc.)
            limit: Number of functions to show
        """
        print(f"\n{'=' * 70}")
        print(f"Profile: {self.name}")
        print(f"Total time: {self.end_time - self.start_time:.4f} seconds")
        print(f"{'=' * 70}")

        stats = pstats.Stats(self.profiler)
        stats.strip_dirs()
        stats.sort_stats(sort_by)
        stats.print_stats(limit)

    def save_stats(self, filepath: Path) -> None:
        """Save profiling statistics to a file.

        Args:
            filepath: Path to save the stats file
        """
        self.profiler.dump_stats(str(filepath))


def compare_benchmarks(
    baseline: BenchmarkResult,
    optimized: BenchmarkResult,
) -> None:
    """Compare two benchmark results and print the comparison.

    Args:
        baseline: Baseline benchmark result
        optimized: Optimized benchmark result
    """
    speedup = baseline.time_per_iteration / optimized.time_per_iteration
    time_saved = baseline.time_per_iteration - optimized.time_per_iteration
    percent_improvement = ((baseline.time_per_iteration - optimized.time_per_iteration)
                          / baseline.time_per_iteration * 100)

    print(f"\n{'=' * 70}")
    print("BENCHMARK COMPARISON")
    print(f"{'=' * 70}")
    print(f"Baseline: {baseline.name}")
    print(f"  Time per iteration: {baseline.time_per_iteration:.4f} seconds")
    print(f"\nOptimized: {optimized.name}")
    print(f"  Time per iteration: {optimized.time_per_iteration:.4f} seconds")
    print(f"\n{'=' * 70}")
    print(f"Speedup: {speedup:.2f}x")
    print(f"Time saved: {time_saved:.4f} seconds per iteration")
    print(f"Improvement: {percent_improvement:.1f}%")
    print(f"{'=' * 70}")


def benchmark_matrix_operations(size: int = 10, iterations: int = 1000) -> None:
    """Benchmark common matrix operations used in ChemEx.

    Args:
        size: Size of the matrices to test
        iterations: Number of iterations to run

    This benchmarks the core operations identified as bottlenecks:
    - Eigenvalue decomposition
    - Matrix inversion
    - Matrix exponential
    - Matrix power
    """
    print(f"\n{'=' * 70}")
    print(f"Matrix Operations Benchmark (size={size}, iterations={iterations})")
    print(f"{'=' * 70}")

    # Create a random complex matrix similar to Liouvillian
    rng = np.random.default_rng(42)
    matrix = rng.standard_normal((size, size)) + 1j * rng.standard_normal((size, size))
    matrix = matrix + matrix.conj().T  # Make Hermitian

    # Benchmark eigenvalue decomposition (current bottleneck)
    bench = Benchmark("Eigenvalue Decomposition (np.linalg.eig)", iterations)
    bench.run(lambda: np.linalg.eig(matrix), profile=False)
    bench.print_results()

    # Benchmark Hermitian eigenvalue decomposition (potential optimization)
    bench = Benchmark("Hermitian Eigenvalue (np.linalg.eigh)", iterations)
    bench.run(lambda: np.linalg.eigh(matrix), profile=False)
    bench.print_results()

    # Benchmark matrix inversion
    bench = Benchmark("Matrix Inversion (np.linalg.inv)", iterations)
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    bench.run(lambda: np.linalg.inv(eigenvectors), profile=False)
    bench.print_results()

    # Benchmark matrix power
    from numpy.linalg import matrix_power

    bench = Benchmark("Matrix Power (matrix_power, n=10)", iterations // 10)
    bench.run(lambda: matrix_power(matrix, 10), profile=False)
    bench.print_results()

    # Benchmark repeated multiplication (alternative to matrix_power)
    def repeated_mult(m: np.ndarray, n: int) -> np.ndarray:
        result = m
        for _ in range(n - 1):
            result = result @ m
        return result

    bench = Benchmark("Matrix Power (repeated @, n=10)", iterations // 10)
    bench.run(lambda: repeated_mult(matrix, 10), profile=False)
    bench.print_results()
