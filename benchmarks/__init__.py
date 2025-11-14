"""ChemEx Performance Benchmark Suite.

This package provides comprehensive benchmarking tools for ChemEx optimization.

Modules:
    benchmark_framework: Core benchmarking utilities
    benchmark_bottlenecks: Specific bottleneck profiling
    benchmark_endtoend: End-to-end workflow benchmarks
    run_benchmarks: Main benchmark runner

Usage:
    # From command line
    cd benchmarks
    python run_benchmarks.py --all

    # From Python
    from benchmarks.benchmark_framework import Benchmark
    bench = Benchmark("Test", iterations=100)
    result = bench.run(my_function)
"""

__version__ = "0.1.0"
__all__ = [
    "Benchmark",
    "BenchmarkResult",
    "ProfilerContext",
    "compare_benchmarks",
]

from .benchmark_framework import (
    Benchmark,
    BenchmarkResult,
    ProfilerContext,
    compare_benchmarks,
)
