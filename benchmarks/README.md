# ChemEx Performance Benchmark Suite

This directory contains comprehensive benchmarking tools for analyzing and improving ChemEx performance.

## Overview

The benchmark suite is organized into three levels:

1. **Quick Benchmarks** (`benchmark_framework.py`) - Fast matrix operation tests
2. **Bottleneck Benchmarks** (`benchmark_bottlenecks.py`) - Specific performance bottleneck profiling
3. **End-to-End Benchmarks** (`benchmark_endtoend.py`) - Complete workflow tests with real data

## Quick Start

### Install ChemEx in Development Mode

```bash
cd /home/user/ChemEx
pip install -e .
```

### Run All Benchmarks

```bash
cd benchmarks
python run_benchmarks.py --all
```

### Run Specific Benchmark Suites

```bash
# Quick matrix operations only (~30 seconds)
python run_benchmarks.py --quick

# Bottleneck analysis (~2-3 minutes)
python run_benchmarks.py --bottlenecks

# End-to-end workflows (~5-10 minutes)
python run_benchmarks.py --e2e
```

### Save Results to File

```bash
python run_benchmarks.py --all --save baseline_results.txt
```

## Benchmark Components

### 1. Quick Benchmarks (`benchmark_framework.py`)

Fast tests of core matrix operations:
- Eigenvalue decomposition (`np.linalg.eig` vs `np.linalg.eigh`)
- Matrix inversion
- Matrix power operations
- Repeated matrix multiplication

**Purpose**: Quickly validate that optimization library changes don't break basic operations.

**Runtime**: ~30 seconds

### 2. Bottleneck Benchmarks (`benchmark_bottlenecks.py`)

Profiles the specific bottlenecks identified in the analysis:

1. **Bottleneck #1**: Eigenvalue decomposition in `calculate_propagators()`
   - Located in: `src/chemex/nmr/spectrometer.py:54`
   - Impact: HIGH - Called thousands of times during fitting
   - Test: Realistic Liouvillian sizes (8x8 to 24x24)

2. **Bottleneck #3**: CPMG `matrix_power` operations
   - Located in: `src/chemex/experiments/catalog/cpmg_*.py`
   - Impact: HIGH - Called for each CPMG cycle count
   - Test: Various ncyc values (5, 10, 20, 50, 100)

3. **Bottleneck #5**: Liouvillian construction with `sum()`
   - Located in: `src/chemex/nmr/liouvillian.py:87`
   - Impact: MEDIUM - Uses Python sum() instead of vectorized ops
   - Test: Realistic basis matrices

**Purpose**: Quantify the impact of each bottleneck and validate optimizations.

**Runtime**: ~2-3 minutes

### 3. End-to-End Benchmarks (`benchmark_endtoend.py`)

Complete workflow tests using real example data:

- **CPMG 15N IP Fitting**: Standard CPMG relaxation dispersion analysis
- **CEST 15N TR Fitting**: CEST profile fitting with model-free analysis
- **Simulation**: Fast forward simulation without optimization

**Purpose**: Measure real-world performance on typical user workflows.

**Runtime**: ~5-10 minutes (depends on convergence)

## Benchmark Results

### Baseline Performance (Pre-Optimization)

Run the benchmarks to establish baseline:

```bash
python run_benchmarks.py --all --save baseline_results.txt
```

### Post-Optimization Comparison

After implementing optimizations, run again:

```bash
python run_benchmarks.py --all --save optimized_results.txt
```

Then compare:

```python
from benchmark_framework import compare_benchmarks, BenchmarkResult

# Load and compare results
# (implement result loading as needed)
```

## Understanding the Results

### Matrix Operation Benchmarks

Key metrics to watch:
- **eig vs eigh**: For Hermitian matrices, `eigh` should be 1.5-3x faster
- **matrix_power vs repeated @**: Crossover point depends on matrix size and exponent

Expected speedups from optimizations:
- Using `eigh` for Hermitian Liouvillians: **2-3x** on propagator calculations
- Optimized `matrix_power` for small ncyc: **1.5-2x** on CPMG calculations
- Vectorized Liouvillian construction: **2-5x** on basis operations

### Bottleneck Benchmarks

These show:
- Time spent in each bottleneck function
- Number of calls during typical operations
- cProfile output showing call stacks

Use this to:
1. Validate that optimizations target the right functions
2. Ensure optimizations don't introduce new bottlenecks
3. Measure actual speedup on isolated components

### End-to-End Benchmarks

These show:
- Total workflow execution time
- Top 30 time-consuming functions
- Real-world performance on user workflows

Target speedups:
- **Quick wins** (eigenvalue caching, vectorization): 2-5x on typical fits
- **Parallelization** (grid search, statistics): Near-linear with CPU cores
- **Full optimization** (all improvements): 5-10x potential

## Profiling Individual Components

### Profile a Specific Function

```python
from benchmark_framework import ProfilerContext

with ProfilerContext("My Test") as profiler:
    # Your code here
    result = my_function()

profiler.print_stats(limit=20)
```

### Run Custom Benchmarks

```python
from benchmark_framework import Benchmark

bench = Benchmark("My Custom Test", iterations=100)
result = bench.run(my_function, arg1, arg2, profile=True)
bench.print_results(result)
```

## Performance Targets

Based on the analysis, expected improvements:

| Optimization | Target Speedup | Difficulty |
|--------------|----------------|------------|
| Eigenvalue caching | 2-3x | Low |
| Use eigh for Hermitian | 2-3x | Low |
| Vectorize Liouvillian | 2-5x | Medium |
| Parallelize grid search | 4-8x (on 8 cores) | Low |
| Optimize matrix_power | 1.5-2x | Medium |
| **Combined (Path 1)** | **5-10x** | **Low-Medium** |
| Replace lmfit (Path 2) | +10-30% | High |
| JAX autodiff (Future) | +2-5x | High |

## Continuous Benchmarking

### Add New Benchmarks

To add a new benchmark:

1. Create a function in `benchmark_bottlenecks.py` or `benchmark_endtoend.py`
2. Use the `Benchmark` class or `ProfilerContext` from `benchmark_framework.py`
3. Add to the appropriate main() function
4. Document expected results

### Regression Testing

Run benchmarks before major changes:

```bash
# Before changes
python run_benchmarks.py --all --save before.txt

# Make changes...

# After changes
python run_benchmarks.py --all --save after.txt

# Compare (manual diff for now)
diff before.txt after.txt
```

## Troubleshooting

### Import Errors

If you get import errors:

```bash
# Ensure ChemEx is installed in development mode
pip install -e /home/user/ChemEx
```

### Missing Examples

End-to-end benchmarks require the examples directory:

```bash
# Verify examples exist
ls -la /home/user/ChemEx/examples/Experiments/
```

### Timeout Issues

If benchmarks timeout (default 5 minutes for E2E):

- Edit `benchmark_endtoend.py` and increase `timeout` values
- Or run with `--quick` or `--bottlenecks` only

## References

- **Analysis Documents**: See `/home/user/ChemEx/LMFIT_INVESTIGATION_SUMMARY.md`
- **Bottleneck Details**: See performance analysis in investigation docs
- **Optimization Strategy**: See 5-phase migration plan

## Next Steps

1. **Establish Baseline**: Run `--all --save baseline_results.txt`
2. **Implement Path 1 Optimizations**: Eigenvalue caching, vectorization, parallelization
3. **Measure Improvements**: Run `--all --save optimized_results.txt`
4. **Compare Results**: Validate 2-5x speedup achieved
5. **Path 2 (Optional)**: lmfit replacement if justified by Path 1 results

---

Generated: 2025-11-14
Purpose: Performance optimization planning and validation
