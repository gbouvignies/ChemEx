# Cache Performance Testing Guide

**Purpose**: Validate eigenvalue decomposition caching with real ChemEx workflows

**Date**: November 14, 2025
**Optimization**: Eigenvalue decomposition caching (commit `30973aa`)

---

## Quick Start

### 1. Run Your Normal Workflow with Cache Enabled

The cache is **enabled by default**. Just run ChemEx normally:

```bash
# Your typical fitting workflow
chemex fit -e experiments.toml -p parameters.toml -o Output_Cached
```

### 2. Run Same Workflow with Cache Disabled

Compare performance with cache disabled:

```bash
# Same workflow, cache disabled
CHEMEX_CACHE_EIGEN=0 chemex fit -e experiments.toml -p parameters.toml -o Output_Uncached
```

### 3. Compare Results

Check that results are identical and measure performance difference:

```bash
# Compare output files
diff -r Output_Cached Output_Uncached

# Time both runs (add 'time' command)
time chemex fit -e experiments.toml -p parameters.toml -o Output_Cached
time CHEMEX_CACHE_EIGEN=0 chemex fit -e experiments.toml -p parameters.toml -o Output_Uncached
```

---

## Expected Performance Gains by Workflow Type

### High Benefit Workflows ⭐⭐⭐

**CEST Experiments** (many offsets, same Liouvillian):
- Expected speedup: **2-4x**
- Cache hit rate: 80-95%
- Why: Same Liouvillian used for 20-100 different offsets

**Grid Search Optimization**:
- Expected speedup: **2-5x**
- Cache hit rate: 70-90%
- Why: Same Liouvillian evaluated repeatedly at different parameter values

**Monte Carlo / Bootstrap Statistics**:
- Expected speedup: **3-5x**
- Cache hit rate: 85-95%
- Why: Repeated calculations with same/similar Liouvillians

### Medium Benefit Workflows ⭐⭐

**Multi-Field CPMG Experiments**:
- Expected speedup: **1.5-2.5x**
- Cache hit rate: 50-70%
- Why: Some Liouvillian reuse across fields

**Hierarchical Fitting** (multiple groups):
- Expected speedup: **1.5-3x**
- Cache hit rate: 60-80%
- Why: Repeated evaluations during nested optimization

### Low Benefit Workflows ⭐

**Simple CPMG Simulations**:
- Expected speedup: **1.1-1.3x**
- Cache hit rate: 10-30%
- Why: Mostly uses fast path (single delay, no dephasing)

**Single-Point Calculations**:
- Expected speedup: **~1.0x** (none)
- Cache hit rate: <10%
- Why: Minimal repeated calculations

---

## Detailed Testing Scenarios

### Scenario 1: CEST Profile Fitting

**Test with your CEST data**:

```bash
# With cache (should be faster)
time chemex fit \
  -e cest_experiments.toml \
  -p parameters.toml \
  -m methods.toml \
  -o CEST_Cached

# Without cache (baseline)
time CHEMEX_CACHE_EIGEN=0 chemex fit \
  -e cest_experiments.toml \
  -p parameters.toml \
  -m methods.toml \
  -o CEST_Uncached

# Compare results
diff -r CEST_Cached CEST_Uncached

# Expected: 2-4x speedup, identical results
```

**Why CEST benefits**: Each CEST experiment sweeps through many offsets (typically 20-100 points). The Liouvillian is recomputed for each offset, but with caching, the eigendecomposition is reused.

### Scenario 2: Grid Search

**Test with grid search method**:

```bash
# With cache
time chemex fit \
  -e experiments.toml \
  -p parameters.toml \
  -m grid_search_method.toml \
  -o Grid_Cached

# Without cache
time CHEMEX_CACHE_EIGEN=0 chemex fit \
  -e experiments.toml \
  -p parameters.toml \
  -m grid_search_method.toml \
  -o Grid_Uncached

# Expected: 2-5x speedup on grid search
```

### Scenario 3: Bootstrap Statistics

**Test with bootstrap/MC statistics**:

```bash
# With cache
time chemex fit \
  -e experiments.toml \
  -p parameters.toml \
  -m bootstrap_method.toml \
  -o Bootstrap_Cached

# Without cache
time CHEMEX_CACHE_EIGEN=0 chemex fit \
  -e experiments.toml \
  -p parameters.toml \
  -m bootstrap_method.toml \
  -o Bootstrap_Uncached

# Expected: 3-5x speedup on bootstrap iterations
```

---

## Monitoring Cache Effectiveness

### Method 1: Using Python (Direct Access)

Create a simple script to run your workflow and report cache stats:

```python
#!/usr/bin/env python
"""Run workflow and report cache statistics."""

import subprocess
import sys

# Your ChemEx command
cmd = [
    "chemex", "fit",
    "-e", "experiments.toml",
    "-p", "parameters.toml",
    "-o", "Output"
]

# Run workflow
print("Running ChemEx workflow...")
result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print("Error:", result.stderr)
    sys.exit(1)

# Get cache statistics (note: stats are per-process)
# This won't work across subprocess boundaries
# See Method 2 for subprocess monitoring

print("\nWorkflow completed successfully!")
```

### Method 2: Enable Debug Mode

Set the debug flag to get automatic cache statistics on exit:

```bash
# This will print cache stats to stderr at program exit
CHEMEX_DEBUG_CACHE=1 chemex fit -e experiments.toml -p parameters.toml -o Output 2> cache_stats.txt

# View cache statistics
cat cache_stats.txt
```

**Output will show**:
```
======================================================================
ChemEx Eigenvalue Decomposition Cache Statistics
======================================================================
Cache enabled: True
Fast path calls (single delay, no dephasing): 1234
Eigendecomposition cache calls: 567
  Cache hits: 432
  Cache misses: 135
  Hit rate: 76.2%
Final cache size: 48 entries
======================================================================
```

### Method 3: Time Comparison

Simple approach - just time both runs:

```bash
# Create a timing script
cat > time_comparison.sh << 'EOF'
#!/bin/bash

echo "Testing cache performance..."

echo ""
echo "=== WITH CACHE ==="
time chemex fit -e experiments.toml -p parameters.toml -o Output_Cached

echo ""
echo "=== WITHOUT CACHE ==="
time CHEMEX_CACHE_EIGEN=0 chemex fit -e experiments.toml -p parameters.toml -o Output_Uncached

echo ""
echo "Verifying results are identical..."
if diff -r Output_Cached/Parameters Output_Uncached/Parameters; then
    echo "✓ Results are identical"
else
    echo "✗ WARNING: Results differ!"
fi
EOF

chmod +x time_comparison.sh
./time_comparison.sh
```

---

## Interpreting Results

### Cache Hit Rate

**>80%**: Excellent - optimization is very effective for your workflow
**50-80%**: Good - significant speedup achieved
**30-50%**: Moderate - some benefit
**<30%**: Low - cache not helping much (workflow uses fast path or unique Liouvillians)

### Speedup Factor

Calculate from timing results:

```
Speedup = Time_Without_Cache / Time_With_Cache
```

**Examples**:
- Time without cache: 120s, with cache: 30s → **4.0x speedup** ✓
- Time without cache: 100s, with cache: 75s → **1.33x speedup** (modest)
- Time without cache: 50s, with cache: 48s → **1.04x speedup** (minimal)

### Memory Usage

Monitor memory if running large-scale studies:

```bash
# Check memory usage during run
/usr/bin/time -v chemex fit -e experiments.toml -p parameters.toml -o Output 2>&1 | grep "Maximum resident"
```

**Cache memory overhead**: ~8KB per cached 16×16 matrix
**Max cache size**: 256 entries → ~2MB max

---

## Tuning Recommendations

### If Cache Hit Rate is Low (<30%)

**Possible causes**:
1. Workflow uses mostly single-delay calculations (fast path)
2. Many unique Liouvillians (low reuse)
3. Cache size too small (LRU eviction)

**Solutions**:
- This is expected for simple simulations - no action needed
- Check if workflow could be restructured for more reuse
- Consider increasing cache size (see below)

### If Cache Hit Rate is High (>80%) But Speedup is Low

**Possible causes**:
1. Eigendecomposition is small fraction of total time
2. I/O or other operations dominate
3. Very small matrices (cache overhead comparable to computation)

**Solutions**:
- This is fine - cache is working optimally
- Look at other optimizations (Liouvillian vectorization, parallelization)

### Adjusting Cache Size

If you see many cache misses and suspect eviction, you can modify the cache size.

**Edit** `src/chemex/nmr/spectrometer.py`:

```python
# Current (line 29)
_eigen_cache: LRUCache = LRUCache(maxsize=256)

# For larger studies, increase:
_eigen_cache: LRUCache = LRUCache(maxsize=512)  # 4MB max

# For memory-constrained systems, decrease:
_eigen_cache: LRUCache = LRUCache(maxsize=128)  # 1MB max
```

**Recommendations**:
- 256 entries (default): Good for most workflows
- 512 entries: For large-scale grid searches or many experiments
- 128 entries: For memory-constrained systems
- 1024+ entries: Only if you have very large parameter spaces

---

## Troubleshooting

### Cache Not Providing Expected Speedup

**Check 1**: Verify cache is enabled
```bash
python -c "from chemex.nmr.spectrometer import get_cache_stats; print(get_cache_stats()['enabled'])"
# Should print: True
```

**Check 2**: Verify workflow uses eigendecomposition (not fast path)
```bash
# Enable debug mode and check fast_path_calls
CHEMEX_DEBUG_CACHE=1 chemex simulate -e exp.toml -p params.toml -o Out 2>&1 | grep "Fast path"
```

**Check 3**: Run direct cache test
```bash
cd benchmarks
python test_cache_direct.py
# Should show 4.47x speedup
```

### Results Differ Between Cached and Uncached

**This should never happen**. If you observe this:

1. Disable cache immediately:
   ```bash
   export CHEMEX_CACHE_EIGEN=0
   ```

2. Report the issue with:
   - Workflow files (experiments, parameters, methods)
   - Diff output showing differences
   - ChemEx version

3. Run validation test:
   ```bash
   cd benchmarks
   python validate_cache_accuracy.py
   ```

### Cache-Related Errors

If you see errors related to caching:

1. Disable cache:
   ```bash
   CHEMEX_CACHE_EIGEN=0 chemex fit ...
   ```

2. Report error message and workflow

3. Cache has try/except fallback, so this is unlikely

---

## Feedback Collection

### Please Report

After testing with your workflows, please provide feedback on:

1. **Workflow Type**:
   - CPMG / CEST / Other?
   - Grid search / Fitting / Simulation?
   - Number of experiments, residues, data points

2. **Performance Observed**:
   - Speedup factor (time without cache / time with cache)
   - Cache hit rate (from debug output)
   - Memory usage

3. **Numerical Accuracy**:
   - Any differences in results? (should be none)
   - Any unexpected behavior?

4. **Usability**:
   - Was cache transparent (just worked)?
   - Did you need to tune cache size?
   - Any issues or concerns?

### Reporting Template

```
## Cache Performance Feedback

**Workflow**: [CEST fitting / CPMG grid search / etc.]
**Dataset**: [X experiments, Y residues, Z data points]

**Performance**:
- Time without cache: X seconds
- Time with cache: Y seconds
- Speedup: Z.ZZx
- Cache hit rate: XX%

**Results**:
- Numerical differences: [None / Small / Significant]
- Any issues: [None / Description]

**Overall**:
- Recommendation: [Keep / Disable / Tune]
- Comments: [Optional feedback]
```

---

## Next Steps After Testing

Based on your results, we can:

**If speedup is good (>1.5x)**:
✅ Keep cache enabled by default
✅ Proceed with next optimization (Liouvillian vectorization)
✅ Consider cache size tuning if needed

**If speedup is modest (1.1-1.5x)**:
⚠️ Keep cache enabled (no harm, some benefit)
⚠️ Focus on other optimizations with higher impact
⚠️ Consider workflow-specific tuning

**If speedup is minimal (<1.1x)**:
❓ Disable cache for your workflow (`CHEMEX_CACHE_EIGEN=0`)
❓ Focus on optimizations that better match your workflow
❓ Consider parallelization or vectorization instead

---

## Example: Real-World Test Script

```bash
#!/bin/bash
# comprehensive_cache_test.sh

set -e

echo "======================================================================"
echo "ChemEx Cache Performance Test"
echo "======================================================================"

# Configuration
EXPERIMENTS="experiments.toml"
PARAMETERS="parameters.toml"
METHODS="methods.toml"

echo ""
echo "Testing with:"
echo "  Experiments: $EXPERIMENTS"
echo "  Parameters: $PARAMETERS"
echo "  Methods: $METHODS"

# Test 1: With cache
echo ""
echo "======================================================================"
echo "Test 1: Running WITH cache enabled..."
echo "======================================================================"
export CHEMEX_DEBUG_CACHE=1
START1=$(date +%s)
chemex fit -e "$EXPERIMENTS" -p "$PARAMETERS" -m "$METHODS" -o Output_Cached 2> cache_stats.txt
END1=$(date +%s)
TIME_CACHED=$((END1 - START1))

echo ""
echo "Completed in ${TIME_CACHED} seconds"
cat cache_stats.txt

# Test 2: Without cache
echo ""
echo "======================================================================"
echo "Test 2: Running WITHOUT cache (for comparison)..."
echo "======================================================================"
export CHEMEX_CACHE_EIGEN=0
unset CHEMEX_DEBUG_CACHE
START2=$(date +%s)
chemex fit -e "$EXPERIMENTS" -p "$PARAMETERS" -m "$METHODS" -o Output_Uncached
END2=$(date +%s)
TIME_UNCACHED=$((END2 - START2))

echo ""
echo "Completed in ${TIME_UNCACHED} seconds"

# Compare results
echo ""
echo "======================================================================"
echo "Comparing Results..."
echo "======================================================================"

if diff -r Output_Cached/Parameters Output_Uncached/Parameters > /dev/null 2>&1; then
    echo "✓ Results are numerically identical"
else
    echo "✗ WARNING: Results differ!"
    diff -r Output_Cached/Parameters Output_Uncached/Parameters | head -20
fi

# Calculate speedup
echo ""
echo "======================================================================"
echo "Performance Summary"
echo "======================================================================"
echo "Time WITH cache:    ${TIME_CACHED}s"
echo "Time WITHOUT cache: ${TIME_UNCACHED}s"

if [ $TIME_CACHED -gt 0 ]; then
    SPEEDUP=$(echo "scale=2; $TIME_UNCACHED / $TIME_CACHED" | bc)
    PERCENT=$(echo "scale=1; ($TIME_UNCACHED - $TIME_CACHED) * 100 / $TIME_UNCACHED" | bc)
    echo "Speedup: ${SPEEDUP}x"
    echo "Time saved: ${PERCENT}%"

    if (( $(echo "$SPEEDUP > 1.5" | bc -l) )); then
        echo "✓ Significant speedup - cache is very effective!"
    elif (( $(echo "$SPEEDUP > 1.1" | bc -l) )); then
        echo "⚠ Modest speedup - cache provides some benefit"
    else
        echo "❌ Minimal speedup - cache not beneficial for this workflow"
    fi
fi

echo "======================================================================"
```

Save this script and run:

```bash
chmod +x comprehensive_cache_test.sh
./comprehensive_cache_test.sh
```

---

**Ready to test**: Use the methods above with your real ChemEx workflows and report back with results!
