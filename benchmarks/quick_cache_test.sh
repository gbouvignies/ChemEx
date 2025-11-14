#!/bin/bash
# Quick cache performance test
# Usage: ./quick_cache_test.sh path/to/experiments.toml path/to/parameters.toml [path/to/methods.toml]

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "======================================================================"
echo "ChemEx Cache Quick Test"
echo "======================================================================"

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <experiments.toml> <parameters.toml> [methods.toml]"
    echo ""
    echo "Example:"
    echo "  $0 ../examples/Experiments/CPMG_15N_IP_0013/Experiments/600mhz.toml \\"
    echo "     ../examples/Experiments/CPMG_15N_IP_0013/Parameters/parameters.toml"
    exit 1
fi

EXPERIMENTS="$1"
PARAMETERS="$2"
METHODS="${3:-}"

# Verify files exist
if [ ! -f "$EXPERIMENTS" ]; then
    echo -e "${RED}Error: Experiments file not found: $EXPERIMENTS${NC}"
    exit 1
fi

if [ ! -f "$PARAMETERS" ]; then
    echo -e "${RED}Error: Parameters file not found: $PARAMETERS${NC}"
    exit 1
fi

if [ -n "$METHODS" ] && [ ! -f "$METHODS" ]; then
    echo -e "${RED}Error: Methods file not found: $METHODS${NC}"
    exit 1
fi

echo ""
echo "Testing with:"
echo "  Experiments: $EXPERIMENTS"
echo "  Parameters: $PARAMETERS"
if [ -n "$METHODS" ]; then
    echo "  Methods: $METHODS"
fi

# Build chemex command
CMD_BASE="chemex fit -e \"$EXPERIMENTS\" -p \"$PARAMETERS\""
if [ -n "$METHODS" ]; then
    CMD_BASE="$CMD_BASE -m \"$METHODS\""
fi

# Test 1: With cache (3 runs for timing stability)
echo ""
echo "======================================================================"
echo "Running with cache ENABLED (3 runs for stable timing)..."
echo "======================================================================"

export CHEMEX_DEBUG_CACHE=1
TOTAL_CACHED=0

for i in 1 2 3; do
    echo -n "Run $i/3... "
    START=$(date +%s.%N)
    eval "$CMD_BASE -o /tmp/chemex_cache_test_$$ 2> /tmp/cache_stats_$$.txt" > /dev/null 2>&1
    END=$(date +%s.%N)
    TIME=$(echo "$END - $START" | bc)
    TOTAL_CACHED=$(echo "$TOTAL_CACHED + $TIME" | bc)
    echo "${TIME}s"
done

TIME_CACHED=$(echo "scale=2; $TOTAL_CACHED / 3" | bc)

echo ""
echo "Cache statistics (from last run):"
if [ -f /tmp/cache_stats_$$.txt ]; then
    cat /tmp/cache_stats_$$.txt
    rm /tmp/cache_stats_$$.txt
fi

# Test 2: Without cache (3 runs)
echo ""
echo "======================================================================"
echo "Running with cache DISABLED (3 runs for comparison)..."
echo "======================================================================"

export CHEMEX_CACHE_EIGEN=0
unset CHEMEX_DEBUG_CACHE
TOTAL_UNCACHED=0

for i in 1 2 3; do
    echo -n "Run $i/3... "
    START=$(date +%s.%N)
    eval "$CMD_BASE -o /tmp/chemex_uncache_test_$$" > /dev/null 2>&1
    END=$(date +%s.%N)
    TIME=$(echo "$END - $START" | bc)
    TOTAL_UNCACHED=$(echo "$TOTAL_UNCACHED + $TIME" | bc)
    echo "${TIME}s"
done

TIME_UNCACHED=$(echo "scale=2; $TOTAL_UNCACHED / 3" | bc)

# Cleanup
rm -rf /tmp/chemex_cache_test_$$
rm -rf /tmp/chemex_uncache_test_$$

# Calculate speedup
echo ""
echo "======================================================================"
echo "Performance Summary"
echo "======================================================================"
echo "Average time WITH cache:    ${TIME_CACHED}s"
echo "Average time WITHOUT cache: ${TIME_UNCACHED}s"

if (( $(echo "$TIME_CACHED > 0" | bc -l) )); then
    SPEEDUP=$(echo "scale=2; $TIME_UNCACHED / $TIME_CACHED" | bc)
    TIME_SAVED=$(echo "scale=2; $TIME_UNCACHED - $TIME_CACHED" | bc)
    PERCENT=$(echo "scale=1; $TIME_SAVED * 100 / $TIME_UNCACHED" | bc)

    echo ""
    echo "Speedup: ${SPEEDUP}x"
    echo "Time saved: ${TIME_SAVED}s (${PERCENT}%)"
    echo ""

    # Interpretation
    if (( $(echo "$SPEEDUP >= 2.0" | bc -l) )); then
        echo -e "${GREEN}✓ Excellent speedup!${NC} Cache is very effective for this workflow."
        echo "  Recommendation: Keep cache enabled (default)"
    elif (( $(echo "$SPEEDUP >= 1.5" | bc -l) )); then
        echo -e "${GREEN}✓ Good speedup.${NC} Cache provides significant benefit."
        echo "  Recommendation: Keep cache enabled (default)"
    elif (( $(echo "$SPEEDUP >= 1.2" | bc -l) )); then
        echo -e "${YELLOW}⚠ Modest speedup.${NC} Cache provides some benefit."
        echo "  Recommendation: Keep cache enabled, or disable if memory constrained"
    elif (( $(echo "$SPEEDUP >= 1.05" | bc -l) )); then
        echo -e "${YELLOW}⚠ Small speedup.${NC} Cache provides minimal benefit."
        echo "  Recommendation: Can disable if preferred (CHEMEX_CACHE_EIGEN=0)"
    else
        echo -e "${RED}❌ No speedup.${NC} Cache not beneficial for this workflow."
        echo "  Recommendation: Disable cache for this workflow (CHEMEX_CACHE_EIGEN=0)"
        echo "  Possible reasons:"
        echo "    - Workflow uses single-delay calculations (fast path)"
        echo "    - Very small matrices"
        echo "    - Few repeated calculations"
    fi
fi

echo "======================================================================"
echo ""
echo "To disable cache for future runs:"
echo "  export CHEMEX_CACHE_EIGEN=0"
echo ""
echo "To enable debug output:"
echo "  export CHEMEX_DEBUG_CACHE=1"
echo ""
echo "For more details, see: CACHE_TESTING_GUIDE.md"
