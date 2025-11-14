# ChemEx lmfit Investigation - Executive Summary

**Date**: November 14, 2025
**Scope**: Comprehensive audit of lmfit usage for migration planning
**Status**: Complete - 3 detailed reports generated

---

## Quick Facts

- **lmfit Version**: >= 1.3.2 (modern stable API)
- **Files Using lmfit**: 7 files across 2 core modules
- **Critical Files**: 2 (minimizer.py, database.py)
- **Total lmfit Integration Points**: 15+ distinct usage patterns
- **Optimization Methods Supported**: 5 (leastsq, brute, differential_evolution, basinhopping, ampgo)

---

## Key Findings

### 1. Well-Encapsulated Integration
lmfit is used through a clean abstraction layer:
- Parameter database (`src/chemex/parameters/database.py`) handles all Parameter creation/updates
- Minimizer wrapper (`src/chemex/optimize/minimizer.py`) handles all optimization calls
- This makes replacement feasible without touching data containers

### 2. Two Critical Integration Points

**Point 1: Parameter Management**
```
ParamSetting objects ←→ lmfit.Parameters ←→ Optimization
    |                         |
    └── database.py ──────────┘
```

The database module must:
- Convert internal ParamSetting → lmfit.Parameters tuples
- Call `Parameters.add_many()` with constraint expressions
- Evaluate `Parameters.update_constraints()` with usersyms
- Extract `.value` and `.stderr` back to internal database

**Point 2: Optimization Loop**
```
Minimizer Setup → minimize(method) → Result → Statistics
    |                                   |
    └─ minimizer.py ──────────────────┘
```

The minimizer module must:
- Support Minimizer class with callbacks
- Support 5 different optimization methods
- Extract result attributes (nfev, chisqr, redchi, residual)
- Support nested optimization in hierarchical fitting

### 3. Complex Constraint System
ChemEx uses Python expression constraints with custom functions:
```python
expr = "__param_id_1 + 2 * rate_function_nh(__param_id_2, a, b)"
usersyms = {
    "rate_function_nh": callable,
    "rate_function_hn": callable,
    ...
}
```

Replacement must handle:
- Expression parsing and evaluation
- User-defined function calls within expressions
- Constraint resolution during optimization

### 4. Hierarchical Fitting is Critical
Two-level nested optimization pattern:
```
Top Level: lmfit.Minimizer(residuals_hierarchical, params, fcn_args=(tree,))
    ↓
Recursion: For each branch: lmfit.minimize(residuals_hierarchical, params, args=(branch,))
    ↓
Leaf: Return experiments.residuals(params)
```

This requires:
- Support for both Minimizer class and minimize() function
- `fcn_args` parameter passing
- `.residual` attribute on result object

### 5. Standard Error Computation
lmfit automatically computes `.stderr` via Hessian-based methods.
Replacement must provide alternative approach:
- Hessian-based: Requires Jacobian computation
- Monte Carlo: Supported (existing code has bootstrap/MC infrastructure)
- Bootstrap: Supported (existing code mentioned)

---

## Detailed Reports Generated

### 1. **LMFIT_ANALYSIS.MD** (17 KB)
Comprehensive architecture and integration guide:
- All files and their usage patterns (CRITICAL → MODERATE priority)
- Data flow diagrams
- lmfit features required (parameter management, minimizer class, callbacks, etc.)
- Integration points and dependencies
- Migration challenges (ranked by difficulty)
- Feature comparison: scipy.optimize vs JAX approaches

**Key Sections**:
- Files 1-3 describe the 7 files using lmfit
- Sections 2-5 detail architecture, features, and integration patterns
- Section 6 analyzes 7 specific migration challenges
- Section 7 evaluates scipy vs JAX approaches
- Section 9 provides summary table of all files

### 2. **LMFIT_CODE_PATTERNS.MD** (15 KB)
Implementation-level code examples:
- 11 detailed code patterns with line references
- Complete implementations of all major functions
- Parameter object interface documentation
- Minimizer constructor and method signatures
- Result object attributes

**Key Content**:
- Patterns 1-3: Optimization strategies (simple, with callbacks, hierarchical)
- Patterns 4-6: Parameter creation and constraint handling
- Pattern 7: Statistics calculation
- Patterns 8-9: Grid search integration
- Patterns 10-11: Residual calculation and parameter extraction
- Reference tables for all critical APIs

### 3. **LMFIT_METHODS_SUMMARY.MD** (9.3 KB)
Optimization methods and algorithms guide:
- All 5 optimization methods used (with specific kwargs)
- Method selection strategy
- Parameter bounds handling
- Result object field documentation
- Hierarchical fitting architecture (detailed)
- Constraint expression system (detailed)
- Callback mechanism
- Statistics computation post-fit
- Grid search integration pattern
- Critical path summary

**Key Content**:
- Complete method configurations (brute, differential_evolution, basinhopping, ampgo)
- Constraint processing workflow
- Integration points checklist

---

## Migration Complexity Assessment

### Difficulty Levels by Component

**HIGH Difficulty (Must Implement Custom Solutions)**:
1. Constraint expression evaluation with usersyms
   - Risk: scipy.optimize has no native support
   - Effort: Custom expression parser + evaluator
   - Lines of code: 200-300

2. Hierarchical fitting framework
   - Risk: Requires nested minimization pattern
   - Effort: Custom optimization controller
   - Lines of code: 150-200

3. Standard error computation
   - Risk: scipy.optimize doesn't compute uncertainties
   - Effort: Hessian computation or Monte Carlo
   - Lines of code: 100-150

4. Callback uniformity
   - Risk: Different APIs across methods
   - Effort: Adapter pattern for callbacks
   - Lines of code: 100-150

5. Result object normalization
   - Risk: scipy results have different structure
   - Effort: Result wrapper/adapter class
   - Lines of code: 50-100

**MEDIUM Difficulty (Direct Mapping Possible)**:
1. leastsq → scipy.optimize.least_squares()
   - Support: Native callback, bounds support
   - Effort: Direct adaptation
   - Lines of code: 50-80

2. differential_evolution → scipy.optimize.differential_evolution()
   - Support: Native bounds, method kwargs
   - Effort: Direct adaptation
   - Lines of code: 20-30

3. Basic parameter bounds handling
   - Effort: Direct mapping
   - Lines of code: 10-20

**LOW Difficulty (Mostly Unchanged)**:
1. Statistics calculation (already in helper.py)
   - Uses scipy.stats, not lmfit
   - Lines of code: 0 (no change)

2. Grid search pattern
   - Can reuse with new minimizer
   - Lines of code: 0 (no change)

3. Residual calculation
   - Uses Parameters in type annotations only
   - Lines of code: 0 (no change)

### Effort Estimate

```
Component                    Effort    Files Changed
─────────────────────────────────────────────────
Expression system            HIGH      database.py (new module)
Hierarchical fitting         HIGH      minimizer.py
Standard error computation   HIGH      helper.py
Callback adapter             MEDIUM    minimizer.py
Result wrapper               MEDIUM    minimizer.py
scipy.optimize methods       MEDIUM    minimizer.py
Parameter interface          MEDIUM    database.py
Testing & validation         HIGH      tests/
─────────────────────────────────────────────────
TOTAL ESTIMATED EFFORT:      400-600 hours development
                              + 200-300 hours testing
```

---

## Recommended Migration Strategy

### Phase 1: Abstraction Layer (2-3 weeks)
1. Create `src/chemex/optimize/optimization_backend.py`
   - Abstract interface for minimizer
   - Adapter pattern for scipy methods
   - Placeholder for lmfit (keep temporarily)

2. Modify `src/chemex/optimize/minimizer.py`
   - Inject backend dependency
   - No logic changes, just refactoring

**Benefit**: Can test new backend without removing lmfit

### Phase 2: Core scipy.optimize Integration (3-4 weeks)
1. Implement scipy least_squares backend
2. Implement scipy differential_evolution backend
3. Basic result normalization
4. Enable for simple fitting (non-hierarchical)

**Benefit**: Get basic functionality working with scipy

### Phase 3: Complex Features (4-5 weeks)
1. Expression system with sympy/custom parser
2. Constraint resolution framework
3. Callback adapter for all methods
4. Hierarchical fitting adaptation

**Benefit**: Full feature parity with lmfit

### Phase 4: Advanced Features (3-4 weeks)
1. Standard error computation (Hessian or bootstrap)
2. basinhopping and ampgo methods
3. brute force method replacement
4. Grid search optimization

**Benefit**: No loss of functionality

### Phase 5: Testing & Optimization (2-3 weeks)
1. Unit tests for all backends
2. Integration tests with real workflows
3. Performance comparison
4. Deprecate lmfit

**Benefit**: Production readiness

---

## JAX vs scipy.optimize Decision Matrix

### JAX Advantages
- Auto-differentiation for all methods (no manual Jacobian)
- Auto Hessian computation (standard errors for free)
- JIT compilation of constraint expressions
- Better for custom nested optimization
- Handles non-differentiable constraints better

### JAX Disadvantages
- Larger learning curve
- Requires array-based everything
- Compilation overhead for small problems
- Less library ecosystem maturity

### scipy.optimize Advantages
- Minimal changes to existing code
- Well-known library, extensive documentation
- Direct compatibility with numpy workflows
- Mature, battle-tested algorithms

### scipy.optimize Disadvantages
- No automatic differentiation
- Standard errors require custom computation
- Limited constraint expression support
- Callback API inconsistency

**Recommendation**: Start with scipy.optimize for fast path to replacement,
consider JAX for future optimization (separate project).

---

## Key Technical Decisions to Make

1. **Expression Evaluation Method**
   - Option A: sympy (symbolic)
   - Option B: Custom parser (lightweight)
   - Option C: eval() with namespace (simple but unsafe)
   - Recommendation: Custom parser (balance of safety and simplicity)

2. **Standard Error Computation**
   - Option A: Numerical Jacobian + Hessian (fast)
   - Option B: Automatic differentiation (accurate, requires JAX)
   - Option C: Monte Carlo bootstrap (existing infrastructure)
   - Recommendation: Option A initially, Option B if migration includes JAX

3. **Hierarchical Fitting Implementation**
   - Option A: Recursive lmfit.minimize wrapper
   - Option B: Custom optimization controller
   - Option C: Sequential optimization (different algorithm)
   - Recommendation: Option B (cleanest implementation)

4. **Callback Mechanism**
   - Option A: Adapter pattern for each method
   - Option B: Generic callback wrapper
   - Option C: Limit callbacks to leastsq/least_squares
   - Recommendation: Option A (maximum compatibility)

---

## Files to Generate During Migration

```
New files:
  src/chemex/optimize/optimization_backend.py
  src/chemex/optimize/backends/base.py
  src/chemex/optimize/backends/scipy_backend.py
  src/chemex/optimize/backends/jax_backend.py (future)
  src/chemex/optimize/constraints.py
  src/chemex/optimize/result_adapter.py
  src/chemex/optimize/hierarchical_optimizer.py
  tests/optimize/test_backends.py
  tests/optimize/test_constraints.py
  tests/optimize/test_hierarchical.py

Modified files:
  src/chemex/optimize/minimizer.py (significant refactoring)
  src/chemex/parameters/database.py (constraint handling)
  src/chemex/optimize/helper.py (error computation)
  pyproject.toml (remove lmfit, add scipy if not present)
```

---

## Risk Mitigation

1. **Risk**: Breaking existing fits
   **Mitigation**: Keep lmfit available during transition, compare results on test cases

2. **Risk**: Expression evaluation failures
   **Mitigation**: Comprehensive expression parsing tests before deployment

3. **Risk**: Standard error computation differences
   **Mitigation**: Validate on known test cases with manual error computation

4. **Risk**: Performance regression
   **Mitigation**: Benchmark suite comparing lmfit vs scipy vs JAX

5. **Risk**: Callback incompleteness
   **Mitigation**: Gradual method rollout (leastsq first, then others)

---

## References

Generated documents provide:

1. **For Architecture Planning**: LMFIT_ANALYSIS.MD
   - Read Sections 1-5 for overview
   - Read Section 2 for data flow
   - Read Section 6 for challenges

2. **For Implementation**: LMFIT_CODE_PATTERNS.MD
   - Use as reference during coding
   - Copy code patterns and adapt
   - Check method signatures in reference tables

3. **For Method Details**: LMFIT_METHODS_SUMMARY.MD
   - Refer when implementing each method
   - Check integration points checklist
   - Use for testing validation

---

## Next Steps

1. Review all three generated documents
2. Create abstraction layer (Phase 1)
3. Implement scipy.optimize least_squares backend
4. Add comprehensive test suite
5. Evaluate JAX for future iteration
6. Plan deprecation timeline for lmfit

---

## Contact Points

For questions about specific features, refer to:

| Feature | File | Lines | Document |
|---------|------|-------|----------|
| Minimizer class | minimizer.py | 60-104 | CODE_PATTERNS |
| Hierarchical fitting | minimizer.py | 107-167 | ANALYSIS + METHODS |
| Parameter building | database.py | 181-213 | CODE_PATTERNS |
| Constraints | database.py | 200-208 | ANALYSIS |
| Statistics | helper.py | 22-46 | CODE_PATTERNS |
| Grid search | gridding.py | 46-92 | CODE_PATTERNS |
| Methods list | minimizer.py | 65-88 | METHODS_SUMMARY |

---

Generated: 2025-11-14
Author: Automated Investigation System
Status: Ready for Migration Planning
