# ChemEx lmfit Optimization Methods Summary

## All Optimization Methods Used in ChemEx

### 1. **"leastsq"** - Levenberg-Marquardt Algorithm
**File**: `src/chemex/optimize/minimizer.py`
**Usage**: 
- Basic fitting method (default)
- Hierarchical fitting sub-optimization
- Best for small to medium problems with fast convergence needed

**Called as**:
```python
result = minimizer.minimize(method="leastsq")
```

**No special kwargs**: Uses scipy's default L-M implementation via lmfit

---

### 2. **"brute"** - Brute Force Grid Search
**File**: `src/chemex/optimize/minimizer.py`
**Usage**: 
- Global optimization for exploration
- Configured with `{"keep": "all"}` to keep all grid results

**Called as**:
```python
kws = {"brute": {"keep": "all"}}
result = minimizer.minimize(method="brute", keep="all")
```

**Purpose**: Exhaustive search of parameter space (uses `brute_step` from ParamSetting)

---

### 3. **"differential_evolution"** - Global Differential Evolution
**File**: `src/chemex/optimize/minimizer.py`
**Usage**: 
- Global optimization algorithm
- Good for non-convex, multi-modal problems

**Called as**:
```python
kws = {"differential_evolution": {"disp": True}}
result = minimizer.minimize(method="differential_evolution", disp=True)
```

**Features**:
- Display progress (disp=True)
- Handles parameter bounds naturally

---

### 4. **"basinhopping"** - Basin-Hopping Global Optimizer
**File**: `src/chemex/optimize/minimizer.py`
**Usage**: 
- Global optimization using local jumps
- Sophisticated global search strategy

**Called as**:
```python
kws = {
    "basinhopping": {
        "disp": True,
        "niter_success": 10,
        "minimizer_kwargs": {"method": "L-BFGS-B"},
    }
}
result = minimizer.minimize(method="basinhopping", **kws["basinhopping"])
```

**Parameters**:
- `disp`: Display optimization progress
- `niter_success`: Number of iterations without improvement before stopping
- `minimizer_kwargs`: Passes `{"method": "L-BFGS-B"}` as the local minimizer

**Local minimizer**: L-BFGS-B (limited memory BFGS with bounds)

---

### 5. **"ampgo"** - Autonomous Metropolis-Hopping Global Optimizer
**File**: `src/chemex/optimize/minimizer.py`
**Usage**: 
- Global optimization algorithm
- Hybrid between local and global approaches

**Called as**:
```python
kws = {"ampgo": {"disp": True}}
result = minimizer.minimize(method="ampgo", disp=True)
```

**Features**:
- Display progress (disp=True)
- Autonomous parameter adjustment

---

## Method Selection Strategy (from code analysis)

### Simple Case: No special configuration needed
```python
method = "leastsq"  # Default, fast convergence
```

### Global Optimization: Multiple global methods supported
```python
methods = ["brute", "differential_evolution", "basinhopping", "ampgo"]
```

### Hierarchical: Only leastsq for sub-problems
```python
method = "leastsq"  # For nested optimization
```

---

## Optimization Method Features Required for Migration

| Method | Scope | Bounds | Constraints | Callback | Notes |
|--------|-------|--------|------------|----------|-------|
| leastsq | Local | Yes | Via expr | No | Via lmfit |
| brute | Global | Yes | Via expr | Limited | Grid-based |
| differential_evolution | Global | Yes | Yes | Limited | Scipy native |
| basinhopping | Global | Yes | Via inner method | Limited | Uses L-BFGS-B inside |
| ampgo | Global | Yes | Via expr | No | Less known |

---

## Parameter Bounds Usage

All methods respect parameter bounds defined in `ParamSetting`:
```python
class ParamSetting:
    min: float = -np.inf
    max: float = np.inf
```

These are passed to lmfit as:
```python
(param_id, value, vary, min, max, expr)
```

---

## Result Object Fields Used by ChemEx

From lmfit MinimizerResult:
```python
result.params          # Updated Parameters object
result.nfev            # Number of function evaluations
result.chisqr          # Chi-square value (sum of residuals^2)
result.redchi          # Reduced chi-square
result.residual        # Residual array (used in hierarchical fitting)
result.success         # Optimization success flag (checked implicitly)
```

---

## Special Features: Hierarchical Fitting

The hierarchical fitting uses a two-level optimization:

```
1. Top level: Optimize all parameters at once
   ├─ Calls residuals_hierarchical(params, param_tree)
   └─ Minimizer wraps this function
   
2. Inside residuals_hierarchical:
   ├─ If leaf node: return experiments.residuals(params)
   └─ If internal: recursively call lmfit.minimize() on each branch
       └─ Uses method="leastsq" only
       └─ Extracts result.residual from each branch
```

This requires:
- `lmfit.minimize()` function (not just Minimizer class)
- Support for `args=()` to pass additional arguments
- `.residual` attribute on result object

---

## Constraint Expression System

Constraints are specified as Python expressions with parameter substitution:

### Internal Format
```python
expr = "__param_id_1 + 2 * __param_id_2"
```

### Processing
```python
# 1. Build lmfit.Parameters with usersyms
usersyms = rate_functions | user_function_registry.get(model.name)
lmfit_params = ParametersLF(usersyms=usersyms)

# 2. Add parameters with constraint expressions
lmfit_params.add_many(
    (param_id_1, value1, True, min1, max1, ""),
    (param_id_2, value2, True, min2, max2, expr),
    ...
)

# 3. Process constraints
lmfit_params.update_constraints()
```

### Result
- Constrained parameters have `.expr` populated
- lmfit evaluates expressions using parameter values and `usersyms`
- Parameter `param_id_2` value is computed from `param_id_1` and user functions

### User Functions Available
- `rate_functions`: Dict of RatesIS callables for model-free analysis
- `user_function_registry`: Custom functions per kinetic model

---

## Callback Mechanism

Progress reporting via iteration callback:

```python
def iter_cb(params, iteration, residuals, **kwargs):
    """Called by lmfit at each iteration."""
    chisqr = (residuals**2).sum()
    nvarys = len([p for p in params.values() if p.vary and not p.expr])
    redchi = chisqr / (len(residuals) - nvarys)
    # Print progress...
```

**Requirements**:
- Called with `(params, iteration, residuals, **kwargs)`
- Can access parameter attributes (`.vary`, `.expr`)
- Can compute statistics from residuals
- Optional: return value ignored

---

## Statistics Post-Fit

After optimization, statistics are calculated:

```python
def calculate_statistics(experiments, params_lf):
    residuals = experiments.residuals(params_lf)
    ndata = len(residuals)
    nvarys = len([p for p in params_lf.values() if p.vary and not p.expr])
    
    chisqr = sum(residuals**2)
    redchi = chisqr / max(1, ndata - nvarys)
    aic = chisqr + 2 * nvarys
    bic = chisqr + np.log(ndata) * nvarys
    
    # Chi-square test
    pvalue = 1.0 - stats.chi2.cdf(chisqr, ndata - nvarys)
    
    # Kolmogorov-Smirnov test (residuals ~ normal distribution)
    _, ks_p_value = stats.kstest(residuals, "norm")
```

**Assumptions**:
- Residuals sum to chi-square
- Reduced chi-square uses (ndata - nvarys) as degrees of freedom
- Information criteria (AIC, BIC) for model comparison
- Statistical tests require fit quality assessment

---

## Grid Search Integration

Grid searches are performed by:

1. Creating parameter combinations
2. For each combination:
   - Set parameter values
   - Run minimization (using same method as normal fitting)
   - Record chi-square
   - Track best result

```python
for grid_point in grid_values:
    _set_param_values(params, param_ids, grid_point)
    optimized = minimize(experiments, params, fitmethod)
    chisqr = calculate_statistics(experiments, optimized)["chisqr"]
    if chisqr < best_chisqr:
        best_params = optimized
```

**Important**: Grid search uses the same method (leastsq by default)

---

## Summary of Integration Points

### 1. Parameter Creation
- Input: ParamSetting with (name, value, vary, min, max, expr)
- Process: Convert to lmfit.Parameters via `add_many()`
- Output: lmfit.Parameters ready for optimization

### 2. Optimization Loop
- Input: lmfit.Parameters + residual function
- Process: lmfit.Minimizer wrapper with optional callback
- Output: MinimizerResult with updated params, statistics

### 3. Parameter Update
- Input: MinimizerResult.params (updated lmfit.Parameters)
- Process: Extract value, stderr from each parameter
- Output: Update internal ParamSetting

### 4. Constraint Handling
- Input: expr string + usersyms dictionary
- Process: lmfit evaluates expressions during minimization
- Output: Constrained parameters have auto-updated values

### 5. Statistics Computation
- Input: Optimized residuals and parameter counts
- Process: Calculate chi-square, AIC, BIC, p-values
- Output: Statistics dict for reporting

---

## Critical Path for Migration

**Must Support**:
1. Minimizer class with fcn, params, fcn_args, iter_cb
2. Multiple optimization methods (at least leastsq, differential_evolution)
3. Parameter object with value, vary, expr, min, max, stderr
4. Constraint expression evaluation via usersyms
5. Result object with params, nfev, chisqr, redchi, residual
6. Nested minimization via lmfit.minimize() function

**Must Replicate**:
1. Callback mechanism for progress reporting
2. Constraint expression system with user functions
3. Standard error attachment to parameters
4. Grid search parameter reuse pattern
5. Hierarchical fitting two-level nesting

