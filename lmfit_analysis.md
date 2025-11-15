# ChemEx lmfit Usage Analysis - Migration Planning Guide

## Summary

ChemEx uses **lmfit >= 1.3.2** across 7 critical files for parameter optimization and management. The integration is well-encapsulated through a parameter database system, which significantly simplifies migration potential.

---

## 1. FILES USING LMFIT (By Importance)

### 1.1 Core Optimization Files

#### **`src/chemex/optimize/minimizer.py`** (CRITICAL)
- **Purpose**: Main optimization engine with three fitting strategies
- **lmfit Features Used**:
  - `lmfit.Minimizer` class (2 main uses)
  - `lmfit.minimize()` function (1 use in hierarchical fitting)
  - `lmfit.Parameters` objects
  - Iteration callbacks (`iter_cb` parameter)
  - Result object attributes: `.result.params`, `.result.nfev`, `.result.chisqr`, `.result.redchi`, `.residual`

**Key Functions**:
```python
def minimize(experiments, params: lmfit.Parameters, fitmethod: str) -> lmfit.Parameters:
    # Simple fit without callbacks
    minimizer = lmfit.Minimizer(experiments.residuals, params)
    minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    return deepcopy(minimizer.result.params)

def minimize_with_report(experiments, params: lmfit.Parameters, fitmethod: str) -> lmfit.Parameters:
    # Fit with progress reporting via iter_cb callback
    reporter = Reporter()
    minimizer = lmfit.Minimizer(experiments.residuals, params, iter_cb=reporter.iter_cb)
    result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    # Uses: result.nfev, result.chisqr, result.redchi
    
def minimize_hierarchical(experiments, params: lmfit.Parameters, fitmethod: str):
    # Hierarchical fitting with fcn_args
    minimizer = lmfit.Minimizer(residuals_hierarchical, params, 
                               fcn_args=(param_tree,), iter_cb=reporter.iter_cb)
    result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))

def residuals_hierarchical(params: lmfit.Parameters, param_tree: ParamTree) -> Array:
    # Uses lmfit.minimize() for sub-optimization
    branch_results = lmfit.minimize(residuals_hierarchical, params, 
                                   args=(branch,), method="leastsq")
    # Uses: branch_results.residual
```

**Optimization Methods Used**:
- `"brute"` - with `{"keep": "all"}` kwargs
- `"differential_evolution"` - with `{"disp": True}`
- `"basinhopping"` - with `{"disp": True, "niter_success": 10, "minimizer_kwargs": {"method": "L-BFGS-B"}}`
- `"ampgo"` - with `{"disp": True}`
- `"leastsq"` - basic least squares

#### **`src/chemex/parameters/database.py`** (CRITICAL)
- **Purpose**: Parameter management and lmfit Parameter construction
- **lmfit Features Used**:
  - `lmfit.Parameters` class (imported as `ParametersLF`)
  - Building Parameters from internal representation: `ParametersLF(usersyms=usersyms)`
  - `lmfit.Parameters.add_many()` method
  - `lmfit.Parameters.update_constraints()` method
  - Parameter object attributes: `.value`, `.stderr`, `.vary`, `.expr`, `.min`, `.max`
  - `lmfit.Parameters.valuesdict()` method for converting to dict

**Key Method**:
```python
def build_lmfit_params(self, param_ids: Iterable[str] | None = None) -> ParametersLF:
    """Constructs lmfit Parameters from catalog's parameters."""
    parameters = self.get_parameters(param_ids)
    parameter_args = (parameter.args for parameter in parameters.values())
    
    # Uses rate_functions and user_function_registry as usersyms
    usersyms = rate_functions | user_function_registry.get(model.name)
    lmfit_params = ParametersLF(usersyms=usersyms)
    lmfit_params.add_many(*parameter_args)
    lmfit_params.update_constraints()
    
    for id_, lf_param in lmfit_params.items():
        lf_param.stderr = parameters[id_].stderr
    
    return lmfit_params

def update_from_lmfit_params(self, parameters: ParametersLF) -> None:
    """Updates internal parameters from lmfit Parameters after fitting."""
    for param_id, parameter in parameters.items():
        self._parameters[param_id].value = parameter.value
        self._parameters[param_id].stderr = parameter.stderr
```

---

### 1.2 Supporting Optimization Files

#### **`src/chemex/optimize/gridding.py`** (IMPORTANT)
- **Purpose**: Grid search optimization
- **lmfit Features Used**:
  - `lmfit.parameter.Parameters` type imported
  - Parameter attribute access: `.value` (setting and reading)
- **Usage**: `_set_param_values()` function sets parameter values in grid loops

#### **`src/chemex/optimize/helper.py`** (IMPORTANT)
- **Purpose**: Post-fit statistics and utilities
- **lmfit Features Used**:
  - `lmfit.Parameters` type annotations
  - Parameter object inspection: `.vary`, `.expr` attributes
  - Counting variable parameters: `len([param for param in params_lf.values() if param.vary and not param.expr])`

**Key Function**:
```python
def calculate_statistics(experiments: Experiments, params_lf: ParametersLF) -> dict[str, int | float]:
    """Calculates fit statistics including chi-square, AIC, BIC."""
    residuals = experiments.residuals(params_lf)
    ndata = len(residuals)
    nvarys = len([param for param in params_lf.values() if param.vary and not param.expr])
    chisqr = sum(residuals**2)
    redchi = chisqr / max(1, ndata - nvarys)
    aic = chisqr + 2 * nvarys
    bic = chisqr + np.log(ndata) * nvarys
    pvalue = 1.0 - stats.chi2.cdf(chisqr, ndata - nvarys)
    return {"ndata": ndata, "nvarys": nvarys, "chisqr": chisqr, "redchi": redchi, "pvalue": pvalue, ...}
```

---

### 1.3 Data Container Files

#### **`src/chemex/containers/experiment.py`** (MODERATE)
#### **`src/chemex/containers/experiments.py`** (MODERATE)
#### **`src/chemex/containers/profile.py`** (MODERATE)

- **Purpose**: Data structures for experiments and profiles
- **lmfit Features Used**:
  - `lmfit.Parameters` type annotations
  - Parameter value extraction: `params[param_id].value`
  - Used in residual calculation functions

**Usage Pattern**:
```python
def residuals(self, params: ParametersLF) -> Array:
    """Calculate residuals given lmfit Parameters."""
    return np.concatenate([profile.residuals(params) for profile in self.profiles])

def _get_parameter_values(self, params: ParametersLF) -> dict[str, float]:
    """Extract parameter values from lmfit Parameters object."""
    return {
        local_name: params[param_id].value
        for local_name, param_id in self.name_map.items()
    }
```

---

## 2. LMFIT INTEGRATION ARCHITECTURE

### Data Flow Diagram:
```
┌─────────────────────────────────────┐
│ Parameter Database (database.py)    │
│ - Maintains internal ParamSetting   │
│ - Builds lmfit.Parameters objects   │
│ - Updates from lmfit.Parameters     │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│ lmfit.Parameters                    │
│ - name, value, vary, expr, min, max │
│ - stderr (set from database)         │
│ - usersyms (rate functions)          │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│ lmfit.Minimizer                     │
│ - Objective: experiments.residuals()│
│ - Methods: brute, leastsq, etc.     │
│ - Callbacks: iter_cb for progress   │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│ lmfit.Result                        │
│ - .params (updated Parameters)       │
│ - .nfev, .chisqr, .redchi           │
│ - .residual (for hierarchical fit)  │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│ Database Update                     │
│ - Extract value and stderr          │
│ - Update internal ParamSetting      │
└─────────────────────────────────────┘
```

---

## 3. CRITICAL LMFIT FEATURES REQUIRED FOR MIGRATION

### 3.1 Parameter Management
- **Parameter object with attributes**: `value`, `vary`, `expr`, `min`, `max`, `stderr`, `name`
- **Parameter collection**: Dictionary-like access by parameter ID
- **Expression constraints**: Support for constraint expressions like `"[param1] + 2*[param2]"`
- **User-defined functions**: `usersyms` dict for constraint expression evaluation
- **Constraint resolution**: `update_constraints()` to process expressions

### 3.2 Minimizer Class
- Constructor: `Minimizer(fcn, params, fcn_args=(), iter_cb=None, ...)`
- Method: `.minimize(method=str, **kwargs)` supporting:
  - `"brute"`: Brute force grid search
  - `"leastsq"`: Levenberg-Marquardt algorithm
  - `"differential_evolution"`: Global optimization
  - `"basinhopping"`: Basin-hopping global optimizer
  - `"ampgo"`: Autonomous Metropolis-hopping global optimizer
- Result object: `.result.params`, `.result.nfev`, `.result.chisqr`, `.result.redchi`, `.result.residual`
- Callbacks: `iter_cb(params, iteration, residuals, **kwargs)` for progress reporting

### 3.3 Free Function API
- `lmfit.minimize(fcn, params, args=(), method='leastsq', ...)` for simple minimization

### 3.4 Objective Function Interface
- **Input**: lmfit.Parameters object
- **Output**: 1D numpy array of residuals
- **Calling**: Typically called as `experiments.residuals(params)`

### 3.5 Statistics Computation
- Chi-square: `sum(residuals**2)`
- Reduced chi-square: `chisqr / (ndata - nvarys)`
- Counting variable parameters: Count params with `.vary=True` and `.expr=""`
- Standard errors: Must be computed and attached to parameters

---

## 4. INTEGRATION POINTS & DEPENDENCIES

### 4.1 Parameter-Lmfit Bidirectional Binding
**In `database.py`**:
```python
# Forward: Internal ParamSetting → lmfit.Parameters
lmfit_params.add_many(*parameter_args)  # parameter.args tuple: (id, value, vary, min, max, expr)

# Backward: lmfit.Parameters → Internal ParamSetting
self._parameters[param_id].value = parameter.value
self._parameters[param_id].stderr = parameter.stderr
```

### 4.2 User Functions for Constraints
```python
usersyms = rate_functions | user_function_registry.get(model.name)
lmfit_params = ParametersLF(usersyms=usersyms)
```
- `rate_functions`: Dict of RatesIS callables (defined in `nmr/rates.py`)
- `user_function_registry`: Custom user functions per kinetic model

### 4.3 Hierarchical Fitting Pattern
**Two-level optimization**:
```python
def residuals_hierarchical(params: lmfit.Parameters, param_tree: ParamTree) -> Array:
    if not param_tree.branches:
        return param_tree.experiments.residuals(params)  # Leaf node
    
    # For branches, recursively minimize each child
    for branch in param_tree.branches:
        branch_results = lmfit.minimize(residuals_hierarchical, params, 
                                       args=(branch,), method="leastsq")
        residuals_list.append(branch_results.residual)  # Needs .residual attribute
```

---

## 5. DETAILED USAGE PATTERNS

### 5.1 Basic Fitting Workflow
```
1. database.build_lmfit_params(param_ids) 
   → Constructs lmfit.Parameters from internal database
   
2. minimize_with_report(experiments, lmfit_params, fitmethod)
   → Creates Minimizer(experiments.residuals, lmfit_params, iter_cb=reporter.iter_cb)
   → Calls .minimize(method=fitmethod, **method_kwargs)
   
3. database.update_from_parameters(result.params)
   → Extracts .value and .stderr from lmfit.Parameters
   → Updates internal ParamSetting objects
```

### 5.2 Callback Pattern
```python
reporter.iter_cb(params: lmfit.Parameters, iteration: int, residuals: Array, **kwargs)
→ Calculates chi-square = sum(residuals**2)
→ Counts nvarys from params with param.vary=True and param.expr==""
→ Calculates redchi = chisqr / (ndata - nvarys)
→ Prints progress table
```

### 5.3 Constraint Expression Syntax
```
Internal representation: "param_id_1 + 2 * param_id_2"
After formatting: "__param_id_1 + 2 * __param_id_2"
Uses lmfit's double-underscore prefix convention for parameter names in expressions
```

### 5.4 Grid Search Integration
```python
# For each grid point:
params["param_id"].value = grid_value
optimized_params = minimize(experiments, params, fitmethod)
# Result parameters returned with optimized .value and .stderr
```

---

## 6. POTENTIAL MIGRATION CHALLENGES

### 6.1 Complex Parameter Expressions
- **Challenge**: Constraint expressions are evaluated by lmfit via `update_constraints()`
- **Risk**: scipy.optimize has no native constraint expression support
- **Solution Required**: Custom expression parser and evaluator

### 6.2 Hierarchical Fitting Structure
- **Challenge**: Two-level recursive optimization with intermediate result extraction
- **Risk**: scipy.optimize doesn't support nested minimization patterns
- **Solution Required**: Custom optimization controller to manage nesting

### 6.3 User-Defined Functions in Constraints
- **Challenge**: Rate function callables passed as `usersyms`
- **Risk**: scipy.optimize doesn't evaluate symbolic expressions
- **Solution Required**: Custom constraint resolver that calls user functions

### 6.4 Multiple Optimization Methods
- **Challenge**: Mixing global optimizers (brute, differential_evolution, basinhopping) with local (leastsq)
- **Risk**: scipy provides different APIs for global vs local
- **Solution Required**: Unified wrapper interface

### 6.5 Iteration Callbacks
- **Challenge**: Progress reporting via `iter_cb` callback
- **Risk**: Not all scipy.optimize methods support callbacks
- **Solution Required**: Method-specific callback implementations

### 6.6 Result Object Structure
- **Challenge**: Accessing `.nfev`, `.chisqr`, `.redchi`, `.residual` from result
- **Risk**: scipy.optimize.Result has different structure
- **Solution Required**: Result wrapper/adapter

### 6.7 Standard Error Estimation
- **Challenge**: lmfit automatically computes `.stderr` for parameters
- **Risk**: scipy.optimize doesn't compute uncertainties
- **Solution Required**: Post-fit uncertainty estimation (using Hessian, Monte Carlo, or bootstrap)

---

## 7. KEY FINDINGS FOR REPLACEMENT STRATEGY

### Quick wins with scipy.optimize:
- ✅ `"leastsq"` → `scipy.optimize.least_squares()` with callback support
- ✅ `"differential_evolution"` → `scipy.optimize.differential_evolution()` with direct mapping
- ✅ Basic parameter bounds mapping

### High-effort migrations with scipy.optimize:
- ❌ Constraint expressions with `usersyms` evaluation
- ❌ Hierarchical fitting framework
- ❌ `"brute"` grid search (scipy has limited brute support)
- ❌ Callback uniformity across methods

### JAX Advantages:
- ✅ Can auto-compute gradients for all methods (unlike scipy)
- ✅ Can auto-compute Hessian for uncertainty estimation
- ✅ Can JIT-compile constraint expressions
- ✅ Easier to implement custom nested optimization

---

## 8. VERSION AND DEPENDENCIES

**Current lmfit requirement**: `"lmfit>=1.3.2"`

### Key lmfit API features in use:
- lmfit 1.3.2+ stable API (Minimizer, Parameters, minimize function)
- No reliance on deprecated features
- All used features stable across lmfit 1.x versions

---

## 9. SUMMARY TABLE

| File | Impact | Key Feature | Difficulty |
|------|--------|-------------|------------|
| minimizer.py | CRITICAL | Minimizer class, callbacks, result attributes | High |
| database.py | CRITICAL | Parameter construction, constraints, usersyms | High |
| helper.py | IMPORTANT | Parameter inspection (vary, expr) | Low |
| gridding.py | IMPORTANT | Parameter value setting | Low |
| experiment.py | MODERATE | Residual function interface | Low |
| experiments.py | MODERATE | Residual aggregation | Low |
| profile.py | MODERATE | Parameter value extraction | Low |

---

## 10. RECOMMENDED NEXT STEPS

1. **Create abstraction layer** in `optimize/minimizer.py` to decouple lmfit API
2. **Implement JAX-based replacement** with autodiff for uncertainty estimation
3. **Test constraint expression evaluation** with model-free analysis
4. **Validate hierarchical fitting** with nested optimization
5. **Benchmark performance** on real fitting workflows
6. **Plan phased migration** (leastsq first, then other methods)

