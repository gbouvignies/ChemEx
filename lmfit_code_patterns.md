# ChemEx lmfit Code Patterns - Detailed Implementation Reference

## File: `/home/user/ChemEx/src/chemex/optimize/minimizer.py`

### Pattern 1: Simple Minimization (Lines 60-71)
```python
def minimize(
    experiments: Experiments,
    params: lmfit.Parameters,
    fitmethod: str,
) -> lmfit.Parameters:
    """Basic optimization without callbacks."""
    kws = {
        "brute": {"keep": "all"},
    }
    
    minimizer = lmfit.Minimizer(experiments.residuals, params)
    minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
    return deepcopy(minimizer.result.params)
```

**Key aspects**:
- Objective function: `experiments.residuals` (callable taking params, returning array)
- Method-specific kwargs: `brute` uses `{"keep": "all"}`
- Result extraction: `minimizer.result.params` (updated Parameters object)

---

### Pattern 2: Minimization with Progress Callbacks (Lines 74-104)
```python
class Reporter:
    last_chisqr: float = +1.0e32
    threshold: float = -1.0e-3
    
    def iter_cb(
        self,
        params: lmfit.Parameters,
        iteration: int,
        residuals: Array,
        *_args,
        **_kwargs,
    ) -> None:
        """Called at each iteration by lmfit.Minimizer."""
        chisqr = (residuals**2).sum()
        change = (chisqr - self.last_chisqr) / self.last_chisqr
        
        if change > self.threshold or iteration < 0:
            return
        
        self.last_chisqr = chisqr
        
        ndata = len(residuals)
        nvarys = len(
            [param for param in params.values() if param.vary and not param.expr],
        )
        redchi = chisqr / max(1, ndata - nvarys)
        
        self.print_line(iteration, chisqr, redchi)

def minimize_with_report(
    experiments: Experiments,
    params: lmfit.Parameters,
    fitmethod: str,
) -> lmfit.Parameters:
    """Optimization with progress reporting."""
    kws = {
        "brute": {"keep": "all"},
        "differential_evolution": {"disp": True},
        "basinhopping": {
            "disp": True,
            "niter_success": 10,
            "minimizer_kwargs": {"method": "L-BFGS-B"},
        },
        "ampgo": {"disp": True},
    }
    
    reporter = Reporter()
    
    minimizer = lmfit.Minimizer(
        experiments.residuals, 
        params, 
        iter_cb=reporter.iter_cb
    )
    
    reporter.print_header()
    
    try:
        result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
        reporter.print_footer(result.nfev, result.chisqr, result.redchi)
    except KeyboardInterrupt:
        print_calculation_stopped_error()
    except ValueError:
        print_value_error()
    
    return deepcopy(minimizer.result.params)
```

**Key aspects**:
- Callback signature: `iter_cb(params, iteration, residuals, *args, **kwargs)`
- Access to result attributes: `.nfev`, `.chisqr`, `.redchi`
- Parameter inspection: `param.vary`, `param.expr` to count variable params

---

### Pattern 3: Hierarchical Fitting (Lines 107-167)
```python
def residuals_hierarchical(
    params: lmfit.Parameters, 
    param_tree: ParamTree
) -> Array:
    """Recursive residual function for hierarchical fitting."""
    if not param_tree.branches:
        # Leaf node: calculate residuals directly
        return param_tree.experiments.residuals(params)
    
    # Internal node: recursively minimize each branch
    for param_id in param_tree.ids_to_fit:
        params[param_id].vary = False
    
    residuals_list: list[Array] = []
    
    for branch in param_tree.branches:
        # Call lmfit.minimize() to optimize each branch
        branch_results = lmfit.minimize(
            residuals_hierarchical,
            params,
            args=(branch,),  # Pass branch as additional argument
            method="leastsq",
        )
        # Extract residuals from result object
        residuals_list.append(branch_results.residual)
    
    residuals = np.concatenate(residuals_list, axis=0)
    
    for param_id in param_tree.ids_to_fit:
        params[param_id].vary = True
    
    return residuals

def minimize_hierarchical(
    experiments: Experiments, 
    params: lmfit.Parameters, 
    fitmethod: str
):
    """Top-level hierarchical fitting orchestrator."""
    kws = {
        "brute": {"keep": "all"},
        "ampgo": {"disp": True},
    }
    
    reporter = Reporter()
    param_tree = create_group_tree(experiments)
    
    new_params = params.copy()
    for param in new_params.values():
        if param.name not in param_tree.ids_to_fit:
            param.vary = False
    
    minimizer = lmfit.Minimizer(
        residuals_hierarchical,
        params,
        fcn_args=(param_tree,),  # Pass param_tree to objective function
        iter_cb=reporter.iter_cb,
    )
    
    reporter.print_header()
    
    try:
        result = minimizer.minimize(method=fitmethod, **(kws.get(fitmethod, {})))
        reporter.print_footer(result.nfev, result.chisqr, result.redchi)
    except KeyboardInterrupt:
        print_calculation_stopped_error()
    except ValueError:
        print_value_error()
    
    return deepcopy(minimizer.result.params)
```

**Key aspects**:
- `Minimizer(fcn, params, fcn_args=(data,), iter_cb=callback)`
- Objective function receives `fcn_args` as additional positional arguments
- `lmfit.minimize()` function used for nested optimization
- Result object has `.residual` attribute for nested access

---

## File: `/home/user/ChemEx/src/chemex/parameters/database.py`

### Pattern 4: Building lmfit Parameters (Lines 181-213)
```python
def build_lmfit_params(
    self,
    param_ids: Iterable[str] | None = None,
) -> ParametersLF:
    """Construct lmfit Parameters from the catalog's parameters.
    
    This method converts the internal parameter representation to the format
    required by the lmfit library for optimization.
    """
    if param_ids is None:
        param_ids = set(self._parameters)
    
    parameters = self.get_parameters(param_ids)
    
    # Generate args tuples: (name, value, vary, min, max, expr)
    parameter_args = (parameter.args for parameter in parameters.values())
    
    # Get user-defined functions for constraint evaluation
    usersyms = rate_functions | user_function_registry.get(model.name)
    
    # Create lmfit Parameters object
    lmfit_params = ParametersLF(usersyms=usersyms)
    
    # Add all parameters at once
    lmfit_params.add_many(*parameter_args)
    
    # Process constraint expressions
    lmfit_params.update_constraints()
    
    # Copy standard errors from internal database
    for id_, lf_param in lmfit_params.items():
        lf_param.stderr = parameters[id_].stderr
    
    return lmfit_params
```

**Key aspects**:
- Parameter args tuple: `(id, value, vary, min, max, expr)`
- `usersyms` dict: Contains callable functions for constraints
- `ParametersLF(usersyms=usersyms)` initialization
- `.add_many(*args_tuples)` to add multiple parameters
- `.update_constraints()` to evaluate expressions
- `.stderr` attribute set manually after creation

---

### Pattern 5: Updating Database from Fitted Parameters (Lines 215-224)
```python
def update_from_lmfit_params(self, parameters: ParametersLF) -> None:
    """Update catalog parameters from lmfit Parameters.
    
    Args:
        parameters (ParametersLF): lmfit Parameters to update from.
    """
    for param_id, parameter in parameters.items():
        self._parameters[param_id].value = parameter.value
        self._parameters[param_id].stderr = parameter.stderr
```

**Key aspects**:
- Iterating over lmfit Parameters: `parameters.items()` or `.values()`
- Extracting: `.value`, `.stderr`
- Writing back to internal database

---

### Pattern 6: Parameter Properties (Lines 128-129 in setting.py)
```python
@property
def args(self) -> tuple[str, float | None, bool | None, float, float, str]:
    """Return args tuple for lmfit Parameter addition."""
    return self.id_, self.value, self.vary, self.min, self.max, self.expr
```

**Key aspects**:
- Format: `(name, value, vary, min, max, expr)`
- Passed to `lmfit.Parameters.add_many()`

---

## File: `/home/user/ChemEx/src/chemex/optimize/helper.py`

### Pattern 7: Statistics Calculation (Lines 22-46)
```python
def calculate_statistics(
    experiments: Experiments,
    params_lf: ParametersLF,
) -> dict[str, int | float]:
    """Calculate fit statistics."""
    residuals = experiments.residuals(params_lf)
    ndata = len(residuals)
    
    # Count variable parameters (excluding constrained ones)
    nvarys = len(
        [param for param in params_lf.values() if param.vary and not param.expr],
    )
    
    chisqr = sum(residuals**2)
    redchi = chisqr / max(1, ndata - nvarys)
    aic = chisqr + 2 * nvarys
    bic = chisqr + np.log(ndata) * nvarys
    _, ks_p_value = stats.kstest(residuals, "norm")
    pvalue: float = 1.0 - stats.chi2.cdf(chisqr, ndata - nvarys)
    
    return {
        "ndata": ndata,
        "nvarys": nvarys,
        "chisqr": chisqr,
        "redchi": redchi,
        "pvalue": pvalue,
        "ks_pvalue": ks_p_value,
        "aic": aic,
        "bic": bic,
    }
```

**Key aspects**:
- Parameter inspection: `param.vary and not param.expr`
- Custom statistics: AIC, BIC, p-values

---

## File: `/home/user/ChemEx/src/chemex/optimize/gridding.py`

### Pattern 8: Grid Parameter Setting (Lines 37-43)
```python
def _set_param_values(
    params: Parameters,
    fnames: Iterable[str],
    values: tuple[float, ...],
) -> None:
    """Set parameter values in the Parameters object."""
    for fname, value in zip(fnames, values, strict=True):
        params[fname].value = value
```

**Key aspects**:
- Dictionary-like access: `params[param_id]`
- Direct `.value` assignment

---

### Pattern 9: Grid Search Workflow (Lines 46-92)
```python
def run_group_grid(
    group: Group,
    grid: dict[str, Array],
    path: Path,
    fitmethod: str,
) -> GridResult:
    """Run grid search for a group of experiments."""
    group_ids = group.experiments.param_ids
    group_params = database.build_lmfit_params(group_ids)
    group_grid = {
        param_id: values for param_id, values in grid.items() 
        if param_id in group_ids
    }
    
    best_chisqr = np.inf
    best_params = group_params
    
    grid_values = product(*group_grid.values())
    
    for values in track(grid_values, total=float(grid_size)):
        # Set parameter values for this grid point
        _set_param_values(group_params, grid_ids, values)
        
        # Optimize at this grid point
        optimized_params = minimize(group.experiments, group_params, fitmethod)
        
        # Calculate statistics
        stats = calculate_statistics(group.experiments, optimized_params)
        chisqr: float = stats.get("chisqr", np.inf)
        
        if chisqr < best_chisqr:
            best_chisqr = chisqr
            best_params = optimized_params
    
    # Update database with best result
    database.update_from_parameters(best_params)
    
    return GridResult(group_grid, chisqr_array)
```

**Key aspects**:
- Build Parameters once: `database.build_lmfit_params()`
- Reuse Parameter object, only change values
- Re-fit at each grid point
- Track best result

---

## File: `/home/user/ChemEx/src/chemex/containers/profile.py`

### Pattern 10: Parameter Value Extraction (Lines 71-76)
```python
def _get_parameter_values(self, params: ParametersLF) -> dict[str, float]:
    """Get the parameter values from the provided parameters."""
    return {
        local_name: params[param_id].value
        for local_name, param_id in self.name_map.items()
    }
```

**Key aspects**:
- Dictionary-like access for parameter values
- Mapping from local names to global parameter IDs

---

### Pattern 11: Residual Calculation (Lines 83-99)
```python
def calculate(self, params: ParametersLF) -> Array:
    """Calculate and return the Array."""
    self.update_spectrometer(params)
    self.data.calc_unscaled = self.pulse_sequence.calculate(
        self.spectrometer, self.data
    )
    if self.is_scaled:
        self.data.calc = self.data.scale * self.data.calc_unscaled
    else:
        self.data.calc = self.data.calc_unscaled
    return self.data.calc

@cachedmethod(attrgetter("cache"), key=_cache_key)
def residuals(self, params: ParametersLF) -> Array:
    """Calculate and return residuals."""
    residuals = (self.calculate(params) - self.data.exp) / self.data.err
    return residuals[self.data.mask]
```

**Key aspects**:
- Parameters passed to calculation functions
- Caching based on parameter values
- Return residuals as numpy array

---

## Key Integration Points

### 1. Parameter Object Interface
```
Attribute          Type      Usage
name               str       Parameter identifier
value              float     Current value
vary               bool      Should be optimized?
expr               str       Constraint expression (if any)
min                float     Lower bound
max                float     Upper bound
stderr             float     Standard error (post-fit)
```

### 2. Minimizer Constructor
```python
Minimizer(fcn, params, fcn_args=(), iter_cb=None, **kwargs)

fcn:        Callable[[Parameters, ...] → Array] residual function
params:     lmfit.Parameters object with initial values
fcn_args:   Tuple of additional arguments passed to fcn
iter_cb:    Callable[[Parameters, int, Array, **kwargs] → None] callback
```

### 3. Minimize Method
```python
minimizer.minimize(method='leastsq', **method_kwargs)

Returns: MinimizerResult with attributes:
  - params: Updated Parameters object
  - nfev: Number of function evaluations
  - chisqr: Chi-square value
  - redchi: Reduced chi-square
  - residual: Residual array (for nested calls)
  - success: Whether optimization succeeded
```

### 4. Free Function API
```python
lmfit.minimize(fcn, params, args=(), method='leastsq', **kwargs)
→ Returns MinimizerResult object
```

### 5. usersyms Dictionary
```python
usersyms = {
    'rate_function_name': callable_that_takes_params,
    'another_function': callable,
    ...
}
# Used in constraint expression evaluation
# E.g., expr = "rate_function_name(a, b, c)"
```

---

## Critical Functions for Migration

| Function | File | Purpose | Difficulty |
|----------|------|---------|------------|
| `Minimizer()` constructor | minimizer.py | Core optimization setup | High |
| `.minimize()` method | minimizer.py | Run optimization | High |
| `lmfit.minimize()` | minimizer.py | Nested optimization | High |
| `build_lmfit_params()` | database.py | Create Parameter objects | High |
| `update_constraints()` | database.py | Process expressions | High |
| `iter_cb()` callback | minimizer.py | Progress reporting | Medium |
| `calculate_statistics()` | helper.py | Post-fit analysis | Low |
| Parameter value setting | gridding.py | Grid search | Low |

