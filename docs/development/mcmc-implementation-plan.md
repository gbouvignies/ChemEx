# MCMC Parameter Uncertainty Implementation Plan

## Goal

Add Markov chain Monte Carlo sampling as a first-class uncertainty analysis in
ChemEx, alongside the existing Monte Carlo and bootstrap statistics. The feature
should reuse the current fitting, parameter, and output infrastructure, keep the
initial implementation small, and leave a clear path for richer posterior
diagnostics later.

## Current Architecture

ChemEx already has the right extension points:

- `chemex.configuration.methods.Statistics` parses the `STATISTICS` method-file
  section and currently supports `MC`, `BS`, and `BSN`.
- `chemex.optimize.fitting._run_statistics()` is called after each normal fit and
  writes one statistics output file per enabled method.
- `chemex.optimize.minimizer.minimize()` and `minimize_with_report()` wrap
  `lmfit.Minimizer`.
- `chemex.optimize.helper.calculate_statistics()` computes fit statistics from
  the weighted residual vector.
- `chemex.parameters.database.ParameterStore.build_lmfit_params()` reconstructs
  the active `lmfit.Parameters`, including bounds, expressions, and current
  best-fit values.

The project already depends on `emcee>=3.1.6`, and `lmfit` provides an
`emcee()` wrapper that samples the posterior from the same residual function that
ChemEx already uses for deterministic fitting.

## Product Behavior

MCMC should be an uncertainty/statistics method, not a replacement for the normal
best-fit step. Users should still run a deterministic fit first, then sample the
posterior around that fitted solution.

Minimal method-file syntax:

```toml
[STEP1]
FITMETHOD = "leastsq"
STATISTICS = { "MCMC" = 5000 }
```

Expanded syntax for advanced users:

```toml
[STEP1.STATISTICS.MCMC]
STEPS = 5000
BURN = 1000
THIN = 10
WALKERS = 64
SEED = 1234
WORKERS = 1
```

The short form keeps parity with `MC`, `BS`, and `BSN`: the integer means
`steps`. The expanded form is needed for reproducibility and basic convergence
control.

## Configuration Design

Introduce a small Pydantic model:

```python
class McmcSettings(BaseModel):
    steps: PositiveInt
    burn: NonNegativeInt | Literal["auto"] = "auto"
    thin: PositiveInt = 1
    walkers: PositiveInt | None = None
    seed: int | None = None
    workers: PositiveInt = 1
    update_parameters: bool = False
```

Then change `Statistics` to accept either short-form or expanded MCMC settings:

```python
class Statistics(BaseModel):
    mc: PositiveInt | None = None
    bs: PositiveInt | None = None
    bsn: PositiveInt | None = None
    mcmc: PositiveInt | McmcSettings | None = None
```

Validation rules:

- Convert `STATISTICS = {"MCMC" = 5000}` to `McmcSettings(steps=5000)`.
- Require numeric `burn < steps`.
- Require enough retained samples after burn/thinning.
- If `walkers` is omitted, default to `max(32, 2 * nvarys)` at runtime, because
  the value depends on the current group.
- Default `burn` to `"auto"` and resolve it from the integrated
  autocorrelation time after the sampler has run. If autocorrelation time is not
  available or implies discarding the whole chain, keep the full chain and record
  a diagnostic warning.
- Default `thin` to `1`. Thinning is kept as an explicit storage/output-size
  control instead of an automatic convergence tool.
- Keep `update_parameters = false` initially so MCMC does not silently replace
  deterministic best-fit output with posterior medians. This makes the new
  feature observational by default and avoids surprising existing workflows.

## Sampling Architecture

Add a new module, `chemex.optimize.mcmc`, with these responsibilities:

- define the immutable runtime settings dataclass or Pydantic model used by the
  sampler;
- validate the active parameter set before sampling;
- run `lmfit.Minimizer(experiments.residuals, params).emcee(...)`;
- transform the `lmfit.MinimizerResult` into ChemEx-owned result objects;
- write summary, samples, and correlation outputs.

Keep MCMC-specific code out of `fitting.py` except for dispatch. A good public
surface is:

```python
def run_mcmc(
    experiments: Experiments,
    params: lmfit.Parameters,
    settings: McmcSettings,
    path: Path,
) -> McmcResult:
    ...
```

Use ChemEx-owned result containers instead of passing `lmfit.MinimizerResult`
through the application:

```python
@dataclass(frozen=True)
class McmcSummary:
    parameter_id: str
    mean: float
    standard_deviation: float
    median: float
    eti_95_lower: float
    eti_95_upper: float
    lower_1sigma: float
    upper_1sigma: float
    stderr: float
    effective_sample_size: float | None
    mcse_mean: float | None

@dataclass(frozen=True)
class McmcResult:
    var_names: tuple[str, ...]
    chain: Array
    lnprob: Array
    summary: tuple[McmcSummary, ...]
    correlations: Array
    acceptance_fraction: Array
    autocorrelation_time: Array | None
    discarded_steps: int
    burn_in_warning: str | None
```

This keeps external-library details localized and makes tests independent of
`lmfit` internals.

## Statistical Model

Use ChemEx's existing weighted residuals and call `lmfit.Minimizer.emcee()` with
`is_weighted=True`. The resulting log likelihood is consistent with the current
error model: residuals are already `(experiment - model) / error`.

Use existing `lmfit.Parameter.min` and `max` values as uniform priors. This is
the least surprising first implementation because users already set bounds in
parameter files:

```toml
[GLOBAL]
PB = [0.05, 0.0, 1.0]
KEX_AB = [200.0, 1.0, 5000.0]
```

Before sampling, warn when any varied parameter has an infinite bound. This
should not be a hard error in the first implementation, because many current
ChemEx defaults are one-sided. The warning should explain that finite bounds are
recommended for meaningful MCMC priors.

## Output Design

Create a dedicated `MCMC` directory under the current group or step output path:

```text
Output/
  STEP1/
    MCMC/
      summary.toml
      samples.out
      correlations.out
      diagnostics.toml
```

`summary.toml` should be human-readable and stable:

```toml
[GLOBAL.KEX_AB]
prior = "uniform"
prior_lower = 1.00000e+00
prior_upper = 5.00000e+03
credible_interval = "95% equal-tailed"
mean = 3.81511e+02
standard_deviation = 9.10234e+00
median = 3.81511e+02
eti_95_lower = 3.64001e+02
eti_95_upper = 3.99020e+02
lower_1sigma = 3.72602e+02
upper_1sigma = 3.90420e+02
stderr = 8.90870e+00
effective_sample_size = 1.25000e+03
mcse_mean = 2.57446e-01

[GLOBAL.PB]
prior = "uniform"
prior_lower = 0.00000e+00
prior_upper = 1.00000e+00
credible_interval = "95% equal-tailed"
mean = 7.02971e-02
standard_deviation = 1.18903e-03
median = 7.02971e-02
eti_95_lower = 6.80343e-02
eti_95_upper = 7.25599e-02
lower_1sigma = 6.91493e-02
upper_1sigma = 7.14449e-02
stderr = 1.14780e-03
effective_sample_size = 1.18000e+03
mcse_mean = 3.46109e-05
```

`samples.out` should use the same lightweight text-table style as existing
statistics output:

```text
# KEX_AB PB [lnprob]
  3.81511e+02  7.02971e-02 -1.23456e+03
```

`correlations.out` should be a square matrix with parameter names in a header.

`diagnostics.toml` should include:

- `steps`
- requested burn-in setting and actual discarded steps
- `thin`
- `walkers`
- number of retained samples
- mean/min/max acceptance fraction
- autocorrelation time when available
- ESS and retained chain length relative to autocorrelation time when available
- a warning field when autocorrelation could not be estimated

Avoid adding plots in the first PR. Chain and corner plots are valuable, but they
introduce heavier dependencies and test surface. Text outputs are enough for the
initial feature.

## Integration Points

Change `_run_statistics()` in `chemex.optimize.fitting` to separate resampling
statistics from posterior sampling:

- keep the existing `MC`, `BS`, and `BSN` loop unchanged;
- dispatch `MCMC` to `chemex.optimize.mcmc.run_mcmc()`;
- do not call `generate_exp_for_statistics()` for MCMC;
- pass the current best-fit parameters from `parameter_store.build_lmfit_params()`.

The call should happen after `execute_post_fit()` and after
`parameter_store.update_from_parameters(best_lmfit_params)`, which is already
the current flow. That ensures the sampler starts from the fitted solution.

## Error Handling

Use specific validation and user-facing errors where possible:

- no varied parameters: skip MCMC with a clear warning;
- invalid burn/thin settings: fail method parsing before fitting starts;
- insufficient finite bounds: warn, do not fail in the first implementation;
- sampler interruption: handle `KeyboardInterrupt` like current statistics;
- `ValueError` from `lmfit` or `emcee`: print ChemEx's existing value-error
  message and leave partial outputs flushed when possible.

Avoid catching broad `Exception`; sampler failures should be visible during
development and CI.

## Testing Plan

Add focused unit tests rather than expensive end-to-end MCMC tests:

- configuration parsing:
  - short-form `STATISTICS = {"MCMC" = 100}`;
  - expanded `[STEP.STATISTICS.MCMC]`;
  - lower-case key coercion;
  - invalid `burn >= steps`;
- runtime defaults:
  - default walkers are derived from number of varied parameters;
  - no varied parameters is handled without invoking the sampler;
- writer tests:
  - `summary.toml`, `samples.out`, `correlations.out`, and `diagnostics.toml`
    formatting from a small synthetic `McmcResult`;
- integration dispatch:
  - `_run_statistics()` calls MCMC without generating bootstrap/Monte Carlo
    experiments;
  - existing `MC`, `BS`, and `BSN` behavior remains unchanged.

For one smoke test, monkeypatch the actual sampler call to return a tiny fixed
chain. This avoids making CI slow or nondeterministic while still validating the
ChemEx integration.

## Documentation Plan

Update the user guide:

- add MCMC to the "Estimating Parameter Uncertainty" section;
- explain that MCMC samples the posterior after a normal fit;
- document the short and expanded syntax;
- recommend finite parameter bounds;
- document output files and diagnostics;
- clarify that MCMC uncertainty is posterior-based, while `MC`/`BS` refit
  synthetic datasets.

Update `messages.py` help text if needed so `MCMC` appears alongside the
statistics methods, not only as `FITMETHOD = "emcee"`.

## Implementation Phases

1. Add configuration models and tests.
2. Add `chemex.optimize.mcmc` result containers, sampler wrapper, validation,
   and output writers.
3. Wire MCMC into `_run_statistics()`.
4. Add focused tests for dispatch and output generation.
5. Update website documentation.
6. Run `pytest` for the touched tests, then the broader suite if runtime is
   acceptable.

## Non-Goals For The First Implementation

- custom priors beyond current bounds;
- posterior plotting;
- resumable chains;
- replacing best-fit parameters by default;
- multiprocessing defaults above one worker;
- adding `pandas` solely for `result.flatchain`.

These can be added later without changing the core architecture if the first
implementation keeps MCMC isolated behind ChemEx-owned result objects.
