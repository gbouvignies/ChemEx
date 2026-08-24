---
sidebar_position: 6
---

# Method Files

## Overview

Method files define the fitting methods used in ChemEx. In these files, you can:

-   Specify which parameters to fit, fix, or constrain.
-   Select the profiles to include in calculations.
-   Activate additional analyses, such as grid search, Monte Carlo, bootstrap,
    or MCMC analyses.

Provide the method file to ChemEx using the `-m` or `--method` option.

`FITMETHOD` is deliberately closed to the native trust-region reflective
solver. Omit it or set `FITMETHOD = "trf"`. The spelling `"least_squares"` is
accepted temporarily as an alias and is canonicalized to `"trf"`; historical
optimizer names fail during method-file validation, before fit outputs are
invalidated. ChemEx no longer forwards arbitrary optimizer names or depends on
lmfit at runtime; fitting and statistics use the native ChemEx parameter and
evaluation stack. When `GRID` and `STATISTICS` occur in the same step, ChemEx
first commits the accepted aggregate grid fit and then runs the requested
statistics from that fit. This canonical workflow applies to both legacy and
version 2 method files.

Method files are structured in sections, each representing a fitting step.
Steps execute in the order they appear. Files without `FORMAT_VERSION` use the
legacy version 1 behavior: profile selection and parameter roles implicitly
carry forward when omitted. Version 2 makes role inheritance explicit and
treats profile selection as local to each step.

:::tip
When performing multi-step fitting, start with a subset of residues with high-quality data (e.g., large CPMG dispersion or clear CEST dips) to estimate global parameters, then fix these parameters in subsequent steps to fit residue-specific parameters. For example, see the method file for `CPMG_CH3_1H_SQ/` in `Examples/Experiments/`.
:::

## Example Method File

An example method file with four fitting steps is shown below:

```toml title="method.toml"
[STEP1]
INCLUDE = [15, 31, 33, 34, 37]
GRID    = [
    "[KEX_AB] = log(100.0, 600.0, 10)",
    "[PB] = log(0.03, 0.15, 10)",
    "[DW_AB] = lin(0.0, 10.0, 5)",
]

[STEP2]
FIT = ["PB", "KEX_AB", "DW_AB"]
STATISTICS = { "MC"=100, "BS"=100, "BSN"=100 }

[STEP3]
INCLUDE = "ALL"
FIX     = ["PB", "KEX_AB"]
GRID    = ["[DW_AB] = lin(0.0, 10.0, 20)"]

[STEP4]
FIT = ["DW_AB"]
```

This method file performs the following steps:

1. A subset of profiles is selected, and a grid search is performed on `"KEX_AB"`, `"PB"`, and `"DW_AB"`.
2. The parameters `"KEX_AB"`, `"PB"`, and `"DW_AB"` are varied, using the same profile selection as in Step 1.
3. All profiles are included, `"KEX_AB"` and `"PB"` are fixed, and a grid search is performed on `"DW_AB"`.
4. `"DW_AB"` is varied.

Results are saved in directories named according to each step.

## Explicit Step Semantics (Version 2)

Set `FORMAT_VERSION = 2` at the top of a method file to use canonical,
step-local semantics. Each step starts from the sealed parameter model's roles
unless `ROLES_FROM` names an earlier step. Starting parameter values always come
from the latest accepted fit; inheriting roles is not required to carry fitted
values into a later step.

`ROLES` is an ordered list. Each entry contains exactly one of `FIT`, `FIX`, or
`CONSTRAIN`. A later matching action replaces the parameter's complete earlier
role, including removal of an earlier constraint.

```toml title="method-v2.toml"
FORMAT_VERSION = 2

[FIRST]
INCLUDE = [15, 31, 33, 34, 37]
ROLES = [
    { FIX = ["PB", "KEX_AB"] },
    { FIT = ["DW_AB"] },
]

[FIRST.SEARCH.GRID]
AXES = ["[DW_AB] = lin(0.0, 10.0, 20)"]

[SECOND]
INCLUDE = "ALL"
ROLES_FROM = "FIRST"
ROLES = [
    { FIT = ["PB", "KEX_AB"] },
    { CONSTRAIN = ["[R2_B, NUC->N] = [R2_A, NUC->N]"] },
]

[SECOND.STATISTICS.MC]
REPLICATES = 100
SEED = 1234
```

Omitting `ROLES_FROM` means no role inheritance, even when the preceding step
changed roles. Omitting `INCLUDE` and `EXCLUDE` selects all profiles for that
step; selection never inherits in version 2. `SEARCH.GRID.AXES` accepts strict
`lin(...)`, `log(...)`, and `values(...)` expressions; the legacy bare tuple
form is not accepted.

### Selected-coordinate DE followed by full TRF

Use `SEARCH.DE` when a small number of global or basin-defining coordinates
need broad stochastic coverage and a Cartesian GRID would require too many
complete local fits:

```toml title="method-de-v2.toml"
FORMAT_VERSION = 2

[STEP]
ROLES = [
    { FIT = ["PB", "KEX_AB", "DW_AB"] },
]

[STEP.SEARCH.DE]
SEED = 597
COORDINATES = [
    "[PB] = log(0.001, 0.1)",
    "[KEX_AB] = log(100.0, 5000.0)",
]
```

DE searches only `PB` and `KEX_AB`. During that search, every other independent
coordinate remains at the value captured at the start of the method step. The
best eligible candidate then initializes one fresh normal TRF fit that releases
the complete final `FIT` set, including every fitted `DW_AB`. Only that final
TRF result can be accepted, committed, or used for deterministic uncertainty
and requested statistics. The current committed value may lie outside its
declared DE range; ChemEx preserves that captured value and initializes the DE
backend inside the declared search box.

Every DE coordinate must resolve to exactly one final independent `FIT`
parameter. Its `lin(low, high)` or `log(low, high)` range must be finite,
strictly ordered, and inside the parameter's physical bounds; logarithmic lower
bounds must be positive. `SEED` is a required unsigned 64-bit integer. DE
population and stopping details are ChemEx policy rather than method-file
controls.

Use GRID for deterministic Cartesian multistart. Use selected-coordinate DE
when Cartesian coverage of a few basin coordinates is impractical. A complete
three-state DCEST example is shipped as
`examples/Experiments/DCEST_15N_3States/Methods/method_de_v2.toml`.

## Setting Parameter Behavior

Parameters can be set to vary, be fixed, or be constrained. You can adjust parameter behavior in each fitting step.

### `FIT`

Parameters in the `FIT` list are varied during the fitting step. Parameters may be specified individually or as groups.

Example:

```toml
FIT = [
    "R2_A, NUC->G23N, B0->800.13MHz, T->23C",
    "R1_A, B0->800.13MHz",
    "DW_AB, NUC->N",
    "R2_B",
]
```

-   `"R2_A, NUC->G23N, B0->800.13MHz, T->23C"` specifies the R<sub>2</sub> of the nitrogen nucleus in state A of Gly23 at 800.13 MHz and 23 ºC.
-   `"R1_A, B0->800.13MHz"` includes all R₁ values for state A at 800.13 MHz.
-   `"DW_AB, NUC->N"` applies to nitrogen chemical shift differences between states A and B.
-   `"R2_B"` includes all R<sub>2</sub> values for state B.

### `FIX`

Parameters in the `FIX` list remain constant during the fitting step. The format matches that of the `FIT` list.

```toml
FIX = [
    "R2_A, NUC->G23N, B0->800.13MHz, T->23C",
    "R1_A, B0->800.13MHz",
    "DW_AB, NUC->N",
    "R2_B",
]
```

### `CONSTRAINTS`

The `CONSTRAINTS` list defines constraints on parameters using mathematical expressions of other parameters.

```toml
CONSTRAINTS = [
    "[R1_B] = 0.5 * [R1_A]",
    "[R2_B, NUC->N] = [R2_A, NUC->N]",
]
```

:::info Parameter Precedence
The settings in method files follow this order of precedence: `CONSTRAINTS` -> `FIX` -> `FIT`. This means:

-   Parameters in `CONSTRAINTS` are initialized first.
-   Parameters listed in `FIX` will be held constant, even if also present in `CONSTRAINTS`.
-   Parameters in `FIT` are set to vary, taking precedence over both `FIX` and `CONSTRAINTS` settings.

For instance, if a parameter is both constrained (in `CONSTRAINTS`) and fixed (in `FIX`), it will remain fixed, ignoring the constraint. Similarly, if a parameter is in both `FIX` and `FIT`, it will ultimately be set to vary.
:::

## Selecting a Subset of Profiles

### `INCLUDE`

The `INCLUDE` key specifies residues for analysis in each fitting step. Residues can be specified by spin-system name (e.g., `"G23N-H"`), group name (e.g., `"G23"`), or residue number (e.g., `23`). The default value, `"ALL"` (or `"*"`), includes all residues.

:::note
When using residue numbers only, provide a list of integers without quotes:

```toml
INCLUDE = ["G2", "A4", "C5", "H6"]
INCLUDE = [2, 4, 5, 6]
```

:::

### `EXCLUDE`

The `EXCLUDE` key excludes specific residues from analysis, using the same format as `INCLUDE`.

## Running a Grid Search

ChemEx supports n-dimensional grid search with 1D and 2D plot outputs for visualizing χ² values. Define grid search with the `GRID` key in any fitting step.

Grid points are defined as:

-   Linear scale:
    ```toml
    "[PB] = lin(<min>, <max>, <nb of points>)"
    ```
-   Log scale:
    ```toml
    "[PB] = log(<min>, <max>, <nb of points>)"
    ```
-   Specific points:
    ```toml
    "[PB] = (<value1>, <value2>, ..., <valuen>)"
    ```

Example:

```toml
GRID = [
    "[PB] = log(0.03, 0.1, 10)",
    "[KEX_AB] = log(200.0, 1000.0, 10)",
    "[DW_AB] = lin(0.0, 10.0, 10)",
]
```

Results are stored in the `grid.toml` file, with 1D and 2D plots generated. If more than two parameters are defined, sub-grids of interdependent parameters are evaluated, yielding 3D grids.

## Estimating Parameter Uncertainty

ChemEx offers additional methods for estimating uncertainty, including Monte Carlo, bootstrap, nucleus-specific bootstrap, and MCMC analyses.

MC, BS, and BSN start from the accepted native Direct TRF fit. Their refits are
evaluation-only successors: they cannot update the committed central parameter
values or populate the generic fitted-parameter error field.

### Monte Carlo Simulations

In Monte Carlo simulations, the fit is run once, Gaussian noise is added to generated profiles, and fitting is repeated. After N simulations, the distribution of fitted parameters provides an uncertainty estimate.

### Bootstrap Analysis

Bootstrap analysis randomly resamples data points from each profile to create synthetic datasets of the same size, and then fits are repeated.

### Nucleus-Specific Bootstrap Analysis

In nucleus-specific bootstrap, profiles are resampled based on the associated nucleus, creating synthetic datasets. Profiles with fewer data points may yield datasets of varying sizes.

:::note
Nucleus-specific bootstrap can create datasets of varying sizes, unlike standard bootstrap analysis.
:::

### MCMC Analysis

MCMC analysis samples the posterior distribution after the native deterministic
fit has converged and committed. It consumes that accepted native problem and
evaluation engine directly, uses ChemEx's weighted residuals as the likelihood,
and uses the fitted-parameter bounds as uniform priors. Initial walkers are
placed reproducibly by applying small seeded perturbations to the committed
fitted values and clipping them into the open bounded intervals.

Every fitted parameter sampled by MCMC must have finite lower and upper bounds,
with the lower bound strictly smaller than the upper bound. ChemEx rejects an
unbounded or empty interval before sampling, writes an incomplete diagnostic,
and leaves the committed central fit unchanged.

### Syntax

Run these analyses at the end of any fitting step using the `STATISTICS` key.

```toml
[STEP1]
STATISTICS = {"MC"= 100}
```

Available types:

-   `"MC"` for Monte Carlo
-   `"BS"` for bootstrap
-   `"BSN"` for nucleus-specific bootstrap
-   `"MCMC"` for posterior sampling with MCMC

To perform multiple analyses:

```toml
[STEP1]
STATISTICS = {"MC"= 100, "BS"= 100}
```

ChemEx runs every requested MC, BS, and BSN replicate. A failed or interrupted
replicate makes that analysis incomplete: successful rows and failure
diagnostics remain available, but ChemEx does not publish a complete summary,
correlation matrix, or plots. The native root seed is recorded in
`diagnostics.toml`; serial and worker execution preserve the same ordered
results.

For MCMC, the compact form sets the number of sampler steps:

```toml
[STEP1]
STATISTICS = {"MCMC"= 5000}
```

Advanced MCMC settings can be provided as a nested table:

```toml
[STEP1.STATISTICS.MCMC]
STEPS = 5000
BURN = "AUTO"
THIN = 1
WALKERS = 64
SEED = 1234
WORKERS = 4
```

`STEPS` is the number of emcee sampler iterations. `WALKERS` defaults to the
larger of 32 or twice the number of fitted parameters. `SEED` fixes the initial
walker positions for reproducible runs.

`BURN` defaults to `"AUTO"`. In automatic mode, ChemEx uses the integrated
autocorrelation time reported by the sampler and discards twice the largest
autocorrelation time when that estimate is available and shorter than the chain.
If the chain is too short for emcee's reliability threshold but still provides a
tentative autocorrelation estimate, ChemEx uses that tentative estimate for
automatic burn-in and records the warning in `diagnostics.toml`. If
autocorrelation time cannot be estimated at all, ChemEx keeps the full chain and
records the reason. A numeric `BURN` value can still be provided to discard a
fixed number of initial sampler steps.

`THIN` defaults to `1`, which keeps every retained sample. Thinning is mainly a
storage and output-size control; it is usually better to keep all post-burn-in
samples unless the output files become too large.

`WORKERS` is optional. When omitted, MCMC inherits the command-line
[`--workers`](multicore_execution.md) setting from `chemex fit`, whose default
is `auto`. Set `WORKERS` in the method file only when this MCMC step needs a
specific positive worker count. Command-line `--workers 0` means all available
CPUs, but method-file `WORKERS` must be a positive integer.

`UPDATE_PARAMETERS = true` is rejected. Posterior summaries are evidence about
the committed central fit and never replace its authoritative parameter values
or generic fitted-parameter errors.

See [Multicore Execution](multicore_execution.md) for performance guidance and
the interaction between `--workers` and `--native-threads`.

Sampling outputs are stored under a `Statistics` directory in the corresponding step or group output directory:

```text
Statistics/
  MonteCarlo/
    summary.toml
    samples.tsv
    correlations.tsv
    diagnostics.toml
    plots.pdf
  Bootstrap/
    summary.toml
    samples.tsv
    correlations.tsv
    diagnostics.toml
    plots.pdf
  BootstrapNS/
    summary.toml
    samples.tsv
    correlations.tsv
    diagnostics.toml
    plots.pdf
  MCMC/
    summary.toml
    samples.tsv
    correlations.tsv
    diagnostics.toml
    plots.pdf
```

For Monte Carlo and bootstrap methods, `samples.tsv` contains one fitted-parameter row per synthetic dataset plus χ². `summary.toml` reports empirical-distribution quantities: mean, median, standard deviation, 95% percentile bounds, 68% percentile bounds, and `half_percentile_68_width`. `correlations.tsv` reports parameter correlations across the fitted synthetic datasets. Missing values are written as `nan`. Their `plots.pdf` reports provide a summary page, one-dimensional sample distributions, a χ² distribution, and two-dimensional sample distributions for parameter pairs with `|r| >= 0.5`.

For MCMC, `summary.toml` reports the uniform prior implied by each parameter's bounds, posterior mean, median, standard deviation, 95% equal-tailed and 68% credible-interval bounds, and `half_credible_interval_68_width`. Effective sample size and `mcse_mean` are included when the autocorrelation estimate passes emcee's reliability threshold; `mcse_mean` is the Monte Carlo standard error of the estimated posterior mean, not a posterior interval width. MCMC diagnostics include sampler versions, retained samples, acceptance fractions, reliable or tentative autocorrelation time, recommended chain lengths, and burn-in decisions. The `plots.pdf` report provides a summary page, one-dimensional posterior distributions, walker traces, the log-probability trace, an autocorrelation monitor, and two-dimensional posterior distributions for parameter pairs with `|r| >= 0.5`.

MC, BS, BSN, and MCMC summaries remain separate statistical results. They do not
mutate the committed central fit or populate a generic error in
`Parameters/fitted.toml`.

If MCMC setup, sampling, processing, or output fails or is interrupted,
`diagnostics.toml` records an explicit incomplete terminal state. ChemEx does not
publish `summary.toml`, `samples.tsv`, `correlations.tsv`, or `plots.pdf` for an
incomplete chain.
