---
sidebar_position: 6
---

# Method Files

## Overview

Method files define the fitting methods used in ChemEx. In these files, you can:

-   Specify which parameters to fit, fix, or constrain.
-   Select the profiles to include in calculations.
-   Activate additional analyses, such as grid search, Monte Carlo, or bootstrap analyses.

Provide the method file to ChemEx using the `-m` or `--method` option.

Method files are structured in sections, each representing a fitting step. Fitting steps are executed in the order they appear, and parameter settings inherit from previous steps if not redefined.

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

ChemEx offers additional methods for estimating uncertainty, including Monte Carlo, bootstrap, and nucleus-specific bootstrap analyses.

### Monte Carlo Simulations

In Monte Carlo simulations, the fit is run once, Gaussian noise is added to generated profiles, and fitting is repeated. After N simulations, the distribution of fitted parameters provides an uncertainty estimate.

### Bootstrap Analysis

Bootstrap analysis randomly resamples data points from each profile to create synthetic datasets of the same size, and then fits are repeated.

### Nucleus-Specific Bootstrap Analysis

In nucleus-specific bootstrap, profiles are resampled based on the associated nucleus, creating synthetic datasets. Profiles with fewer data points may yield datasets of varying sizes.

:::note
Nucleus-specific bootstrap can create datasets of varying sizes, unlike standard bootstrap analysis.
:::

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

To perform multiple analyses:

```toml
[STEP1]
STATISTICS = {"MC"= 100, "BS"= 100}
```

Outputs are stored in a file in the corresponding step directory. When data is unavailable (e.g., in nucleus-specific bootstrap), placeholders (`"--"`) are used.
