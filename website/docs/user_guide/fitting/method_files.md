---
sidebar_position: 6
---

# Method files

## Description

The method file contains the fitting methods to be used during the fitting
process. This is where you can:

-   define which parameters are fit, fixed, or constrained,
-   select the profiles to include in the calculation,
-   activate additional calculations, such as grid search, Monte Carlo or
    bootstrap analyses.

This file is provided to ChemEx using the `-m` or `--method` option.

The method file is structured in sections. Each section corresponds to a fitting
step. The section name is the name of the fitting step. The steps are run in the
order they are defined.

If the fitting method contains multiple fitting steps, the value and behavior of
each parameter always inherit from the previous fitting step if not set in the
current step. This means that if a parameter is fixed in one step, it remains
fixed in the following steps as long as its state is not changed.

:::tip

The set of residues to be included in global fits should be chosen carefully. A
commonly used multi-step fitting strategy is to select a subset of residues with
relatively large CPMG dispersion or good quality CEST minor dips to first get
global parameters (p<sub>B</sub>, k<sub>ex</sub>), and then carry out
single-residue fits with (p<sub>B</sub>, k<sub>ex</sub>) fixed to get
residue-specific parameters (e.g., Δϖ) in the next step. In CPMG experiments, to
get reasonable initial estimates of Δϖ for each residue, an additional
single-residue fitting step can be carried out at the very beginning, see the
method file for `CPMG_CH3_1H_SQ/` under `Examples/Experiments/` for such an
example.

:::

## Example

Here is an example method file demonstrating the different possibilities:

```toml title="method.toml"
[STEP1]
INCLUDE = [15, 31, 33, 34, 37]                                                                                              #, 33, 34, 37]
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

This fit contains 4 distinct steps:

1. A subset of profiles are selected, and a grid search is performed on the
   parameters `"KEX_AB"`, `"PB"` and `"DW_AB"`.
2. The parameters `"KEX_AB"`, `"PB"` and `"DW_AB"` are set to vary. The
   selection of profiles remains the same as in the first step.
3. All profiles are included in the fit, the parameters `"KEX_AB"` and `"PB"`
   are fixed, and a grid search is performed on the parameter `"DW_AB"`.
4. The parameter `"DW_AB"` is set to vary.

Results are written on separate directories named according to the corresponding
section name.

## Setting parameter behavior

All model parameters can either be varied, fixed or constrained. The status of
each parameter is defined by default at the beginning of the fitting process.
However, you can change this behavior in each fitting step of the method file.

### `FIT`

Parameters in the `FIT` list are varied during the fitting process. Here,
parameter names can either designate a unique parameter or a group of
parameters. For the later, simply mention the attributes identifying the group.

For example:

```toml
FIT = [
    "R2_A, NUC->G23N, B0->800.13MHz, T->23C",
    "R1_A, B0->800.13MHz",
    "DW_AB, NUC->N",
    "R2_B",
]
```

In this example:

-   `"R2_A, NUC->G23N, B0->800.13MHz, T->23C"` corresponds to the amide nitrogen
    R<sub>2</sub> of state A of Gly23 measured at 23 ºC, 800.13 MHz.
-   `"R1_A, B0->800.13MHz"` selects all the state A R<sub>1</sub> values measured
    at 800.13 MHz, independently of the residue number and temperature.
-   `"DW_AB, NUC->N"` corresponds to all the amide nitrogen chemical shift
    differences between states A and B .
-   `"R2_B"` selects all R<sub>2</sub> values os state B.

### `FIX`

Parameters in the `FIX` list are fixed during the fitting process. The format is
similar to the `FIT` list.

```toml
FIX = [
    "R2_A, NUC->G23N, B0->800.13MHz, T->23C",
    "R1_A, B0->800.13MHz",
    "DW_AB, NUC->N",
    "R2_B",
]
```

### `CONSTRAINTS`

The `CONSTRAINTS` list defines the constraints to be applied to the parameters.
Constraints are mathematical expression of other parameters. The value of the
constrained parameter is calculated using this expression.

Parameters in the mathematical expression given in the `CONSTRAINTS` list should
be put in brackets.

```toml
CONSTRAINTS = [
    "[R1_B] = 0.5 * [R1_A]",
    "[R2_B, NUC->N] = [R2_A, NUC->N]",
]
```

:::info

Keys are read in that order: `CONSTRAINTS` -> `FIX` -> `FIT`.

:::

## Selecting a subset of profiles

### `INCLUDE`

The `INCLUDE` key in method file allows selecting a subset of residues for
analysis during each fitting step. The residue name should match the spin-system
assignment provided in the experiment file(s). You can use the full spin-system
name (e.g. **"G23N-H"**) or the group name (e.g. **"G23"**) or the residue
number (e.g. **23**). `"ALL"` (or `"*"`) is the default value, which indicates
that all residues are to be included in the current fitting step.

:::note

When only the residue number is used, use a list of integer, that is omit the
quotes. These two formulations are equivalent:

```toml
INCLUDE = ["G2", "A4", "C5", "H6"]
INCLUDE = [2, 4, 5, 6]
```

:::

### `EXCLUDE`

The `EXCLUDE` key in method file allows excluding a subset of residues from
analysis during each fitting step. its usage is similar to the `INCLUDE` key.

## Running a grid search

ChemEx has a built-in grid search method that offers the possibility to run an
nD grid search and plot the resulting χ<sup>2</sup> values as 1D and 2D plots.
Grid search can be defined and run using the key `GRID` in any section of the
method file.

The grid is defined on a parameters basis. The parameters defining the nD grid
are fixed to the value of the grid that is evaluated, while the other parameters
are set as defined by the `FIX`, `FIT` and `CONSTRAINTS` options. Points of the
grid can be defined using a linear scale, a log scale or point by point:

-   Linear scale:

    ```toml
    "[PB] = lin(<min>, <max>, <nb of points>)"
    ```

-   Log scale:

    ```toml
    "[PB] = log(<min>, <max>, <nb of points>)"
    ```

-   Point by point:

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

Parameters are selected as usual. For example, `[R2_A]` would select all the
**R2_A**. `[R2_A, NUC->G43N]` would select all the **R2_A** parameters of the
nucleus **G43N** (multiple values if there are several temperatures, B0 field).

At the end of the grid search the best point is selected and the corresponding
parameters are used in the next step of the fitting procedure. The χ<sup>2</sup>
values are reported in the `grid.toml` file as well as in the form of 1D and 2D
plots.

For an nD grid with n > 1, multiple 1D grids are plotted, one per parameter,
corresponding to the parameter values on the x axis and the minimum
χ<sup>2</sup> value along the other dimension.

Similarly, if n > 2, then a series of 2D χ<sup>2</sup> surface corresponding to
each pair of independent parameters is plotted. 2D surface correspond to 2D
projection in which the best fit values in the other dimensions are used to
evaluate each point.

:::info

When two parameters are independent from each others, the corresponding
χ<sup>2</sup> surface is entirely flat. It then becomes possible to define
sub-grids of parameters that are all dependent on each others. The algorithm
used in ChemEx starts by defining these minimal individual grids and then
evaluates them separately. The evaluated sub-grids may share common parameters.
They are therefore combined together at the end to recover the global minimum
and to plot 1D χ<sup>2</sup> curve for each parameter and 2D plots when
possible.

This algorithm is much faster in the sense that a search involving 2 global
parameters and 8 independent residue specific parameters does not generate a
10-dimensional grid, but 10 3-dimensional grids, which is much faster and
importantly retains all the information.

:::

## Evaluating the uncertainty on the fitted parameters

The uncertainty on fitted parameters is, in general, estimated through the
covariance matrix obtained from the Levenberg-Marquardt optimization. However,
ChemEx offers additional methods to evaluate the parameter uncertainties, that
is, Monte Carlo simulation, bootstrap analysis and nucleus-specific bootstrap.

### Monte Carlo simulations

For the Monte Carlo simulation, the fit is run once and Gaussian noise is added
to the back-calculated values based on the error. Fits are subsequently run on
these generated profiles. After N simulations, the distribution of the fitted
parameters provides an estimate of the uncertainty on the fitted parameters.

### Bootstrap analysis

The bootstrap analysis is similar to the Monte Carlo simulations, except that
the synthetic profiles are realized by randomly picking data points from each
profile to generate new ones with the same number of points as the original.

### Nucleus-specific bootstrap analysis

For the nucleus-specific bootstrap analysis, full profiles are randomly selected
based on their associated nucleus to generate the synthetic datasets. In other
words, if we have a dataset that depends on the nuclei \{G2N, H8N, R9N, R9H\},
potential new datasets could include the profiles of the following sets of
nuclei \{H8N, H8N, R9N, R9H\} or \{G2N, G2N, G2N, R9H\}.

:::note

Contrary to standard bootstrap analysis, nucleus-specific bootstrap analysis can
produce datasets with different number of data points in them. For examples, if
profiles depending on **H8N** appears in multiple experiments and the ones
depending on **G2N** in only one, then the two datasets mentioned above would
have different number of data points. This goes against the main principle
underlying bootstrap analysis that normally requires that all the newly sampled
datasets are of the same size.

:::

### Syntax

These calculations can be run at the end of any fitting step by using the key
`STATISTICS`.

The syntax is the following:

```toml
    [STEP1]
    STATISTICS = {"MC"= 100}
```

where MC is the type of simulation and 100 is the number of simulations.

Types can be:

-   "MC" for Monte Carlo
-   "BS" for bootstrap
-   "BSN" for nucleus-specific bootstrap

To run two or more types of simulation just add additional pairs of values:

```toml
[STEP1]
STATISTICS = {"MC"= 100, "BS"= 100}
```

The output for each kind of simulation is stored in a single file stored in the
directory corresponding to the step it belongs to. Parameter values are stored
in different columns. When no values are available, which can be the case for
nucleus-specific bootstrap analysis, the characters "--" are used to fill the
space.
