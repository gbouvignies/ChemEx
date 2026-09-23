---
sidebar_position: 2.75
sidebar_label: Fit ¹⁵N CEST data
title: Fit ¹⁵N CEST Data
description: Fit real two-B1 15N CEST profiles with a two-state exchange model, inspect A/B-labeled chemical shifts, and compare the calculated profiles.
---

# Fit ¹⁵N CEST data

This tutorial continues the practical path from
[the first CPMG fit](./first_analysis.md). CPMG detects exchange through
relaxation dispersion as the pulse-train frequency changes. CEST instead
applies a weak radio-frequency field at a series of offsets and measures how
that irradiation changes the detected signal.

When the exchange regime and experiment are suitable, irradiation near another
exchanging state's resonance transfers saturation into the observed signal.
The resulting CEST profile can therefore report both on exchange and on the
chemical-shift positions assigned to the exchanging states. That extra
information is useful in this example; it is not guaranteed for every system or
every CEST dataset.

## What you will fit

You will fit the complete shipped `CEST_15N` example:

- real ¹⁵N CEST data for the A39G FF domain referenced on the
  [pure in-phase ¹⁵N CEST experiment page](./experiments/cest/cest_15n.md);
- 16 residues, each measured with nominal B1 fields of 13.0 and 26.3 Hz;
- 32 profiles fitted simultaneously with the two-state `2st` model;
- one population and exchange rate shared by every profile, together with
  residue-specific chemical-shift and relaxation information.

ChemEx uses stable A/B state labels throughout this analysis. The shipped
parameters assign the directly observed resonance to A and the other
exchange-related position to B. Those labels are an analyst-defined description
of this dataset, not a population or energy ordering imposed by ChemEx. See
[Exchange States and Parameters](./user_guide/fitting/exchange_states_parameters.md)
for the canonical conventions.

## Open the example

If you completed the CPMG tutorial from a repository checkout, reuse that same
checkout; there is no need to clone ChemEx again. From its root, change into:

```shell
cd examples/Experiments/CEST_15N
```

If you do not yet have the examples, follow
[Obtain the example](./first_analysis.md#obtain-the-example), then use the
directory above. A PyPI installation provides the application but does not
install the repository's example directory.

The files needed here are:

```text
CEST_15N/
├── Data/
│   ├── 13Hz/                 # 30 raw profile files; 16 mapped below
│   └── 26Hz/                 # 30 raw profile files; 16 mapped below
├── Experiments/
│   ├── 13hz.toml
│   └── 26hz.toml
├── Parameters/
│   └── parameters.toml
├── pick_cest.sh              # optional interactive initialization
├── run.sh                    # complete shipped fit
└── simulate.sh               # optional simulation
```

Each data directory contains 30 raw profile files; the corresponding experiment
TOML explicitly maps the 16 profiles used in this tutorial. ChemEx does not
automatically load every file merely because it is present in the directory.

The generic roles of experiment, parameter, and data files are already covered
by the first tutorial. Here, focus on what is specific to CEST.

## Read one CEST profile

Open `Data/13Hz/52N-HN.out`. Its three columns are the irradiation offset from
the configured carrier in Hz, the measured intensity, and a file uncertainty:

```text title="Data/13Hz/52N-HN.out"
#Offset (Hz)        Intensity    Uncertainty
  -1.000e+05    5.0370030e+04  4.5534923e+02
  -7.500e+02    2.1683766e+04  2.3381834e+02
  -7.250e+02    2.1929222e+04  2.3550312e+02
```

The first row is a reference measurement made far outside the CEST sweep. The
remaining rows scan the weak irradiation field across offsets from -750 to
+750 Hz.

ChemEx's generated CEST plots convert those offsets to the absolute ¹⁵N
irradiation position in ppm. The x-axis is labeled **B1 position (ppm)** and is
displayed in the usual decreasing NMR direction. The y-axis is `I/I0`: the
measured or calculated profile intensity divided by the mean experimental
reference intensity for that profile. Reference rows define `I0` but are not
drawn in the PDF profile.

The points with error bars are measurements, and the red curve is calculated
from the fitted model. A pronounced dip occurs when irradiation directly
affects the observed A-labeled resonance. In this dataset, a second, shallower
feature occurs near the B-labeled position because saturation is transferred
through exchange.

## Why there are two B1 fields

`b1_frq` is the nominal amplitude of the weak irradiation field, expressed as a
frequency in Hz. It is not the offset being scanned: each data row supplies an
offset, while `b1_frq` controls the irradiation strength used at every point in
that experiment file.

Changing B1 changes the widths and depths of the features in the measured
profile. The 13.0 and 26.3 Hz profiles therefore make different observations of
the same exchange process. ChemEx fits them with the same population, exchange
rate, and residue-specific shifts and relaxation parameters. A parameter set
must explain both responses at once, which can constrain it more strongly than
either profile alone, although it does not guarantee that every parameter is
uniquely determined.

## Inspect the experiment files

Open `Experiments/13hz.toml`. Its CEST-specific settings are:

```toml title="Experiments/13hz.toml"
[experiment]
name = "cest_15n"
time_t1 = 0.5
carrier = 118.987
b1_frq = 13.0

b1_distribution = { type = "dephasing" }

[conditions]
h_larmor_frq = 499.243

[data]
path = "../Data/13Hz/"
error = "scatter"
filter_offsets = [[0.0, 13.0]]
```

- `time_t1 = 0.5` sets the CEST saturation period to 0.5 seconds.
- `carrier = 118.987` places the ¹⁵N carrier at 118.987 ppm. The offsets in each
  data file are in Hz relative to this carrier.
- `b1_frq = 13.0` sets the nominal B1 field to 13.0 Hz.
- `b1_distribution = { type = "dephasing" }` selects ChemEx's complete-dephasing
  treatment for this example's B1 behavior; it is a named calculation mode, not
  a numerical distribution that needs additional parameters.
- `h_larmor_frq` is the ¹H Larmor frequency in MHz and supplies the field needed
  for frequency conversion and the spin calculation.

`Experiments/26hz.toml` has the same carrier, saturation time, magnetic field,
and 16 profile names. It changes `b1_frq` to 26.3 Hz, points to `Data/26Hz/`, and
uses `filter_offsets = [[0.0, 26.0]]`.

See the [experiment-file reference](./user_guide/fitting/experiment_files.md)
for the common tables and the
[¹⁵N CEST experiment reference](./experiments/cest/cest_15n.md) for all
supported settings.

## Inspect the starting parameters

The most relevant parts of `Parameters/parameters.toml` are:

```toml title="Parameters/parameters.toml"
[GLOBAL]
PB = 0.015
KEX_AB = 70.0

[CS_A]
52N = 111.358
55N = 128.301

[DW_AB]
52N = 8.5
55N = -6.5
```

The global `PB` and `KEX_AB` values initialize the shared population and
exchange rate. `CS_A` gives the fixed A-labeled chemical-shift position for each
residue. Each `DW_AB` initializes that residue's B-minus-A shift difference, so
its sign determines which side of A the initial B position occupies. The full
file supplies all 16 residues and a `TAUC_A` value used to initialize relaxation
parameters.

These numbers are starting information, not fitted truth. The current parameter
file was prepared for this dataset, so the command-line fit is reproducible
without an interactive picking step. For exact definitions and label semantics,
use [Exchange States and Parameters](./user_guide/fitting/exchange_states_parameters.md),
and see the [parameter-file reference](./user_guide/fitting/parameter_files.md)
for initialization and bounds.

## Filtering, references, scaling, and uncertainties

The two files deliberately filter a narrow band centered at zero Hz relative to
the A-labeled resonance: 13 Hz wide in `13hz.toml` and 26 Hz wide in
`26hz.toml`. On the sampled grids, that excludes 9 observations in the 13.0 Hz
dataset and 16 observations in the 26.3 Hz dataset. They do not contribute to
profile scaling or the fit objective. They remain in `Data/13hz.dat` and
`Data/26hz.dat` marked `NOT USED IN THE FIT`. They also remain in the `.exp`
plot-data files and are drawn in pale red in the PDFs.

Each raw profile also contains one far-off-resonance reference row at
`-100000 Hz`. This example leaves the default `filter_ref_planes = false`, so
those 32 reference observations remain active in profile scaling and fitting.
They remain in the machine-readable `Data/*.dat` files but are omitted from the
CEST PDFs and `Plots/*.exp` files, whose x-axis covers only the saturation
sweep.

The example also uses the default scaled-profile behavior and `error =
"scatter"`; default experiment-level pooling applies separately within each of
the two experiment files. The formulas, mask ordering, uncertainty semantics,
fit-residual convention, and fit-statistic counts are defined in the canonical
[Scaling, Uncertainties, Residuals, and Fit Statistics](./user_guide/fitting/scaling_uncertainties_residuals.md)
page.

## Run the fit

From `ChemEx/examples/Experiments/CEST_15N`, run:

```shell
chemex fit -e Experiments/13hz.toml Experiments/26hz.toml -p Parameters/parameters.toml -d 2st -o OutputTutorial
```

This is the shipped workflow with explicit experiment filenames and an explicit
model name. It performs one direct bounded trust-region-reflective fit of all 32
profiles. No Method file, grid search, stochastic search, resampling, or MCMC is
used. Method files remain useful when an analysis needs selection, staged roles,
or additional searches; this example simply does not need one for its direct
fit.

## What is being fitted

The current run resolves 66 optimizer-controlled parameters:

- 2 parameters shared by all profiles: `PB` and `KEX_AB`;
- 16 residue-specific `DW_AB` values;
- 16 residue-specific `R1_A`, 16 `R2_A`, and 16 `R2_B` values.

Both experiment files have the same field and conditions, so a residue's
relaxation parameters are shared between its 13.0 and 26.3 Hz profiles. The 16
`CS_A` values remain fixed. `PA`, `KAB`, `KBA`, `CS_B`, and `R1_B` are derived
from the independent parameters rather than varied separately.

In addition, ChemEx analytically profiles one amplitude scale for each of the 32
profiles. Those scales are not part of the 66 optimizer-controlled coordinates.

## Confirm success

The command returns to the shell after reporting `Making plots...` and listing
`13hz.pdf` and `26hz.pdf`. The authoritative completion record is
`OutputTutorial/run_info/outcome.toml`:

```toml
schema_version = 2
status = "complete"
```

The key outputs are:

```text
OutputTutorial/
├── run_info/outcome.toml
├── Parameters/
│   ├── fitted.toml
│   ├── fixed.toml
│   └── constrained.toml
├── Data/
│   ├── 13hz.dat
│   └── 26hz.dat
├── Plots/
│   ├── 13hz.pdf
│   ├── 13hz.exp
│   ├── 13hz.fit
│   ├── 26hz.pdf
│   ├── 26hz.exp
│   └── 26hz.fit
└── statistics.toml
```

A current run reports covariance available with a full-rank 66-parameter
Jacobian and no boundary warnings. It gives a reduced χ² of about 0.762 under
the configured scatter-uncertainty model. These diagnostics supplement, rather
than redefine, `status = "complete"`; follow the
[output reference](./user_guide/fitting/outputs.mdx) and canonical statistics
page when assessing a scientific fit.

The run also reports `"Kolmogorov-Smirnov test" = 9.02501e-06`. ChemEx computes
this value by comparing the normalized residual distribution with a standard
normal distribution. Because the residuals come from a fitted model and are not
necessarily independent, this p-value is best treated as a diagnostic rather
than as a formal model-acceptance or rejection test. Small values can highlight
departures from the idealized residual distribution and are worth considering
alongside the residual plots, uncertainty model, and other fit diagnostics;
they do not imply that the optimization failed or invalidate an otherwise
complete fit.

## Inspect the exchange and shift results

Open `OutputTutorial/Parameters/fitted.toml`. Current results are close to:

| Quantity       |   Approximate result | Output role                        |
| -------------- | -------------------: | ---------------------------------- |
| `PB`           | 0.01618 (about 1.6%) | fitted, shared across all profiles |
| `KEX_AB`       |            58.29 s⁻¹ | fitted, shared across all profiles |
| `DW_AB`, `52N` |           +8.410 ppm | fitted for residue 52N             |
| `DW_AB`, `55N` |           -6.401 ppm | fitted for residue 55N             |

The fixed and derived shift information is split across output files:

| Residue | Fixed `CS_A` | Fitted `DW_AB` | Derived `CS_B` |
| ------- | -----------: | -------------: | -------------: |
| 52N     |  111.358 ppm |     +8.410 ppm |    119.768 ppm |
| 55N     |  128.301 ppm |     -6.401 ppm |    121.900 ppm |

Read the A positions in `Parameters/fixed.toml`, the fitted differences in
`Parameters/fitted.toml`, and the resulting B positions in
`Parameters/constrained.toml`. This separation is useful: it shows which
chemical-shift information was supplied, optimized, or derived. Do not infer a
major/minor or ground/excited ordering from the letters alone.

## Compare the 52N profiles

Open the page titled `52N` in both `Plots/13hz.pdf` and `Plots/26hz.pdf`.

In each lower panel:

- red points with error bars are retained measurements;
- pale-red points, when present, were excluded by `filter_offsets`;
- the red curve is the fitted calculation;
- the solid vertical line is the A-labeled position;
- the next, densely dashed vertical line is the B-labeled position.

The line styles follow ChemEx label order; they do not encode physical identity,
population, or energy. The small upper panel is a visual difference diagnostic.
Use the canonical scaling and residuals page for the residual vector actually
minimized by the fit.

For 52N, the deep feature near 111.36 ppm aligns with `CS_A`. The shallower
exchange-related feature near 119.77 ppm aligns with the derived `CS_B`. In the
current calculated curves, the B-labeled feature reaches about `I/I0 = 0.33` at
13.0 Hz and `I/I0 = 0.29` at 26.3 Hz, while the off-resonance level is about
0.44. The 26.3 Hz features are also visibly broader. These are observations of
this fitted pair of profiles, not universal rules for how every CEST feature
must change with B1.

This side-by-side comparison shows why the joint fit is informative: the same
`PB`, `KEX_AB`, and 52N shift and relaxation parameters must reproduce two
distinct B1 responses.

## Optional: initialize shifts with `pick_cest`

The reproducible fit above does not require the GUI because the shipped
parameter file already contains suitable starting shifts. When preparing a new
CEST analysis, the existing `pick_cest.sh` provides an optional initialization
workflow:

```shell
chemex pick_cest -e Experiments/13hz.toml -o Sandbox
```

The window lets you inspect each CEST profile, click the position to assign to A
and then the position to assign to B, and swap or clear those choices. ChemEx
writes `Sandbox/cs_a.toml` and `Sandbox/dw_ab.toml` as you work.

`pick_cest` helps assign and initialize A/B-labeled chemical shifts; those
labels remain modeling choices, and ChemEx does not subsequently reorder them
by fitted population. The tool does not prove which physical conformation is
major, minor, ground, or excited. See the
[`pick_cest` reference](./user_guide/additional_modules.mdx#pick-cest) and the
canonical exchange-state page before using the generated files.

## Try simulation next

`simulate.sh` uses the same two experiment files and parameter file with
`chemex simulate`. It is a useful next experiment if you want to see profiles
calculated from supplied parameter values, but simulation is not required for
this fit.

## What ChemEx just did

```text
32 experimental ¹⁵N CEST profiles at two B1 fields
                         ↓
                 two-state model
                         ↓
shared PB and KEX_AB + residue-specific shifts and relaxation
                         ↓
profile scaling and uncertainty-weighted global optimization
                         ↓
calculated CEST profiles, fitted parameters, diagnostics, and provenance
```

## Connect CEST to CPMG

| CPMG                                                | CEST                                                                                          |
| --------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| Measures exchange-sensitive relaxation dispersion   | Measures response across saturation offsets                                                   |
| Varies the refocusing pulse-train frequency         | Scans saturation offset and can compare B1 fields                                             |
| Can constrain kinetics in suitable exchange regimes | Can reveal state-specific shift positions and constrain kinetics in suitable exchange regimes |

Neither experiment guarantees a unique model or parameter set. Their
information can be complementary for some systems, but the useful combination
depends on the exchange regime, pulse sequence, field conditions, and data
quality.

## Next steps

- Revisit [Exchange States and Parameters](./user_guide/fitting/exchange_states_parameters.md).
- Review [Scaling, Uncertainties, Residuals, and Fit Statistics](./user_guide/fitting/scaling_uncertainties_residuals.md).
- Configure your own [experiment files](./user_guide/fitting/experiment_files.md)
  and [parameter files](./user_guide/fitting/parameter_files.md).
- Read the complete [pure in-phase ¹⁵N CEST experiment reference](./experiments/cest/cest_15n.md).
- Learn the full [`pick_cest` workflow](./user_guide/additional_modules.mdx#pick-cest).
- Interpret the generated [outputs and evidence](./user_guide/fitting/outputs.mdx).
- Browse the [representative examples](./examples/index.mdx).
