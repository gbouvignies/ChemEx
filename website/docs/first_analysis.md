---
sidebar_position: 2
sidebar_label: Run your first analysis
title: Run Your First ChemEx Analysis
description: Fit real two-field 15N CPMG relaxation-dispersion data with a two-state chemical-exchange model and inspect the fitted exchange parameters and profiles.
---

# Run your first ChemEx analysis

In this tutorial, you will fit real ¹⁵N CPMG relaxation-dispersion data at two
magnetic fields with a two-state exchange model. The result will contain shared
exchange parameters and fitted dispersion profiles for five residues.

The goal is a first successful, scientifically meaningful fit. The full shipped
example contains a more robust staged analysis, but you do not need that extra
machinery yet.

## Prerequisite

[Install ChemEx](./welcome_to_chemex.md#installation), then check that the
command is available:

```shell
chemex --version
```

## Obtain the example

A PyPI installation contains the application, not the repository's example
directory. Clone the ChemEx repository to obtain the example files:

```shell
git clone --depth 1 https://github.com/gbouvignies/ChemEx.git
cd ChemEx/examples/Experiments/CPMG_15N_IP
```

These commands work in macOS and Linux shells and in Windows PowerShell when
Git is installed. If you do not use Git, download and extract the
[repository ZIP](https://github.com/gbouvignies/ChemEx/archive/refs/heads/main.zip),
then open `ChemEx-main/examples/Experiments/CPMG_15N_IP` in a terminal.

## The analysis directory

The checkout contains additional files for the full example. These are the
inputs used by this tutorial:

```text
CPMG_15N_IP/
├── Data/
│   ├── 500MHz/             # one measured CPMG profile per .out file
│   └── 800MHz/
├── Experiments/
│   ├── 500mhz.toml         # 500 MHz experiment and data locations
│   └── 800mhz.toml         # 800 MHz experiment and data locations
├── Parameters/
│   └── parameters.toml     # starting parameter values
└── Methods/
    └── tutorial.toml       # one fit using five informative profiles
```

- **Experiment files** describe how the NMR measurements were performed and
  map spin-system names to data files.
- **Data files** contain measured CPMG intensities as a function of the number
  of CPMG cycles. ChemEx converts these to effective transverse relaxation
  rates for plotting.
- The **parameter file** supplies starting values for the exchange model and
  the measured state-A chemical shifts.
- The **method file** selects profiles 15, 31, 33, 34, and 37 at both fields for
  one direct fit. These are the same informative profiles used to establish the
  shared exchange parameters in the full shipped analysis.

The [fitting reference](./user_guide/fitting/index.mdx) describes these file
types in detail.

## The two-state model

The model calls the exchanging states A and B. These are stable labels; deciding
which physical conformations they represent is the analyst's responsibility.
For this fit, the most important parameters are:

- `PB`: the equilibrium population fraction of state B;
- `KEX_AB`: the total A↔B exchange rate, in s⁻¹;
- `DW_AB`: the profile-specific ¹⁵N chemical-shift difference
  (state B minus state A), in ppm.

`PB` and `KEX_AB` are shared by all ten profiles: five residues measured at
500 and 800 MHz. Each residue has its own `DW_AB`, and each residue and field
has its own baseline transverse relaxation rate. Fitting the profiles together
lets both fields and all five residues constrain the shared exchange process.

## Run the fit

From `ChemEx/examples/Experiments/CPMG_15N_IP`, copy and run this single
command:

```shell
chemex fit -e Experiments/500mhz.toml Experiments/800mhz.toml -p Parameters/parameters.toml -m Methods/tutorial.toml -d 2st -o OutputTutorial
```

Here, `-e` names the experiment files, `-p` the parameter file, `-m` the Method
file, `-d 2st` selects the two-state kinetic model, and `-o` names the output
directory.

ChemEx reports four main stages in the terminal: loading both datasets, reading
the method and starting parameters, running the minimization, and writing the
result files and plots. This tutorial method performs one bounded
trust-region-reflective fit; it does not use a grid search, stochastic search,
or resampling.

## Confirm success

The run has succeeded when the shell prompt returns without an error after
ChemEx reports `Making plots...` and lists `500mhz.pdf` and `800mhz.pdf`.
Confirm the machine-readable outcome by opening
`OutputTutorial/run_info/outcome.toml`:

```toml
schema_version = 2
status = "complete"
```

You should also have this single-step result layout:

```text
OutputTutorial/
├── run_info/
│   └── outcome.toml
├── Parameters/
│   ├── fitted.toml
│   ├── fixed.toml
│   └── constrained.toml
├── Data/
│   ├── 500mhz.dat
│   └── 800mhz.dat
├── Plots/
│   ├── 500mhz.pdf
│   └── 800mhz.pdf
└── statistics.toml
```

`status = "complete"` is the authoritative indication that every requested
fit step and output operation finished. See [Outputs](./user_guide/fitting/outputs.mdx)
for the complete result layout and provenance fields.

## Inspect the fitted parameters

Open `OutputTutorial/Parameters/fitted.toml`. A current run gives values close
to the following; the last digits can vary slightly across platforms:

| Parameter         |  Approximate result | Meaning                                        |
| ----------------- | ------------------: | ---------------------------------------------- |
| `PB`              | 0.0703 (about 7.0%) | Equilibrium population of state B              |
| `KEX_AB`          |             382 s⁻¹ | Total rate of exchange between A and B         |
| `DW_AB` for `15N` |            2.00 ppm | B-minus-A ¹⁵N shift difference for profile 15N |
| `DW_AB` for `34N` |            3.63 ppm | B-minus-A ¹⁵N shift difference for profile 34N |

The population and exchange rate are global results. The chemical-shift
differences vary by residue because each nucleus reports on its own local
environment. ChemEx also fits a baseline `R2_A` for each residue at each field;
those values let the shared exchange model account for the non-exchange part of
each dispersion profile.

## Inspect a fitted dispersion profile

Open `OutputTutorial/Plots/800mhz.pdf` and find the page titled `34N`. In the
lower panel, the points with error bars are experimental effective transverse
relaxation rates and the smooth line is the profile calculated from the fitted
model. The upper panel shows the residuals.

For this profile, the effective rate drops strongly as the CPMG pulse frequency
increases. That frequency dependence is CPMG relaxation dispersion: faster
refocusing increasingly suppresses the exchange contribution to transverse
relaxation. The corresponding 500 MHz profile is in `500mhz.pdf`; ChemEx fitted
both fields simultaneously, so the field dependence helps constrain the shared
population and exchange rate rather than treating each curve as an unrelated
fit.

The numerical experimental and calculated values behind the plots are also
available in `OutputTutorial/Plots/800mhz.exp`,
`OutputTutorial/Plots/800mhz.fit`, and `OutputTutorial/Data/800mhz.dat`.

## What ChemEx just did

ChemEx loaded the measured CPMG intensities, evaluated them under the selected
two-state exchange model, and globally optimized shared exchange parameters
together with profile-specific chemical-shift differences and relaxation
rates. It then wrote the fitted parameters, back-calculated data, plots,
fit-quality statistics, and run provenance.

## Next steps

- Learn how to configure [experiments](./user_guide/fitting/experiment_files.md)
  and [data files](./user_guide/fitting/data_files.mdx).
- Learn how to set [starting parameters and bounds](./user_guide/fitting/parameter_files.md).
- Learn how [Method v2](./user_guide/fitting/method_files.md) controls profile
  selection, parameter roles, and multi-step fits.
- Learn how to interpret and continue from [ChemEx outputs](./user_guide/fitting/outputs.mdx).
- See the complete [pure in-phase ¹⁵N CPMG experiment reference](./experiments/cpmg/cpmg_15n_ip.md).

The full `CPMG_15N_IP` example deliberately follows the five-profile global fit
with a second all-profile step that inherits parameter roles and fixes the
shared exchange parameters. That more robust staged strategy is useful in real
analyses, but is best learned after this first one-step fit.
