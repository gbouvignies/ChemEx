---
sidebar_position: 2.5
sidebar_label: Build a full CPMG analysis
title: From Your First Fit to a Full CPMG Analysis
description: Extend the first two-field 15N CPMG fit to ChemEx's staged, all-profile Method-v2 workflow while preserving shared exchange parameters.
---

# From your first fit to a full CPMG analysis

This tutorial continues from [Run Your First ChemEx Analysis](./first_analysis.md).
That first fit used both magnetic fields for five residues and allowed the shared
exchange parameters and local profile parameters to vary together. Here you will
run the complete staged analysis already shipped in the same `CPMG_15N_IP`
example.

The staged method separates two practical tasks:

1. estimate shared exchange parameters from a selected, exchange-sensitive
   subset;
2. keep those shared values fixed while fitting the local response of the
   complete dataset.

This is one useful analysis strategy, not a rule that every ChemEx analysis must
follow.

## The full Method file

Open `Methods/method.toml` in
`examples/Experiments/CPMG_15N_IP`. The complete file is:

```toml title="Methods/method.toml"
FORMAT_VERSION = 2

[STEP1]
INCLUDE = ["15", "31", "33", "34", "37"]

[STEP2]
INCLUDE = "ALL"
ROLES_FROM = "STEP1"
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
]
```

`FORMAT_VERSION = 2` selects the current Method-v2 syntax. ChemEx then executes
the two tables in declaration order: `STEP1`, followed by `STEP2`. Each step has
its own profile selection and resolved parameter roles, while fitted numerical
values continue from one successfully committed step to the next.

```text
STEP1
5 selected residues × 2 fields = 10 profiles
        ↓
fit shared PB and KEX_AB
plus residue- and field-specific parameters
        ↓ committed fitted values
STEP2
54 residues × 2 fields = 108 profiles
        ↓
hold PB and KEX_AB at the STEP1 values
fit the remaining local parameters
```

## STEP1: establish the shared exchange parameters

```toml
[STEP1]
INCLUDE = ["15", "31", "33", "34", "37"]
```

`INCLUDE` selects spin systems for this step. Because residue 15, for example,
is present in both experiment files, ChemEx includes its 500 and 800 MHz
profiles. The five residue numbers therefore select ten profiles and 230 data
points in total.

STEP1 does not contain a `ROLES` array, so the model and experiment use their
baseline parameter roles. In this example ChemEx varies 17 parameters:

- the two global parameters `PB` and `KEX_AB`;
- one residue-specific `DW_AB` for each of the five residues;
- one residue- and field-specific `R2_A` for each of the ten profiles.

Other required quantities keep their baseline roles. For example, the supplied
state-A chemical shifts (`CS_A`) and calculated `R1_A` values are fixed, while
model-owned quantities such as `PA`, `KAB`, and `KBA` remain derived.

### Global and local parameters in this example

| Parameter | Scope in this analysis                              | Role in STEP1 | Practical meaning                                       |
| --------- | --------------------------------------------------- | ------------- | ------------------------------------------------------- |
| `PB`      | Shared across both fields and all selected residues | Fitted        | Population fraction of state B                          |
| `KEX_AB`  | Shared across both fields and all selected residues | Fitted        | Total A↔B exchange rate, in s⁻¹                         |
| `DW_AB`   | One value per residue                               | Fitted        | Residue-specific B-minus-A ¹⁵N shift difference, in ppm |
| `R2_A`    | One value per residue and magnetic field            | Fitted        | Profile baseline transverse relaxation rate, in s⁻¹     |

The shared parameters can be informed simultaneously by several profiles
because the analysis assumes that they report on the same exchange process.
`DW_AB` describes how an individual nucleus responds to that process, while
`R2_A` also reflects its particular residue and dataset.

The five chosen residues all show clear exchange-sensitive dispersion in the
shipped data, so they provide a practical subset from which to establish the
shared parameters. The repository does not define a formal selection criterion
or demonstrate that these are objectively the best five residues. Treat the
list as a deliberate choice for this example, not a universal recipe.

A current run gives approximately:

| STEP1 result |     Value |
| ------------ | --------: |
| `PB`         |   0.07030 |
| `KEX_AB`     | 381.7 s⁻¹ |
| χ²           |     434.6 |
| Reduced χ²   |     2.040 |

### How STEP1 relates to the first tutorial

STEP1 and `Methods/tutorial.toml` resolve the same ten profiles and the same 17
fitted parameters. Starting from the same parameter file, they produce the same
fit. The visible differences are organizational: the beginner Method calls its
only step `FIT` and writes results directly under `OutputTutorial`, whereas the
two-step Method writes this result under `OutputFull/STEP1` before continuing.

## STEP2: extend the model to every profile

```toml
[STEP2]
INCLUDE = "ALL"
ROLES_FROM = "STEP1"
ROLES = [
  { FIX = ["PB", "KEX_AB"] },
]
```

Each line has a separate job:

- `INCLUDE = "ALL"` selects the full profile set for this step. Selection is
  step-local; it is not inherited from STEP1. The current example has 54 spin
  systems at each of two fields, so STEP2 uses 108 profiles—98 more than STEP1.
- `ROLES_FROM = "STEP1"` reuses only STEP1's effective parameter-role and
  constraint setup. It does not copy profile selection, search settings,
  statistics requests, or numerical values.
- The `FIX` action is applied after that inherited setup. It changes the complete
  role of `PB` and `KEX_AB` from fitted to fixed for STEP2.

Numerical continuity is related but distinct: every later step starts from the
latest successfully committed parameter values whether or not `ROLES_FROM` is
present. In this run, STEP2 therefore starts with STEP1's fitted `PB = 0.0702971`
and `KEX_AB = 381.709 s⁻¹`; `FIX` prevents either value from moving during the
second fit.

There is an important detail in this compact Method file. STEP1 has no explicit
`ROLES` actions, so its effective setup is simply the baseline setup. In this
specific case, `ROLES_FROM` introduces no additional custom action that STEP2
would otherwise lack. Its exact meaning remains role inheritance—not result
inheritance—and it makes the intended relationship between the steps explicit.

### What still varies in STEP2

The baseline fitted roles for `DW_AB` and `R2_A` are preserved. With all profiles
active, ChemEx fits:

- 54 residue-specific `DW_AB` values;
- 108 residue- and field-specific `R2_A` values.

That gives 162 varying parameters. `PB` and `KEX_AB` appear in
`STEP2/Parameters/fixed.toml` at their STEP1 values. The `CS_A` and `R1_A`
parameters remain fixed as before, and model-owned derived relationships remain
derived.

For the five residues already present in STEP1, the local parameters also start
from their committed fitted values and remain free to adjust. Parameters for the
49 newly included residues were not changed by STEP1, so they start from the
original invocation values. In the shipped parameter setup, that means
`DW_AB = 2.0 ppm`, with initial `R2_A` values derived from `TAUC_A` (about
4.56 s⁻¹ at 500 MHz and 5.27 s⁻¹ at 800 MHz), before STEP2 optimizes them.

## Run the full workflow

From `ChemEx/examples/Experiments/CPMG_15N_IP`, run:

```shell
chemex fit -e Experiments/500mhz.toml Experiments/800mhz.toml -p Parameters/parameters.toml -m Methods/method.toml -d 2st -o OutputFull
```

The terminal first reports `Selecting profiles -> 10 profiles` for STEP1 and
then `Selecting profiles -> 108 profiles` for STEP2. A completed two-step run
has this top-level layout:

```text
OutputFull/
├── run_info/
│   └── outcome.toml
├── STEP1/
│   ├── Parameters/
│   ├── Data/
│   ├── Plots/
│   ├── Statistics/
│   └── statistics.toml
└── STEP2/
    ├── Parameters/
    ├── Data/
    ├── Plots/
    ├── Statistics/
    └── statistics.toml
```

Confirm that `OutputFull/run_info/outcome.toml` contains
`status = "complete"`. Then open `OutputFull/STEP2/Parameters/fixed.toml` to
verify the fixed global values and `OutputFull/STEP2/Parameters/fitted.toml` to
inspect the fitted local values.

## What STEP2 adds scientifically

STEP2 asks whether each newly included spin system can be described using the
exchange population and rate established in STEP1, while allowing that spin
system's `DW_AB` and its two field-specific `R2_A` values to adjust.

For example, residues 5 and 32 are absent from STEP1 but have clear dispersion
in the STEP2 plots. A current run fits the following local values while keeping
the same fixed `PB` and `KEX_AB` for both:

| Newly included residue | `DW_AB` (ppm) | `R2_A`, 500 MHz (s⁻¹) | `R2_A`, 800 MHz (s⁻¹) |
| ---------------------- | ------------: | --------------------: | --------------------: |
| 5N                     |         0.465 |                  4.32 |                  5.03 |
| 32N                    |          7.15 |                  5.23 |                  7.97 |

Inspect pages `5N` and `32N` in
`OutputFull/STEP2/Plots/500mhz.pdf` and `800mhz.pdf`. Their different local
parameters describe different profile responses without changing the shared
exchange process. This is the practical purpose of the second stage.

## Compare the two workflows

|                          | First tutorial                            | Full shipped workflow                                             |
| ------------------------ | ----------------------------------------- | ----------------------------------------------------------------- |
| Method file              | `Methods/tutorial.toml`                   | `Methods/method.toml`                                             |
| Profiles                 | Five residues at two fields (10 profiles) | STEP1 uses those 10; STEP2 uses all 108                           |
| Method steps             | One                                       | Two                                                               |
| Shared `PB` and `KEX_AB` | Fitted                                    | Fitted in STEP1, fixed at those values in STEP2                   |
| Local `DW_AB` and `R2_A` | Fitted for the selected residues          | Fitted for every included residue/profile                         |
| Purpose                  | Learn and inspect one global fit          | Apply an established exchange description to the complete dataset |

## Read the diagnostics proportionately

A current run completes both steps. STEP1 reports covariance available, with a
full-rank 17-parameter Jacobian. STEP2 gives approximately χ² = 2694.7 and
reduced χ² = 1.160 for 2484 data points and 162 varying parameters. Its
covariance is also available and full rank, but ChemEx reports boundary
warnings: the fitted `DW_AB` for 27N lies near a fitted-coordinate boundary, so
its symmetric covariance uncertainty may be misleading.

These statements answer different questions. `status = "complete"` means every
requested fit step and output operation finished. The covariance and boundary
evidence describes how the local uncertainty approximation behaved at the
accepted result. Completion does not erase a diagnostic that deserves
inspection, and a non-fatal boundary warning does not by itself mean that the
staged workflow failed.

See [Outputs](./user_guide/fitting/outputs.mdx#statistics) for the generated
`Statistics/Covariance/evidence.json` and the meaning of parameter annotations.

## Adapt the strategy to your data

Before applying this pattern elsewhere, reconsider the scientific decisions it
encodes:

- Which profiles, if any, form an informative initial subset for the exchange
  parameters you want to estimate?
- Does your dataset benefit from staging, or is one simultaneous global fit more
  appropriate?
- Which parameters genuinely describe one shared process, and which must remain
  residue- or dataset-specific?
- Is it scientifically justified to fix selected shared parameters in a later
  step, given their uncertainty and the behavior of the added profiles?

Do not automatically choose the strongest-looking dispersions or copy this
`FIX` list. Profile selection and parameter fixing are analysis decisions that
should follow the experiment, model, and question being tested.

## Next steps

- Return to [Run Your First ChemEx Analysis](./first_analysis.md).
- Read the complete [Method-v2 reference](./user_guide/fitting/method_files.md).
- Learn how to set [starting parameters and bounds](./user_guide/fitting/parameter_files.md).
- Learn how to inspect [ChemEx outputs and evidence](./user_guide/fitting/outputs.mdx).
- See the [pure in-phase ¹⁵N CPMG experiment reference](./experiments/cpmg/cpmg_15n_ip.md).
