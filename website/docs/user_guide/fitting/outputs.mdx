---
sidebar_position: 8
---

# Outputs

## Location

The results of the fits are written in the directory `Output/` by default.
However, you can change this location using the command-line option `--output`
(or `-o`) followed by the desired path name.

```bash
chemex fit -o <path> [other options]
```

### Multi-step fit

If the fitting method has multiple fitting steps, each step will create its own
output subdirectory with the name of the fitting step.

```shell
Output/
├── STEP1/
└── STEP2/
```

### Multi-group fit

During any fitting step, if the dataset can be divided in multiple groups
depending on distinct sets of fitting parameters, ChemEx will fit each group of
data separately. The results of these individual fits are then stored in
separate subfolders placed in the the `Groups/` folder. Subfolders are named
with a number and an ID, which depends on the parameters that have been
optimized for the specific group. The results are also put together in the
`All/` folder for convenience.

For example, if all global parameters (p<sub>B</sub>, k<sub>ex</sub>, etc.) are
fixed or residue-specific fitting model is used (e.g. the model
[`2st_rs`](kinetic_models.md)), then all parameters are allowed to vary are
residue-specific. Residue-specific fits will then be run and two separate
subdirectories `All/` and `Groups/` will be created, which contain fitting
results for all residues and each individual residue, respectively.

```shell
Output/
└── STEP2/
    ├── All/
    └── Groups/
        ├── 10_11N/
        ├── 11_12N/
        ├── 12_13N/
```

## Content

The fitting output typically contains the following files and directories:

### `Parameters/`

Contains fitting results as three separate files `fitted.toml`, `fixed.toml` and
`constrained.toml`, which contain output parameters that are fitted, fixed and
constrained during the fitting process, respectively.

#### Example files

import Tabs from '@theme/Tabs'; import TabItem from '@theme/TabItem';

<Tabs>
<TabItem value="fitted" label="fitted.toml" default>

```toml
[GLOBAL]
KEX_AB =  3.81511e+02 # ±8.90870e+00
PB     =  7.02971e-02 # ±1.14784e-03

[DW_AB]
15N =  2.00075e+00 # ±2.30817e-02
31N =  1.98968e+00 # ±1.90842e-02
33N =  1.82003e+00 # ±2.44821e-02
34N =  3.63170e+00 # ±3.62801e-02
37N =  1.69183e+00 # ±2.41070e-02

["R2_A, B0->500.0MHZ"]
15N =  3.98674e+00 # ±2.55793e-01
31N =  5.85923e+00 # ±1.94734e-01
33N =  4.02099e+00 # ±2.78003e-01
34N =  4.16615e+00 # ±2.05190e-01
37N =  3.67705e+00 # ±3.04621e-01

["R2_A, B0->800.0MHZ"]
15N =  6.33712e+00 # ±3.86894e-01
31N =  7.99927e+00 # ±2.81972e-01
33N =  6.23967e+00 # ±5.82366e-01
34N =  5.99535e+00 # ±3.17832e-01
37N =  5.37295e+00 # ±5.69929e-01
```

:::note

If uncertainties on fitted parameters have been calculated – this depends on the
fitting algorithm used –, they are reported as comments at the end of the line
and are preceded by the sign "±". Uncertainties on fitted parameters are,
typically, estimated through the covariance matrix obtained from the
Levenberg-Marquardt optimization.

:::

</TabItem>
<TabItem value="fixed" label="fixed.toml">

```toml
[CS_A]
15N =  1.19849e+02 # (fixed)
31N =  1.26388e+02 # (fixed)
33N =  1.18762e+02 # (fixed)
34N =  1.14897e+02 # (fixed)
37N =  1.21618e+02 # (fixed)

["R1_A, B0->500.0MHZ"]
15N =  1.50000e+00 # (fixed)
31N =  1.50000e+00 # (fixed)
33N =  1.50000e+00 # (fixed)
34N =  1.50000e+00 # (fixed)
37N =  1.50000e+00 # (fixed)

["R1_A, B0->800.0MHZ"]
15N =  1.50000e+00 # (fixed)
31N =  1.50000e+00 # (fixed)
33N =  1.50000e+00 # (fixed)
34N =  1.50000e+00 # (fixed)
37N =  1.50000e+00 # (fixed)
```

</TabItem>
<TabItem value="constrained" label="constrained.toml">

```toml
[GLOBAL]
KAB =  2.68192e+01 # ±3.06068e-01 ([KEX_AB] * [PB])
KBA =  3.54692e+02 # ±8.28245e+00 ([KEX_AB] * [PA])
PA  =  9.29703e-01 # ±1.14784e-03 (1.0 - [PB])

[CS_B]
15N =  1.21850e+02 # ±2.30817e-02 ([CS_A, NUC->15N] + [DW_AB, NUC->15N])
31N =  1.28378e+02 # ±1.90842e-02 ([CS_A, NUC->31N] + [DW_AB, NUC->31N])
33N =  1.20582e+02 # ±2.44821e-02 ([CS_A, NUC->33N] + [DW_AB, NUC->33N])
34N =  1.18529e+02 # ±3.62801e-02 ([CS_A, NUC->34N] + [DW_AB, NUC->34N])
37N =  1.23310e+02 # ±2.41070e-02 ([CS_A, NUC->37N] + [DW_AB, NUC->37N])

["R1_B, B0->500.0MHZ"]
15N =  1.50000e+00 # ([R1_A, NUC->15N, B0->500.0MHZ])
31N =  1.50000e+00 # ([R1_A, NUC->31N, B0->500.0MHZ])
33N =  1.50000e+00 # ([R1_A, NUC->33N, B0->500.0MHZ])
34N =  1.50000e+00 # ([R1_A, NUC->34N, B0->500.0MHZ])
37N =  1.50000e+00 # ([R1_A, NUC->37N, B0->500.0MHZ])

["R1_B, B0->800.0MHZ"]
15N =  1.50000e+00 # ([R1_A, NUC->15N, B0->800.0MHZ])
31N =  1.50000e+00 # ([R1_A, NUC->31N, B0->800.0MHZ])
33N =  1.50000e+00 # ([R1_A, NUC->33N, B0->800.0MHZ])
34N =  1.50000e+00 # ([R1_A, NUC->34N, B0->800.0MHZ])
37N =  1.50000e+00 # ([R1_A, NUC->37N, B0->800.0MHZ])

["R2_B, B0->500.0MHZ"]
15N =  3.98674e+00 # ±2.55793e-01 ([R2_A, NUC->15N, B0->500.0MHZ])
31N =  5.85923e+00 # ±1.94734e-01 ([R2_A, NUC->31N, B0->500.0MHZ])
33N =  4.02099e+00 # ±2.78003e-01 ([R2_A, NUC->33N, B0->500.0MHZ])
34N =  4.16615e+00 # ±2.05190e-01 ([R2_A, NUC->34N, B0->500.0MHZ])
37N =  3.67705e+00 # ±3.04621e-01 ([R2_A, NUC->37N, B0->500.0MHZ])

["R2_B, B0->800.0MHZ"]
15N =  6.33712e+00 # ±3.86894e-01 ([R2_A, NUC->15N, B0->800.0MHZ])
31N =  7.99927e+00 # ±2.81972e-01 ([R2_A, NUC->31N, B0->800.0MHZ])
33N =  6.23967e+00 # ±5.82366e-01 ([R2_A, NUC->33N, B0->800.0MHZ])
34N =  5.99535e+00 # ±3.17832e-01 ([R2_A, NUC->34N, B0->800.0MHZ])
37N =  5.37295e+00 # ±5.69929e-01 ([R2_A, NUC->37N, B0->800.0MHZ])
```

:::note

Propagated uncertainties – when available – and applied constrained are reported
at the end of the line as comments.

:::

</TabItem>
</Tabs>

### `Plot/`

Contains fitting results as plots (in `.pdf` format) and also the raw datasets
(both the original input and fitted data points) for creating the plots. Example
fitting results for CPMG and CEST experiments are shown below:

import CestProfile from '@site/static/img/cest_26hz_fit.png'; import CpmgProfile
from '@site/static/img/cpmg_800mhz_fit.png';

<figure>
  <img src={CestProfile} alt="CEST profile" width="50%"/>
  <img src={CpmgProfile} alt="CPMG profile" width="50%"/>
  <figcaption align="center"><b>Examples of CEST and CPMG fitting results</b></figcaption>
</figure>

:::note

In plots of (D-/cos-)CEST fitting results, the positions for ground and excited
states are indicated by solid and dashed vertical lines respectively. Besides,
data points that are filtered out from the fit are shown in lighter color. In
plots of D-CEST/COS-CEST fitting results, the "folded" positions for ground and
excited states are indicated by "\*" at the vertical lines.

:::

### `Data/`

Contains all the data values used for the fitting along with the back-calculated
values. These files can be used to calculate the $χ^2$ value.

```toml title=Data/500mhz.dat
[15N]
#         NCYC   INTENSITY (EXP)       ERROR (EXP)  INTENSITY (CALC)
             0    3.47059800e+04    1.77491406e+02    3.47055362e+04
            30    3.05930380e+04    1.77491406e+02    3.06111963e+04
             1    1.81234230e+04    1.77491406e+02    1.80856300e+04
            28    3.07144730e+04    1.77491406e+02    3.05838767e+04
             2    2.02155120e+04    1.77491406e+02    2.02110184e+04
            26    3.05222020e+04    1.77491406e+02    3.05501255e+04
             3    2.23056070e+04    1.77491406e+02    2.23601656e+04
            24    3.05381830e+04    1.77491406e+02    3.05077686e+04
             4    2.43783050e+04    1.77491406e+02    2.44111932e+04
            22    3.06981570e+04    1.77491406e+02    3.04536377e+04
             5    2.58673980e+04    1.77491406e+02    2.59884541e+04
            20    3.06069180e+04    1.77491406e+02    3.03829761e+04
             6    2.70660870e+04    1.77491406e+02    2.71148924e+04
            18    3.00982700e+04    1.77491406e+02    3.02883916e+04
             7    2.81512990e+04    1.77491406e+02    2.79178594e+04
            16    3.02515700e+04    1.77491406e+02    3.01579211e+04
             8    2.84045570e+04    1.77491406e+02    2.84985729e+04
            14    2.97528530e+04    1.77491406e+02    2.99712519e+04
             9    2.84536650e+04    1.77491406e+02    2.89264158e+04
            13    2.98219180e+04    1.77491406e+02    2.98461022e+04
            10    2.96936410e+04    1.77491406e+02    2.92495591e+04
            12    2.95269900e+04    1.77491406e+02    2.96918706e+04
            11    2.92540210e+04    1.77491406e+02    2.94972506e+04
             2    2.04476470e+04    1.77491406e+02    2.02110184e+04
            28    3.05466550e+04    1.77491406e+02    3.05838767e+04
             8    2.86183900e+04    1.77491406e+02    2.84985729e+04

[31N]
#         NCYC   INTENSITY (EXP)       ERROR (EXP)  INTENSITY (CALC)
             0    4.71577550e+04    1.77491406e+02    4.71537007e+04
            30    4.00023250e+04    1.77491406e+02    3.94823953e+04
             1    2.30167260e+04    1.77491406e+02    2.34136180e+04
            28    4.00615990e+04    1.77491406e+02    3.94547546e+04
             2    2.61934530e+04    1.77491406e+02    2.61653780e+04
            26    3.96898190e+04    1.77491406e+02    3.94190638e+04
             3    2.86882570e+04    1.77491406e+02    2.88729551e+04
            24    3.98238690e+04    1.77491406e+02    3.93722573e+04
             4    3.21031440e+04    1.77491406e+02    3.13942278e+04
            22    3.96733310e+04    1.77491406e+02    3.93097577e+04
             5    3.35871430e+04    1.77491406e+02    3.34628362e+04
            20    3.86957600e+04    1.77491406e+02    3.92245398e+04
             6    3.53600830e+04    1.77491406e+02    3.49641428e+04
            18    3.88115960e+04    1.77491406e+02    3.91054958e+04
             7    3.58938150e+04    1.77491406e+02    3.59453347e+04
            16    3.85843960e+04    1.77491406e+02    3.89345050e+04
             8    3.61450060e+04    1.77491406e+02    3.66900206e+04
            14    3.85711880e+04    1.77491406e+02    3.86811094e+04
             9    3.72815770e+04    1.77491406e+02    3.72429094e+04
            13    3.82358560e+04    1.77491406e+02    3.85092366e+04
            10    3.76426640e+04    1.77491406e+02    3.76803042e+04
            12    3.78470460e+04    1.77491406e+02    3.82929922e+04
            11    3.76000060e+04    1.77491406e+02    3.80220076e+04
             2    2.62301400e+04    1.77491406e+02    2.61653780e+04
            28    3.97731200e+04    1.77491406e+02    3.94547546e+04
             8    3.63510140e+04    1.77491406e+02    3.66900206e+04

[...]
```

### `statistics.toml`

Contains all goodness-of-fit statistics, such as $χ^2$.

```toml title='"statistics.toml"'
"number of data points"                = 230
"number of variables"                  = 17
"chi-square"                           =  4.34824e+02
"reduced-chi-square"                   =  2.04143e+00
"chi-squared test"                     =  0.00000e+00
"Kolmogorov-Smirnov test"              =  6.72236e-02
"Akaike Information Criterion (AIC)"   =  1.80479e+02
"Bayesian Information Criterion (BIC)" =  2.38926e+02

```
