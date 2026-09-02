---
sidebar_position: 5
---

# Parameter Files

## Overview

Parameter files contain initial estimates for the parameters used in the fitting process. These files are specified in ChemEx using the `-p` or `--parameters` option:

```shell
chemex fit [...] -p <parameter_file> [...]
```

Parameter files are organized into sections:

-   The `[GLOBAL]` section applies parameters universally to all residues.
-   Residue-specific parameters are defined in sections named after each parameter, such as `[CS_A]`. Multiple parameter files can be provided if needed.

:::warning
To ensure accurate results and avoid local minima, set appropriate initial values for each parameter, as the χ<sup>2</sup> minimization process involves multidimensional searching.
:::

:::info
If no initial value is provided in the parameter files, a default value will be assigned.
:::

## Example Parameter File

Below is an example of a parameter file:

```toml title="parameters.toml"
[GLOBAL]
PB     = 0.015
KEX_AB = 70.0
TAUC_A = 10.0

[CS_A]
13N = 108.207
26N = 115.711
28N = 113.882
29N = 115.318
33N = 115.636
37N = 116.159
41N = 114.635
42N = 113.525
43N = 108.876
50N = 107.855
52N = 111.358
55N = 128.301
59N = 116.388
66N = 119.429
67N = 114.454
68N = 120.595

[DW_AB]
13N = 4.0
26N = 5.5
28N = 6.5
29N = 6.0
33N = 4.5
37N = 5.0
41N = 6.0
42N = 6.0
43N = 12.5
50N = 8.0
52N = 8.5
55N = -6.5
59N = 6.5
66N = 4.0
67N = 8.0
68N = 4.5
```

:::tip
Setting model-free parameters (e.g., `TAUC_A`) can provide a good initial estimate for relaxation parameters (e.g., `R1_A`, `R2_A`). For biomolecules in H<sub>2</sub>O at T = 300 K, the overall tumbling time is roughly 1 ns per 2.6 kDa of molecular weight. Tumbling time is proportional to η/T, where η is solution viscosity and T is temperature in Kelvin.
:::

## Setting Parameter Bounds

To set upper and lower bounds for any fitting parameter, replace the initial value with a list of three elements:

```toml title="parameters.toml"
PARAMETER_WITH_NO_BOUNDS = <initial_value>
PARAMETER_WITH_BOUNDS = [<initial_value>, <lower_bound>, <upper_bound>]
```

Setting bounds helps prevent parameters from reaching unrealistic values during χ<sup>2</sup> minimization. However, avoid overly strict bounds as they can hinder convergence. Some minimization algorithms, such as `differential_evolution`, require finite bounds for all fitted parameters.
