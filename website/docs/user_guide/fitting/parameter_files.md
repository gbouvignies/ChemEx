---
sidebar_position: 5
---

# Parameter files

## Description

The parameter files contain initial estimates of the parameters to be used
during the fitting process. These files are provided to ChemEx using the `-p` or
`--parameters` option.

```shell
chemex fit [...] -p <parameter_file> [...]
```

Parameter files are divided in multiple sections:

- The parameter values under the section `[GLOBAL]` apply to all residues.
- Residue-specific parameters are specified under sections named after the
  parameter name, such as `[CS_A]`. Multiple parameter files can be provided if
  necessary.

:::caution

Due to the multidimensional searching feature of $χ^2$ minimization process, it
is essential to set a suitable initial value for each parameter to avoid getting
trapped in a local minimum.

:::

:::info

When no starting value is provided in the parameter files, a default value is
assigned as the initial value.

:::

## Example

Here is an example parameter file:

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

Setting model-free parameters (e.g. `TAUC_A`) is a simple way to obtain initial
estimates of relaxation parameters (e.g., `R1_A`, `R2_A`, etc.). For every 2.6
kDa molecular weight, the overall tumbling time is approximately 1 ns at T = 300
K for biomolecules in H<sub>2</sub>O. Assuming similar molecular structure at
different conditions, the overall tumbling time is proportional to η/T, where η
is solution viscosity and T is temperature in Kelvin.

:::

## Setting bounds

You can set upper and lower bounds to any fitting parameters by replacing the
initial value by a list of three elements:

```toml title="parameters.toml"
PARAMETER_WITH_NO_BOUNDS = <initial_value>
PARAMETER_WITH_BOUNDS = [<initial_value>, <lower_bound>, <upper_bound>]
```

Such boundaries can help to prevent parameters wandering off to unrealistic
values to minimize $χ^2$. However, one should be careful not to set too
stringent boundaries either, as this can result in convergence problems. Certain
minimization algorithms (e.g. `differential_evolution`) require finite bounds on
all fitted parameters.
