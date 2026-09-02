---
sidebar_label: ¹⁵N exchange induced shifts with ¹⁵N–¹H HSQC/HMQC
sidebar_position: 1
description: '"shift_15n_sqmq"'
---

# ¹⁵N exchange induced shifts with ¹⁵N–¹H HSQC/HMQC

## Module name

`"shift_15n_sqmq"`

## Description

Analyzes exchange induced ¹⁵N chemical shift changes measured in (¹⁵N–¹HN) HMQC
and HSQC data sets.

:::note

Since this experiment is used for determining the sign of Δϖ, it is usually
combined with other CPMG experiments.

:::

## References

- N.R. Skrynnikov, F.W. Dahlquist, L.E. Kay. _J. Am. Chem. Soc._ **124**,
  12352-12360 (2002)
- P. Vallurupalli, G. Bouvignies, and L.E. Kay. _J. Phys. Chem. B_ **115**,
  14891-14900 (2011)

## Example

An example use of the module associated with ¹⁵N and ¹H CPMG datasets is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/Shifts/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'shift_15n_sqmq'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "shift_15n_sqmq"

## State of the observed resonance [optional, default: "a"]
# observed_state = "a"

[conditions]

## 1H Larmor frequency, in MHz
h_larmor_frq = 800.0

## Sample temperature, in Celsius [optional, depending on the kinetic model]
# temperature = 25.0

## Protein concentration, in M [optional, depending on the kinetic model]
# p_total = 500.0e-6

## Ligand concentration, in M [optional, depending on the kinetic model]
# l_total = 50.0e-6

[data]

## Directory containing the profiles [optional, default: "./"]
# path = "./"

## Filename of the file containing the list of the shifts in ppb.
## The file should be formatted as follow:
##
## #     name   shift  error
##     G2N-HN    10.9    0.5
##     H3N-HN    32.1    0.5
##     K4N-HN   -54.3    1.5
##     S5N-HN     0.7    0.5
##     L6N-HN   -15.2    0.5
##
## The name of the spin systems should follow the Sparky-NMR
## conventions.
shifts = "sqmq.txt"
```
