---
sidebar_label: ¹⁵N longitudinal relaxation
sidebar_position: 1
description: '"relaxation_nz"'
---

# ¹⁵N longitudinal relaxation

## Module name

`"relaxation_nz"`

## Description

Analyzes ¹⁵N T1 experiments. This keeps the spin system purely in-phase
throughout, and is calculated using the (1n)×(1n), single-spin matrix, where n
is the number of states:

    \{ Iz(a), Iz(b), ... \}

## Reference

-   L.E. Kay, L.K. Nicholson, F. Delaglio, A. Bax, and D.A. Torchia. _J. Mag.
    Reson._ **97**, 359-375 (1992)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/RELAXATION_NZ/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'relaxation_nz'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "relaxation_nz"

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

## Labeling scheme of the sample, for deuterated samples "2H" should
## be used to obtain accurate initial estimates of relaxation rates
## based on model-free parameters [optional, default: []]
# label = ["2H"]

[data]

## Directory containing the profiles [optional, default: "./"]
# path = "./"

## Option defining how intensity uncertainties are estimated.
## "file": uncertainties are taken from the profile files
## "duplicates": uncertainties are calculated from the duplicate points
## [optional, default: "file"]
# error = "file"

  ## List of the profile names and their associated filenames.
  ## The name of the spin systems should follow the Sparky-NMR
  ## conventions.
  [data.profiles]
  G2N = "G2N-HN.out"
  H3N = "H3N-HN.out"
  K4N = "K4N-HN.out"
  S5N = "S5N-HN.out"
  L6N = "L6N-HN.out"
```
