---
sidebar_label: Amide ¹H-¹⁵N longitudinal two-spin order relaxation
sidebar_position: 2
description: '"relaxation_hznz"'
---

# Amide ¹H-¹⁵N longitudinal two-spin order relaxation

## Module name

`"relaxation_hznz"`

## Description

Analyzes ¹H-¹⁵N longitudinal two-spin order relaxation experiments. Decay is
calculated using the (2n)×(2n), two-spin matrix, where n is the number of
states:

    \{ Iz(a), 2IzSz(a),
      Iz(b), 2IzSz(b), ... \}

## Reference

-   D.F. Hansen, D. Yang, H. Feng, Z. Zhou, S. Wiesner, Y. Bai, and L.E. Kay. _J.
    Am. Chem. Soc._ **129**, 11468-11479 (2007)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/RELAXATION_HZNZ/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'relaxation_hznz'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "relaxation_hznz"

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
  G2N-HN = "G2N-HN.out"
  H3N-HN = "H3N-HN.out"
  K4N-HN = "K4N-HN.out"
  S5N-HN = "S5N-HN.out"
  L6N-HN = "L6N-HN.out"
```
