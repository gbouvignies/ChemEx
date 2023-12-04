---
sidebar_label: Pure in-phase ¹⁵N with [0013] phase cycle
sidebar_position: 2
description: '"cpmg_15n_ip_0013"'
---

# Pure in-phase ¹⁵N CPMG with [0013] phase cycle

## Module name

`"cpmg_15n_ip_0013"`

## Description

Analyzes ¹⁵N chemical exchange in the presence of high power ¹H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)×(3n), single-spin matrix, where n is the number
of states:

    \{ Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... \}

This version is modified such that CPMG pulses are applied with [0013] phase
cycle as shown in Daiwen's paper. The CW decoupling on ¹H is assumed to be
strong enough (> 15 kHz) such that perfect ¹H decoupling can be achieved.

## Reference

-   B. Jiang, B. Yu, X. Zhang, M. Liu, and D. Yang. _J. Magn. Reson._ **257**, 1-7
    (2015)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_15N_IP_0013/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_15n_ip_0013'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_15n_ip_0013"

## CPMG relaxation delay, in seconds
time_t2 = 0.04

## Position of the ¹⁵N carrier during the CPMG period, in ppm
carrier = 118.0

## ¹⁵N 90 degree pulse width of CPMG pulses, in seconds
pw90 = 35.0e-6

## Maximum number of cycles
ncyc_max = 40

## Equilibration delay at the end of the CPMG period, in seconds
## [optional, default: 0.0]
# time_equil = 0.0

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
