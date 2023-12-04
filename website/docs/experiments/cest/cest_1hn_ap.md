---
sidebar_label: Pure anti-phase amide ¹H
sidebar_position: 5
description: '"cest_1hn_ap"'
---

# Pure anti-phase amide ¹H CEST

## Module name

`"cest_1hn_ap"`

## Description

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

## Reference

-   A. Sekhar, R. Rosenzweig, G. Bouvignies, and L.E. Kay. _Proc. Natl. Acad. Sci.
    USA_ **113**, E2794-E2801 (2016)

## Example

An application for measuring PRE values in a protein excited state is available
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_1HN_AP/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cest_1hn_ap'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cest_1hn_ap"

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the carrier during the CEST period, in ppm
carrier = 8.3

## B1 radio-frequency field strength, in Hz
b1_frq = 25.0

## B1 inhomogeneity expressed as a fraction of 'b1_frq'. If set to "inf",
## a faster calculation takes place assuming full dephasing of the
## magnetization components orthogonal to the effective field.
## [optional, default: 0.1]
# b1_inh_scale = 0.1

## Number of points used to simulate B1 inhomogeneity, the larger
## the longer the calculation. [optional, default: 11]
# b1_inh_res = 11

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
## "scatter": uncertainties are calculated from the baseline
## [optional, default: "file"]
# error = "file"

## List of offsets relative to the main resonance position
## (nu) and bandwidths (delta_nu) defining regions where
## points are excluded from the calculation (nu +/- 0.5 * delta_nu),
## both are in Hz [optional, default: [[0.0, 0.0]] ]
# filter_offsets = [
#   [0.0, 0.0],
# ]

  ## List of the profile names and their associated filenames.
  ## The name of the spin systems should follow the Sparky-NMR
  ## conventions.
  [data.profiles]
  G2HN-N = "G2N-HN.out"
  H3HN-N = "H3N-HN.out"
  K4HN-N = "K4N-HN.out"
  S5HN-N = "S5N-HN.out"
  L6HN-N = "L6N-HN.out"
```
