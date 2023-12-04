---
sidebar_label: ¹⁵N TROSY
sidebar_position: 4
description: '"cest_15n_tr"'
---

# ¹⁵N TROSY CEST

## Module name

`"cest_15n_tr"`

## Description

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

## Reference

-   D. Long, G. Bouvignies, and L.E. Kay. _Proc. Natl. Acad. Sci. USA_ **111**,
    8820-8825 (2014)

## Example

An application to the
[measurement of solvent exchange rates in invisible excited states](../../examples/trosy_cest.md)
is available
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N_TR/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cest_15n_tr'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cest_15n_tr"

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the ¹⁵N carrier during the CEST period, in ppm
carrier = 118.0

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

## Perform anti-trosy CEST experiment [optional, default: false]
# antitrosy = false

## State of the observed resonance [optional, default: "a"]
# observed_state = "a"

[conditions]

## ¹H Larmor frequency, in MHz
h_larmor_frq = 800.0

## Sample temperature, in Celsius [optional, depending on the kinetic model]
# temperature = 25.0

## Protein concentration, in M [optional, depending on the kinetic model]
# p_total = 500.0e-6

## Ligand concentration, in M [optional, depending on the kinetic model]
# l_total = 50.0e-6

## Labeling scheme of the sample, for deuterated samples "2H" should
## be used to obtain accurate initial estimates of relaxation rates
## based on model-free parameters, for uniformly ¹³C-labeled samples "13C"
## should be used to account for 1JCC properly [optional, default: []]
# label = ["2H", "13C"]

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
  G2N-HN = "G2N-HN.out"
  H3N-HN = "H3N-HN.out"
  K4N-HN = "K4N-HN.out"
  S5N-HN = "S5N-HN.out"
  L6N-HN = "L6N-HN.out"
```
