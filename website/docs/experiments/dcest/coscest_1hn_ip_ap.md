---
sidebar_label: In-phase/anti-phase amide ¹H cos-CEST
sidebar_position: 3
description: '"coscest_1hn_ip_ap"'
---

# In-phase/anti-phase amide ¹H cos-CEST

## Module name

`"coscest_1hn_ip_ap"`

## Description

Analyzes chemical exchange during the COS-CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

## References

-   T. Yuwen, G. Bouvignies, and L.E. Kay. _J. Mag. Reson._ **292**, 1-7 (2018)

## Example

An example for studying ¹⁵N-labeled sample is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/COSCEST_1HN_IP_AP/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'coscest_1hn_ip_ap'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "coscest_1hn_ip_ap"

## Recycle delay, in seconds
d1 = 0.1

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the carrier during the CEST period, in ppm
carrier = 8.3

## Cosine "spectral width", in Hz
sw = 800.0

## Number of excitation frequencies
cos_n = 3

## Number of points used to simulate the cosine-modulated shape
# cos_res = 10

## B1 radio-frequency field strength, in Hz
b1_frq = 25.0

## B1 inhomogeneity expressed as a fraction of 'b1_frq'. If set to "inf",
## a faster calculation takes place assuming full dephasing of the
## magnetization components orthogonal to the effective field. The "inf" value
## should not be used with an "eta_block" value larger than 0.
## [optional, default: 0.1]
# b1_inh_scale = 0.1

## Number of points used to simulate B1 inhomogeneity, the larger
## the longer the calculation. [optional, default: 11]
# b1_inh_res = 11

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

## List of plane indices to be excluded from the calculation.
## The first plane has index 0 and is usually the reference plane,
## this plane is always excluded by default [optional, default: [] ]
# filter_planes = []

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
