---
sidebar_label: Pure in-phase ¹³C cos-CEST
sidebar_position: 4
description: Analyze pure in-phase ¹³C cos-CEST chemical exchange data with ¹H composite decoupling.
---

# Pure in-phase ¹³C cos-CEST

## Module name

`"coscest_13c"`

## Description

Analyzes chemical exchange during a cosine-modulated multi-site ¹³C CEST
(cos-CEST) block in the presence of ¹H composite decoupling. Magnetization
evolution is calculated using the (3n)×(3n), single-spin matrix, where n is the
number of states:

    \{ Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... \}

## Reference

-   T. Yuwen, G. Bouvignies, and L.E. Kay. _J. Magn. Reson._ **292**, 1-7 (2018)

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'coscest_13c'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "coscest_13c"

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Equilibration delay at the end of the CEST period, in seconds
## [optional, default: 0.0]
# time_equil = 0.0

## Position of the ¹³C carrier during the CEST period, in ppm
carrier = 55.0

## Cosine "spectral width", in Hz
sw = 800.0

## Number of excitation frequencies
cos_n = 3

## Number of points used to simulate the cosine-modulated shape
# cos_res = 10

## B1 radio-frequency field strength, in Hz (or provide pw90 instead, not both)
b1_frq = 25.0

## B1 inhomogeneity distribution (required; "dephasing" is not valid for cos-CEST)
[experiment.b1_distribution]
type = "gaussian"
scale = 0.1
# res = 15  # optional override (automatic by default)

## Initial condition: equilibrium preparation (all states populated according to
## their equilibrium populations). Set start_state for non-equilibrium
## preparation when chemical shift evolution occurs before the measured
## relaxation/exchange element (see Yuwen et al., J. Biomol. NMR 2016,
## 65:143-156).
## [optional, default: equilibrium]
# start_state = "a"

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

## Labeling scheme of the sample, for uniformly ¹³C-labeled samples "13C"
## should be used to account for 1JCC properly, note that ¹⁵N labeling
## is always assumed in this experiment [optional, default: []]
# label = ["13C"]

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
## both are in Hz [optional, default: [] ]
# filter_offsets = [
#   [0.0, 0.0],
# ]

  ## List of the profile names and their associated filenames.
  ## The name of the spin systems should follow the Sparky-NMR
  ## conventions.
  [data.profiles]
  G2CA = "G2CA-HA.out"
  H3CA = "H3CA-HA.out"
  K4CA = "K4CA-HA.out"
  S5CA = "S5CA-HA.out"
  L6CA = "L6CA-HA.out"
```
