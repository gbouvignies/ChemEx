---
sidebar_label: In-phase/anti-phase methyl ¹H
sidebar_position: 7
description: '"cest_ch3_1h_ip_ap"'
---

# In-phase/anti-phase methyl ¹H CEST

## Module name

`"cest_ch3_1h_ip_ap"`

## Description

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

:::note

The same module works for both <sup>13</sup>CH<sub>3</sub>- and <sup>13</sup>CHD<sub>2</sub>-labeled methyl groups.

:::

## References

-   T. Yuwen, and L.E. Kay. _J. Biomol. NMR_ **68**, 215-224 (2017)
-   T. Yuwen, and L.E. Kay. _J. Biomol. NMR_ **70**, 93-102 (2018)

## Example

An example recorded on a <sup>13</sup>CH<sub>3</sub>-labeled sample is available
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_CH3_1H_IP_AP/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cest_ch3_1h_ip_ap'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cest_ch3_1h_ip_ap"

## Recycle delay, in seconds
d1 = 1.0

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the carrier during the CEST period, in ppm
carrier = 0.5

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
  L3HD1-CD1  = "L3CD1-QD1.out"
  L3HD2-CD2  = "L3CD2-QD2.out"
  I12HD1-CD1 = "I12CD1-QD1.out"
  V25HG1-CG1 = "V25CG1-QG1.out"
  V25HG2-CG2 = "V25CG2-QG2.out"
  M36HE-CE   = "M36CE-QE.out"
```
