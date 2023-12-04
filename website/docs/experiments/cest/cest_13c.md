---
sidebar_label: Pure in-phase ¹³C
sidebar_position: 2
description: '"cest_13c"'
---

# Pure in-phase ¹³C CEST

## Module name

`"cest_13c"`

## Description

Analyzes chemical exchange in the presence of ¹H composite decoupling during the
CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)×(3n), single-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... \}

## References

-   P. Vallurupalli, G. Bouvignies, and L.E. Kay. _ChemBioChem_ **14**, 1709-1713
    (2014)
-   G. Bouvignies, P. Vallurupalli, and L.E. Kay. _J. Mol. Biol._ **426**, 763-774
    (2014)
-   P. Vallurupalli, and L.E. Kay. _Angew. Chem. Int. Ed._ **52**, 4156-4159
    (2013)
-   D.F. Hansen, G. Bouvignies, and L.E. Kay. _J. Biomol. NMR_ **55**, 279-289
    (2013)
-   G. Bouvignies, and L.E. Kay. _J. Biomol. NMR_ **53**, 303-310 (2012)
-   E. Rennella, R. Huang, A. Velyvis, and L.E. Kay. _J. Biomol. NMR_ **63**,
    187-199 (2015)

## Example

-   An example for studying side-chain methyl groups in selectively ¹³C-labeled
    sample is available
    [here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_13C/).
-   An example for studying
    [uniformly ¹³C, ¹⁵N-labeled sample](../../examples/cest_13c_15n.md) is
    available
    [here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_13C_LABEL_CN/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cest_13c'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cest_13c"

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the ¹³C carrier during the CEST period, in ppm
carrier = 45.0

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
## both are in Hz [optional, default: [[0.0, 0.0]] ]
# filter_offsets = [
#   [0.0, 0.0],
# ]

  ## List of the profile names and their associated filenames.
  ## The name of the spin systems should follow the Sparky-NMR
  ## conventions.
  [data.profiles]
  G2CA  = ["G2CA-HA1.out", "G2CA-HA2.out"]
  H3CB  = ["H3CB-HB2.out", "H3CB-HB3.out"]
  K4CD  = ["K4CD-HD2.out", "K4CD-HD3.out"]
  S5CB  = ["S5CB-HB2.out", "S5CB-HB3.out"]
  L6CA  = "L6CA-HA.out"
  L6CD1 = "L6CD1-QD1.out"
```
