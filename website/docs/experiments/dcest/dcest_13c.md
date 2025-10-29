---
sidebar_label: Pure in-phase ¹³C D-CEST
sidebar_position: 2
description: '"dcest_13c"'
---

# Pure in-phase ¹³C DANTE-CEST

## Module name

`"dcest_13c"`

## Description

Analyzes chemical exchange in the presence of ¹H composite decoupling during the
D-CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)×(3n), single-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... \}

## References

-   T. Yuwen, L.E. Kay, and G. Bouvignies. _ChemPhysChem_ **19**, 1707-1710 (2018)
-   T. Yuwen, G. Bouvignies, and L.E. Kay. _J. Mag. Reson._ **292**, 1-7 (2018)

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'dcest_13c'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "dcest_13c"

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the ¹³C carrier during the CEST period, in ppm
carrier = 176.0

## Pulse width of a 90 pulse at the power used during the DANTE, in seconds
pw90 = 15e-6

## DANTE "spectral width", in Hz
sw = 800.0

## Effective B1 field for the equivalent continuous-wave CEST irradiation, in Hz
b1_eff = 25.0  # alias: b1_frq (still supported)

## B1 inhomogeneity distribution (replaces b1_inh_scale/b1_inh_res)
[experiment.b1_distribution]
type = "gaussian"  # distribution shape
scale = 0.1        # fractional standard deviation relative to b1_frq
res = 11           # number of sampling points across the distribution

## Notes on pw90 vs b1_eff (D-CEST specific):
## - The B1 distribution is centered on the nominal B1 implied by pw90 (1/(4*pw90)).
## - b1_eff (when provided) is used to compute the DANTE pulse width (pw_dante).
## You must provide both for D-CEST. 'b1_frq' remains an alias for b1_eff.

## Notes on pw90 vs b1_frq (D-CEST specific):
## - The B1 distribution is centered on the nominal B1 implied by pw90 (1/(4*pw90)).
## - b1_frq (when provided) is used to compute the DANTE pulse width (pw_dante).

## Equilibration delay at the end of the CEST period, in seconds
## [optional, default: 0.0]
# time_equil = 0.0

## Initial condition: equilibrium (all states populated according to their
## equilibrium populations). Set to true for non-equilibrium initial condition
## if chemical shift evolution occurs before the CEST element in your pulse
## sequence (see Yuwen et al., J. Biomol. NMR 2016, 65:143-156).
## [optional, default: false]
# cs_evolution_prior = false

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

## The amount of D2O in the sample, in mass fraction
## [optional, only used for HD exchange measurement]
# d2o = 0.1

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
  A1C = "G2N-HN.out"
  G2C = "H3N-HN.out"
  H3C = "K4N-HN.out"
  K4C = "S5N-HN.out"
  S5C = "L6N-HN.out"
```
