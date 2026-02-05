---
sidebar_label: Pure in-phase ¹⁵N D-CEST
sidebar_position: 1
description: '"dcest_15n"'
---

# Pure in-phase ¹⁵N DANTE-CEST

## Module name

`"dcest_15n"`

## Description

Analyzes chemical exchange in the presence of ¹H composite decoupling during the
D-CEST block. This keeps the spin system purely in-phase throughout, and is
calculated using the (3n)×(3n), single-spin matrix, where n is the number of
states:

    \{ Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... \}

## References

-   T. Yuwen, L.E. Kay, and G. Bouvignies. _ChemPhysChem_ **19**, 1707-1710 (2018)
-   T. Yuwen, A. Bah, Brady, F. Ferrage, G. Bouvignies, and L.E. Kay. _J. Phys.
    Chem. B_ **122**, 11206-11217 (2018)

## Examples

-   An example for studying a two-state exchange system is given
    [here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N/).
-   An example for studying a three-state exchange system is given
    [here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N_3States/).
-   An example for H/D solvent exchange measurement using ¹⁵N D-CEST can be found
    [here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N_HD_EXCH/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'dcest_15n'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "dcest_15n"

## CEST relaxation delay, in seconds
time_t1 = 0.5

## Position of the ¹⁵N carrier during the CEST period, in ppm
carrier = 118.0

## Pulse width of a 90 pulse at the power used during the DANTE, in seconds
pw90 = 45e-6

## DANTE "spectral width", in Hz
sw = 800.0

## Effective B1 field for the equivalent continuous-wave CEST irradiation, in Hz
b1_eff = 25.0  # alias: b1_frq (still supported)

## Notes on pw90 vs b1_eff (D-CEST specific):
## - pw90 defines the hardware B1 amplitude; the B1 inhomogeneity distribution
##   is centered on the nominal B1 derived from pw90 (1/(4*pw90)).
## - b1_eff sets the effective B1 used to time the DANTE pulses (pw_dante).
## You must provide both for D-CEST: pw90 (hardware) and b1_eff (effective).
## Backward-compatibility: the key 'b1_frq' is accepted as an alias for 'b1_eff'.

## D-CEST timing relationships (from Supporting Information):
## Let T be the total D-CEST duration, k the number of pulses, τp the pulse
## width, and τ' the inter-pulse interval. With spectral width sw = 1/τ',
## Eq. [S12] gives: τp = (b1_eff / b1_D-CEST) * τ'.

## B1 inhomogeneity distribution (replaces b1_inh_scale/b1_inh_res)
## type:  distribution family (e.g., "gaussian")
## scale: relative standard deviation (e.g., 0.1 means 10% width)
## res:   number of quadrature points used in the distribution
[experiment.b1_distribution]
type = "gaussian"
scale = 0.1
res = 11

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
  G2N = "G2N-HN.out"
  H3N = "H3N-HN.out"
  K4N = "K4N-HN.out"
  S5N = "S5N-HN.out"
  L6N = "L6N-HN.out"
```
