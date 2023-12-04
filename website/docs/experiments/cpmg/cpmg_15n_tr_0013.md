---
sidebar_label: ¹⁵N TROSY with [0013] phase cycle
sidebar_position: 4
description: '"cpmg_15n_tr_0013"'
---

# ¹⁵N TROSY CPMG with [0013] phase cycle

## Module name

`"cpmg_15n_tr_0013"`

## Description

Analyzes ¹⁵N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    \{ Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
      Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b), ... \}

This version is modified such that CPMG pulses are applied with [0013] phase
cycle as used in ¹⁵N pure in-phase experiments.

## References

-   Jiang, Yu, Zhang, Liu, and Yang. J Magn Reson **257**, 1-7 (2015)
-   P. Vallurupalli, D.F. Hansen, E. Stollar, E. Meirovitch, and L.E. Kay. _Proc.
    Natl. Acad. Sci. USA_ **104**, 18473-18477 (2007)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_15N_TR_0013/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_15n_tr_0013'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_15n_tr_0013"

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

## P-element delay = 1/4J, in seconds [optional, default: 2.68e-3]
# taub = 2.68e-3

## Perform anti-trosy CPMG RD experiment [optional, default: false]
# antitrosy = false

## S3E trosy selection [optional, default: true]
# s3e = true

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
