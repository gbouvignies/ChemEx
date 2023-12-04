---
sidebar_label: Pure anti-phase amide ¹H with [0013] phase cycle
sidebar_position: 6
description: '"cpmg_1hn_ap_0013"'
---

# Pure anti-phase amide ¹H CPMG with [0013] phase cycle

## Module name

`"cpmg_1hn_ap_0013"`

## Description

Analyzes amide proton chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use the
(6n)×(6n), two-spin matrix, where n is the number of states:

    \{ Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),
      Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b), ... \}

This version is modified such that CPMG pulses are applied with [0013] phase
cycle to help better overcome off-resonance effects.

## Reference

-   T. Yuwen, and L.E. Kay. _J. Biomol. NMR_ **73**, 641-650 (2019)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_1HN_AP_0013/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_1hn_ap_0013'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_1hn_ap_0013"

## CPMG relaxation delay, in seconds
time_t2 = 0.02

## Position of the 1H carrier during the CPMG period, in ppm
carrier = 8.3

## 1H 90 degree pulse width of CPMG pulses, in seconds
pw90 = 10.0e-6

## Maximum number of cycles
ncyc_max = 40

## 1/4J, in seconds [optional, default: 2.38e-3]
# taua = 2.38e-3

## Apply IPAP scheme for IP/AP differential relaxation suppression
## [optional, default: false]
# ipap_flg = false

## Use the EBURP scheme instead of the standard 180 central pulse
## [optional, default: false]
# eburp_flg = false

## Use the REBURP scheme instead of the standard 180 central pulse
## [optional, default: false]
# reburp_flg = false

## 1H EBURP pulse width of CPMG pulses, in seconds [optional, default: 1.4e-3]
# pw_eburp = 1.4e-3

## 1H REBURP pulse width of CPMG pulses, in seconds [optional, default: 1.52e-3]
# pw_reburp = 1.52e-3

## Equilibration delay at the beginning of the CPMG period, in seconds
## [optional, default: 0.0]
# time_equil_1 = 0.0

## Equilibration delay at the end of the CPMG period, in seconds
## [optional, default: 0.0]
# time_equil_2 = 0.0

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
  G2HN-N = "G2N-HN.out"
  H3HN-N = "H3N-HN.out"
  K4HN-N = "K4N-HN.out"
  S5HN-N = "S5N-HN.out"
  L6HN-N = "L6N-HN.out"
```
