---
sidebar_label: Pure anti-phase carbonyl ¹³C
sidebar_position: 9
description: '"cpmg_13co_ap"'
---

# Pure anti-phase carbonyl ¹³C CPMG

## Module name

`"cpmg_13co_ap"`

## Description

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a
(6n)×(6n), two-spin matrix, where n is the number of states:

    \{ COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
      COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b), ... \}

Because of the length of the shaped pulses used during the CPMG blocks,
off-resonance effects are taken into account only for the 90-degree pulses that
create COxNz before the CPMG and COzNz after the CPMG.

The calculation can be run with or without C–C J-coupling refocusing element via
the "refocusing" flag, with such option ncyc_cp should be set as even.

# Reference:

-   Lundström, D.F. Hansen, and L.E. Kay. _J. Biomol. NMR_ **42**, 35-47 (2008)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_13CO_AP/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_13co_ap'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_13co_ap"

## CPMG relaxation delay, in seconds
time_t2 = 0.03

## Position of the ¹³C carrier during the CPMG period, in ppm
carrier = 176.0

## ¹³C 90 degree pulse width of CPMG pulses, in seconds
pw90 = 25.0e-6

## Equilibration delay at the end of the CPMG period, in seconds
## [optional, default: 0.0]
# time_equil = 0.0

## Set to true when C-C J-coupling refocusing element is applied
## [optional, default: false]
# refocusing = false

## 1/2*JCC, in seconds [optional, default: 9.09e-3]
# taucc = 9.09e-3

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
  M1CO-G2N = "G2N-HN.out"
  G2CO-H3N = "H3N-HN.out"
  H3CO-K4N = "K4N-HN.out"
  K4CO-S5N = "S5N-HN.out"
  S5CO-L6N = "L6N-HN.out"
```
