---
sidebar_label: Methyl ¹³C (“H-to-C”) with [0013] phase cycle
sidebar_position: 17
description: '"cpmg_ch3_13c_h2c_0013"'
---

# Methyl ¹³C CPMG (“H-to-C”) with [0013] phase cycle

## Module name

`"cpmg_ch3_13c_h2c_0013"`

## Description

Measures methyl carbon chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as in-phase. Resulting magnetization
intensity after the CPMG block is calculated using the (6n)×(6n), two-spin
matrix, where n is the number of states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

This version is modified such that CPMG pulses are applied with [0013] phase
cycle, and additional purge gradients are applied around the P-element, which
help to overcome off-resonance effects, variations of 1JCH coupling constants,
and pulse miscalibration. In this version ncyc can be set as even or odd (for
odd ncyc the numbers of CPMG pulses on both sides of the P-element are unequal).

## Reference

-   T. Yuwen, J. Liu, Z. Xia, Y. Xia, P. Rossi, and C. G. Kalodimos. _J. Biomol.
    NMR_ **79** (2007) https://doi.org/10.1007/s10858-025-00474-x

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_13C_H2C_0013/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_ch3_13c_h2c_0013'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_ch3_13c_h2c_0013"

## CPMG relaxation delay, in seconds
time_t2 = 0.04

## Position of the ¹³C carrier during the CPMG period, in ppm
carrier = 20.0

## ¹³C 90 degree pulse width of CPMG pulses, in seconds
pw90 = 15.0e-6

## Equilibration delay at the end of the CPMG period, in seconds
## [optional, default: 5.0e-3]
# time_equil = 0.0

## Sum of purge gradient and recovery delay on each side of the P-element, in seconds
## [optional, default: 1.2e-3]
# time_grad = 0.0

## P-element delay = 1/4J, in seconds [optional, default: 2.0e-3]
# taub = 2.0e-3

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
  L3CD1-HD1  = "L3CD1-QD1.out"
  L3CD2-HD2  = "L3CD2-QD2.out"
  I12CD1-HD1 = "I12CD1-QD1.out"
  V25CG1-HG1 = "V25CG1-QG1.out"
  V25CG2-HG2 = "V25CG2-QG2.out"
  M36CE-HE   = "M36CE-QE.out"
```
