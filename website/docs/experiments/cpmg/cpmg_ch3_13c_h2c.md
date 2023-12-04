---
sidebar_label: Methyl ¹³C (“H-to-C”)
sidebar_position: 10
description: '"cpmg_ch3_13c_h2c"'
---

# Methyl ¹³C CPMG (“H-to-C”)

## Module name

`"cpmg_ch3_13c_h2c"`

## Description

Measures methyl carbon chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as in-phase. Because of the P-element only
even ncyc should be recorded. Resulting magnetization intensity after the CPMG
block is calculated using the (6n)×(6n), two-spin matrix, where n is the number
of states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

## Reference

-   P. Lundström, P. Vallurupalli, T.M. Religa, F.W. Dahlquist, and L.E. Kay. _J.
    Biomol. NMR_ **38**, 79-88 (2007)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_13C_H2C/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_ch3_13c_h2c'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_ch3_13c_h2c"

## CPMG relaxation delay, in seconds
time_t2 = 0.04

## Position of the ¹³C carrier during the CPMG period, in ppm
carrier = 20.0

## ¹³C 90 degree pulse width of CPMG pulses, in seconds
pw90 = 15.0e-6

## Equilibration delay at the end of the CPMG period, in seconds
## [optional, default: 0.0]
# time_equil = 0.0

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
