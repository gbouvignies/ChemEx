---
sidebar_label: Pure anti-phase methyl ¹H
sidebar_position: 16
description: '"cpmg_chd2_1h_ap"'
---

# Pure anti-phase methyl ¹H CPMG

## Module name

`"cpmg_chd2_1h_ap"`

## Description

Measures methyl proton chemical exchange recorded on site-specifically
13CHD2-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as anti-phase prior to ¹³C evolution.
Resulting magnetization intensity after the CPMG block is calculated using the
(6n)×(6n), two-spin matrix, where n is the number of states:

    \{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... \}

## Reference

-   A.J. Baldwin, T.M. Religa, D.F. Hansen, G. Bouvignies, and L.E. Kay. _J. Am.
    Chem. Soc._ **132**, 10992-10995 (2010)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CHD2_1H_AP/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_chd2_1h_ap'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_chd2_1h_ap"

## CPMG relaxation delay, in seconds
time_t2 = 0.02

## Position of the 1H carrier during the CPMG period, in ppm
carrier = 0.5

## 1H 90 degree pulse width of CPMG pulses, in seconds
pw90 = 10.0e-6

## Equilibration delay at the end of the CPMG period, in seconds
## [optional, default: 0.0]
# time_equil = 0.0

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
  L3HD1-CD1  = "L3CD1-QD1.out"
  L3HD2-CD2  = "L3CD2-QD2.out"
  I12HD1-CD1 = "I12CD1-QD1.out"
  V25HG1-CG1 = "V25CG1-QG1.out"
  V25HG2-CG2 = "V25CG2-QG2.out"
  M36HE-CE   = "M36CE-QE.out"
```
