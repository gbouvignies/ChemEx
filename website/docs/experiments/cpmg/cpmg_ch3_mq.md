---
sidebar_label: Methyl ¹³C–¹H multiple-quantum
sidebar_position: 11
description: '"cpmg_ch3_mq"'
---

# Methyl ¹³C–¹H multiple-quantum CPMG

## Module name

`"cpmg_ch3_mq"`

## Description

Analyzes HyCx methyl group multiple quantum CPMG measured on site-specific
13CH3-labeled methyl groups in a highly deuterated background. Resulting
magnetization intensity after the CPMG block is calculated using the (4n)×(4n),
two-spin matrix, where n is the number of states:

    \{ HxCx(a), HyCx(a), HxCy(a), HyCy(a),
      HxCx(b), HyCx(b), HxCy(b), HyCy(b), ... \}

This is a simplified basis set, which assumes ¹³C is on-resonance (i.e.,
off-resonance effects are not taken into account).

This calculation is designed specifically to analyze data from the experiment
found in the reference and can be run with either small_protein_flag = "y" or
"n".

## Reference

-   D.M. Korzhnev, K. Kloiber, V. Kanelis, V. Tugarinov and L.E. Kay _J. Am. Chem.
    Soc._ **126**, 3964-3973 (2004)

## Example

An example use of the module is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_MQ/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_ch3_mq'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_ch3_mq"

## CPMG relaxation delay, in seconds
time_t2 = 0.02

## Perform the calculation with the purge element [optional, default: false])
# small_protein = false

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
