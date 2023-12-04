---
sidebar_label: Methyl ¹H double quantum
sidebar_position: 13
description: '"cpmg_ch3_1h_dq"'
---

# Methyl ¹H double quantum CPMG

## Module name

`"cpmg_ch3_1h_dq"`

## Description

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as anti-phase. The calculation uses a
simplified scheme that includes only (6n)×(6n) basis set, two-spin matrix, where
n is the number of states:

    \{ 2HDQxCz(a), 2HDQyCz(a), 2HzCz(a), HDQx(a), HDQy(a), Hz(a),
      2HDQxCz(b), 2HDQyCz(b), 2HzCz(b), HDQx(b), HDQy(b), Hz(b), ... \}

## Reference

-   Gopalan, T. Yuwen, L.E. Kay, and P. Vallurupalli. _J. Biomol. NMR_ **72**,
    79-91 (2018)

## Example

An example, where joint fit of
[methyl ¹H double quantum CPMG (cpmg_ch3_1h_dq)](cpmg_ch3_1h_dq.md) and
[methyl ¹H triple quantum CPMG (cpmg_ch3_1h_tq)](cpmg_ch3_1h_tq.md) experiments
is performed, is available
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/CPMG_CH3_1H_DQ_TQ/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_ch3_1h_dq'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_ch3_1h_dq"

## CPMG relaxation delay, in seconds
time_t2 = 0.04

## Position of the 1H carrier during the CPMG period, in ppm
carrier = 0.5

## 1H 90 degree pulse width of CPMG pulses, in seconds
pw90 = 12.0e-6

## Performs CPMG using 90-180-90 composite pulses [optional, default: true]
# comp180_flg = true

## Enters the CPMG period with equal amount of initial IP and
## AP DQ magnetizations [optional, default: false]
# ipap_flg = false

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
  L3HD1-CD1  = "L3CD1-QD1.out"
  L3HD2-CD2  = "L3CD2-QD2.out"
  I12HD1-CD1 = "I12CD1-QD1.out"
  V25HG1-CG1 = "V25CG1-QG1.out"
  V25HG2-CG2 = "V25CG2-QG2.out"
  M36HE-CE   = "M36CE-QE.out"
```
