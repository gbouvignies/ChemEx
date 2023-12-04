---
sidebar_label: Pure in-phase ¹³C
sidebar_position: 8
description: '"cpmg_13c_ip"'
---

# Pure in-phase ¹³C CPMG

## Module name

`"cpmg_13c_ip"`

## Description

Analyzes ¹³C chemical exchange in the presence of high power ¹H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)×(3n), single-spin matrix, where n is the number
of states:

    \{ Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... \}

The CW decoupling on ¹H is assumed to be strong enough (> 15 kHz) such that
perfect ¹H decoupling can be achieved. In the case of CHD2 experiment, CW
decoupling on 2H should be applied properly.

## References

-   D.F. Hansen, P. Vallurupalli, Lundström, Neudecker, and L.E. Kay. _J. Am.
    Chem. Soc._ **130**, 2667-2675 (2008)
-   E. Rennella, Schuetz, and L.E. Kay. _J. Biomol. NMR_ **65**, 59-64 (2016)

## Example

An example use of the module, based on datasets measured for CHD2-labeled
side-chain methyl groups, is given
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_13C_IP/).

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_13c_ip'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_13c_ip"

## CPMG relaxation delay, in seconds
time_t2 = 0.04

## Position of the ¹³C carrier during the CPMG period, in ppm
carrier = 18.0

## ¹³C 90 degree pulse width of CPMG pulses, in seconds
pw90 = 15.0e-6

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
  L3CD1  = "L3CD1-QD1.out"
  L3CD2  = "L3CD2-QD2.out"
  I12CD1 = "I12CD1-QD1.out"
  V25CG1 = "V25CG1-QG1.out"
  V25CG2 = "V25CG2-QG2.out"
  M36CE  = "M36CE-QE.out"
```
