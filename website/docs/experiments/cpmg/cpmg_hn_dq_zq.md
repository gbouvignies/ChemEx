---
sidebar_label: Amide ¹⁵N–¹H double-quantum/zero-quantum
sidebar_position: 7
description: '"cpmg_hn_dq_zq"'
---

# Amide ¹⁵N–¹H double-quantum/zero-quantum CPMG

## Module name

`"cpmg_hn_dq_zq"`

## Description

Analyzes ¹⁵N and ¹H chemical exchange by applying CPMG pulses on ¹⁵N and 1H
simultaneously. The spin system is maintained as DQ or ZQ during Trelax, and is
calculated using the (15n)×(15n), two-spin matrix, where n is the number of
states:

    \{        Ix(a),   Iy(a),   Iz(a),   Sx(a), IxSx(a), IySx(a), IzSx(a),
      Sy(a), IxSy(a), IySy(a), IzSy(a), Sz(a), IxSz(a), IySz(a), IzSz(a),
             Ix(b),   Iy(b),   Iz(b),   Sx(b), IxSx(b), IySx(b), IzSx(b),
      Sy(b), IxSy(b), IySy(b), IzSy(b), Sz(b), IxSz(b), IySz(b), IzSz(b), ... \}

The phase cycle of CPMG pulses is chosen based on νCPMG as described in the
reference, which is a mixture of constant-phase and XY-family phase cycles.

## References

-   Orekhov, Korzhnev, and L.E. Kay _J. Am. Chem. Soc._ **126**, 1886-1891 (2004)

## Sample configuration file

```toml title="experiment.toml"
## This is a sample configuration file for the module 'cpmg_hn_dq_zq'

[experiment]

## Name of the chemex module corresponding to the experiment
name = "cpmg_hn_dq_zq"

## CPMG relaxation delay, in seconds
time_t2 = 0.04

## Position of the ¹⁵N carrier during the CPMG period, in ppm
carrier_n = 118.0

## Position of the 1H carrier during the CPMG period, in ppm
carrier_h = 8.3

## ¹⁵N 90 degree pulse width of CPMG pulses, in seconds
pw90_n = 40.0e-6

## 1H 90 degree pulse width of CPMG pulses, in seconds
pw90_h = 15.0e-6

## Perform DQ CPMG RD experiment, otherwise perform ZQ CPMG RD experiment
dq_flg = true

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
