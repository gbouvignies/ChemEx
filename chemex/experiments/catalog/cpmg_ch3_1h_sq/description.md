# 1H (Methyl-13CH3) Single-Quantum Proton CPMG

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as anti-phase prior to 1H detection.
Resulting magnetization intensity after the CPMG block is calculated using
the (6n)×(6n), two-spin matrix, where n is the number of states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  - Yuwen, Huang, Vallurupalli and Kay. Angew Chem Int Ed (2019) 58:6250-6254

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_ch3_1h_sq
