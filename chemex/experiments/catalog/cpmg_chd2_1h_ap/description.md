# 1H (Methyl-13CHD2) Pure Anti-Phase Proton CPMG

Measures methyl proton chemical exchange recorded on site-specifically
13CHD2-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as anti-phase prior to 13C evolution.
Resulting magnetization intensity after the CPMG block is calculated using
the (6n)×(6n), two-spin matrix, where n is the number of states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  - Baldwin, Religa, Hansen, Bouvignies and Kay. J Am Chem Soc (2010) 132:10992-10995

## Note

A sample configuration  file for this module is available using the command:

    $ chemex config cpmg_chd2_1h_ap
