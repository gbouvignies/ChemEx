# 13C (Methyl) H-to-C CPMG

Measures methyl carbon chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization is
initially anti-phase and is read out as in-phase. Because of the P-element
only even ncyc should be recorded. Resulting magnetization intensity after
the CPMG block is calculated using the (6n)×(6n), two-spin matrix, where n
is the number of states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

Lundström, Vallurupalli, Religa, Dahlquist and Kay. J Biomol NMR (2007) 38, 79-88

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_ch3_13c_h2c
