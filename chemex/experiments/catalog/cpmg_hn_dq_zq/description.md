# 15N–1HN DQ/ZQ CPMG

Analyzes 15N and 1H chemical exchange by applying CPMG pulses on 15N and 1H
simultaneously. The spin system is maintained as DQ or ZQ during Trelax, and
is calculated using the (15n)×(15n), two-spin matrix, where n is the number
of states:

    {        Ix(a),   Iy(a),   Iz(a),   Sx(a), IxSx(a), IySx(a), IzSx(a),
      Sy(a), IxSy(a), IySy(a), IzSy(a), Sz(a), IxSz(a), IySz(a), IzSz(a),
             Ix(b),   Iy(b),   Iz(b),   Sx(b), IxSx(b), IySx(b), IzSx(b),
      Sy(b), IxSy(b), IySy(b), IzSy(b), Sz(b), IxSz(b), IySz(b), IzSz(b), ... }

The phase cycle of CPMG pulses is chosen based on νCPMG as described in the
reference, which is a mixture of constant-phase and XY-family phase cycles.

## References

  - Orekhov, Korzhnev and Kay. J Am Chem Soc (2004) 126:1886-1891

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_hn_dq_zq
