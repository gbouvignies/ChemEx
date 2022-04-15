# 1H (Methyl-13CH3/13CHD2) In-phase/Anti-phase Proton CEST

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  - Yuwen and Kay. J Biomol NMR (2017) 68:215-224
  - Yuwen and Kay. J Biomol NMR (2018) 70:93-102

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cest_ch3_1h_ip_ap
