# 1HN In-phase/Anti-phase Proton CEST

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  - Yuwen, Sekhar and Kay. Angew Chem Int Ed (2017) 56:6122-6125
  - Yuwen and Kay. J Biomol NMR (2017) 67:295-307
  - Yuwen and Kay. J Biomol NMR (2018) 70:93-102


## Note

A sample configuration file for this module is available using the command:

    $ chemex config cest_1hn_ip_ap
