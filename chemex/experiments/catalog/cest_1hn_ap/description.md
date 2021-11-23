# 1HN Pure Anti-phase Proton CEST

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  -Sekhar, Rosenzweig, Bouvignies and Kay. PNAS (2016) 113:E2794-E2801


## Note

A sample configuration file for this module is available using the command:

    $ chemex config cest_1hn_ap
