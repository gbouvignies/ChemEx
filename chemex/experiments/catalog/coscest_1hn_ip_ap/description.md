# 1HN In-phase/Anti-phase Proton COS-CEST

Analyzes chemical exchange during the COS-CEST block. Magnetization evolution
is calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  - Yuwen, Bouvignies and Kay. J Mag Reson (2018) 292:1-7

## Note

A sample configuration file for this module is available using the command:

    $ chemex config coscest_1hn_ip_ap
