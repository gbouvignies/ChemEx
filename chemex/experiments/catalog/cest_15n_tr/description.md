# 15N–1HN TROSY CEST

Analyzes chemical exchange during the CEST block. Magnetization evolution is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),
      Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

Long, Bouvignies and Kay. PNAS (2014) 111:8820-8825

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cest_15N_tr
