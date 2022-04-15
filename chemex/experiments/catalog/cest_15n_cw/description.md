# 15N CEST with CW Decoupling

Analyzes chemical exchange in the presence of 1H CW decoupling during the
CEST block. Magnetization evolution is calculated using the (15n)×(15n),
two-spin matrix, where n is the number of states::

    {        Ix(a),   Iy(a),   Iz(a),   Sx(a), IxSx(a), IySx(a), IzSx(a),
      Sy(a), IxSy(a), IySy(a), IzSy(a), Sz(a), IxSz(a), IySz(a), IzSz(a),
             Ix(b),   Iy(b),   Iz(b),   Sx(b), IxSx(b), IySx(b), IzSx(b),
      Sy(b), IxSy(b), IySy(b), IzSy(b), Sz(b), IxSz(b), IySz(b), IzSz(b), ... }

## References

  - Bouvignies and Kay. J Phys Chem B (2012) 116:14311-14317


## Note

A sample configuration file for this module is available using the command::

    $ chemex config cest_15n_cw
