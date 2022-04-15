# 15N Pure In-phase CPMG

Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)×(3n), single-spin matrix, where n is the
number of states:

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

The CW decoupling on 1H is assumed to be strong enough (> 15 kHz) such that
perfect 1H decoupling can be achieved.

## References

  - Hansen, Vallurupalli and Kay. *J Phys Chem B* (2008) 112:5898-5904


## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_15n_ip
