# 13C Pure In-phase CPMG

Analyzes 13C chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the (3n)×(3n), single-spin matrix, where n is the
number of states:

    { Ix(a), Iy(a), Iz(a),
      Ix(b), Iy(b), Iz(b), ... }

The CW decoupling on 1H is assumed to be strong enough (> 15 kHz) such that
perfect 1H decoupling can be achieved. In the case of CHD2 experiment, CW
decoupling on 2H should be applied properly.

## References

  - Hansen, Vallurupalli, Lundström, Neudecker and Kay. J Am Chem Soc (2008) 130:2667-2675
  - Rennella, Schuetz and Kay. J Biomol NMR (2016) 65:59-64

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_13c_ip
