# 1HN Pure Anti-phase Proton CPMG with [0013] Phase Cycle

Analyzes amide proton chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use
the (6n)×(6n), two-spin matrix, where n is the number of states:

    { Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),
      Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b), ... }

This version is modified such that CPMG pulses are applied with [0013]
phase cycle in order to help better overcome off-resonance effects.

## References

  - Yuwen and Kay. J Biomol NMR (2019) 73:641-650

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_1hn_ap_0013
