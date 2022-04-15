# 1HN Pure Anti-phase Proton CPMG

Analyzes amide proton chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use
a (6n)×(6n), two-spin matrix, where n is the number of states:

    { Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),
      Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b), ... }

## References

Adapted from:

  - Ishima and Torshia. J Biomol NMR (2003) 25:243-248

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_1hn_ap
