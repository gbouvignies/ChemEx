# 15N–1HN TROSY CPMG with [0013] Phase Cycle

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
      Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b), ... }

This version is modified such that CPMG pulses are applied with [0013] phase
cycle as used in 15N pure in-phase experiments.

## References

  - Jiang, Yu, Zhang, Liu and Yang. J Magn Reson (2015) 257:1-7
  - Vallurupalli, Hansen, Stollar, Meirovitch and Kay. PNAS (2007) 104:18473-18477


## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_15n_tr_0013
