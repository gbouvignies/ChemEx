# 15N–1HN TROSY CPMG

Analyzes 15N constant-time TROSY CPMG relaxation dispersion experiments for
measurement of ΔD NH in protein systems undergoing millisecond-time-scale
exchange dynamics. Resulting magnetization intensity after the CPMG block is
calculated using the (6n)×(6n), two-spin matrix, where n is the number of
states:

    { Nx(a), Ny(a), Nz(a), 2HzNx(a), 2HzNy(a), 2HzNz(a),
      Nx(b), Ny(b), Nz(b), 2HzNx(b), 2HzNy(b), 2HzNz(b), ... }

## References

  - Vallurupalli, Hansen, Stollar, Meirovitch and Kay. PNAS (2007) 104:18473-18477

## Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_15n_tr
