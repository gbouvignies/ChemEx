# 13C (Carbonyl) Pure Anti-phase CPMG

Analyzes carbonyl chemical exchange that is maintained as anti-phase
magnetization throughout the CPMG block. This results in lower intrinsic
relaxation rates and therefore better sensitivity. The calculations use a
(6n)×(6n), two-spin matrix, where n is the number of states:

    { COx(a), COy(a), COz(a), 2COxNz(a), 2COyNz(a), 2COzNz(a),
      COx(b), COy(b), COz(b), 2COxNz(b), 2COyNz(b), 2COzNz(b), ... }

Because of the length of the shaped pulses used during the CPMG blocks,
off-resonance effects are taken into account only for the 90-degree pulses
that create COxNz before the CPMG and COzNz after the CPMG.

The calculation can be run with or without C–C J-coupling refocusing element
via the "refocusing" flag, with such option ncyc_cp should be set as even.

# References

  - Lundström, Hansen and Kay. J Biomol NMR (2008) 42:35-47

# Note

A sample configuration file for this module is available using the command:

    $ chemex config cpmg_13co_ap
