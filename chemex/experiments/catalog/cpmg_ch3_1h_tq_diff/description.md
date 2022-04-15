# 1H (Methyl-13CH3) Triple-Quantum Proton CPMG with Gradients

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. The experiment
is specifically designed for measuring translational diffusion constants
of invisible states using a pulsed-field gradient approach that exploits
methyl 1H triple-quantum coherences. Note that the diffusion constants
are given in μm²/s. Magnetization is initially anti-phase and is read out
as anti-phase. The calculation uses a simplified scheme that includes only
(6n)×(6n) basis set, two-spin matrix, where n is the number of states:

    { 2HTQxCz(a), 2HTQyCz(a), 2HzCz(a), HTQx(a), HTQy(a), Hz(a),
      2HTQxCz(b), 2HTQyCz(b), 2HzCz(b), HTQx(b), HTQy(b), Hz(b), ... }

## References

  - Yuwen, Sekhar, Baldwin, Vallurupalli and Kay. Angew Chem Int Ed (2018) 57, 16777-16780

## Note

A sample configuration  file for this module is available using the command:

    $ chemex config cpmg_ch3_1h_tq_diff
