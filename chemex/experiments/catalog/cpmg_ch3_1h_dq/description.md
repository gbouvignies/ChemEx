# 1H (Methyl-13CH3) Double-Quantum Proton CPMG

Measures methyl proton chemical exchange recorded on site-specifically
13CH3-labeled proteins in a highly deuterated background. Magnetization
is initially anti-phase and is read out as anti-phase. The calculation
uses a simplified scheme that includes only (6n)×(6n) basis set,
two-spin matrix, where n is the number of states:

    { 2HDQxCz(a), 2HDQyCz(a), 2HzCz(a), HDQx(a), HDQy(a), Hz(a),
      2HDQxCz(b), 2HDQyCz(b), 2HzCz(b), HDQx(b), HDQy(b), Hz(b), ... }

## References

  - Gopalan, Yuwen, Kay and Vallurupalli. J Biomol NMR (2018) 72:79-91

## Note

A sample configuration  file for this module is available using the command:

    $ chemex config cpmg_ch3_1h_dq
