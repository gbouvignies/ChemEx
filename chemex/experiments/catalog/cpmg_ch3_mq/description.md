# 13C–1H (Methyl-13CH3) Multiple-Quantum CPMG

Analyzes HyCx methyl group multiple quantum CPMG measured on site-specific
13CH3-labeled methyl groups in a highly deuterated background. Resulting
magnetization intensity after the CPMG block is calculated using the
(4n)×(4n), two-spin matrix, where n is the number of states:

    { HxCx(a), HyCx(a), HxCy(a), HyCy(a),
      HxCx(b), HyCx(b), HxCy(b), HyCy(b), ... }

This is a simplified basis set, which assumes 13C is on-resonance
(i.e., off-resonance effects are not taken into account).

This calculation is designed specifically to analyze data from the experiment
found in the reference and can be run with either small_protein_flag = "y"
or "n".

## References

  - Korzhnev, Kloiber, Kanelis, Tugarinov and Kay. J Am Chem Soc (2004) 126:3964-3973

## Note

A sample configuration  file for this module is available using the command:

    $ chemex config cpmg_ch3_mq
