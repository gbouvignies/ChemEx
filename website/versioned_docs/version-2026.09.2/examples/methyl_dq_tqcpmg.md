---
sidebar_position: 2
sidebar_label: "Methyl ¹H DQ- and TQ-CPMG"
---

# Methyl ¹H DQ-CPMG and TQ-CPMG experiments

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/CPMG_CH3_1H_DQ_TQ).

:::

This example demonstrates the application of methyl ¹H DQ-CPMG
([`cpmg_ch3_1h_dq`](experiments/cpmg/cpmg_ch3_1h_dq.md)) and TQ-CPMG
([`cpmg_ch3_1h_tq`](experiments/cpmg/cpmg_ch3_1h_tq.md)) experiments for
studying <sup>13</sup>CH<sub>3</sub>-labeled methyl groups. Both experiments are carried out on the
same sample at single magnetic field. The major advantage of combining these two
types of experiments is that accurate exchange parameters can already be
obtained at single magnetic field. Typically it is necessary to carry out CPMG
experiments at multiple static fields to obtain accurate chemical exchange
parameters (p<sub>B</sub> and k<sub>ex</sub>), while in DQ- and TQ-CPMG
experiments the effective field strength behaves as two and three times of the
static field, respectively. If datasets from ¹H SQ-CPMG experiment
([`cpmg_ch3_1h_sq`](experiments/cpmg/cpmg_ch3_1h_sq.md)) are available it is
also possible to perform joint fits for all of these three types of experiments.

Both DQ-CPMG [^1] and TQ-CPMG [^2] experiments collect signals from I = 3/2
manifold, the selection of different terms of coherences is mainly achieved by
different types of phase cycling. In both experiments CPMG pulses are applied
with XY-16 phase cycle to help overcome off-resonance effects, besides, CPMG
pulses have the option to be applied as 90x180y90x composite variety, which is
indicated by `comp180_flg`.

[^1]:
    A. Gopalan, T. Yuwen, L. E. Kay, and P. Vallurupalli, _A Methyl ¹H Double
    Quantum CPMG Experiment to Study Protein Conformational Exchange_, _J Biomol
    NMR_, 72:79-91, 2018. https://doi.org/10.1007/s10858-018-0208-z

[^2]:
    T. Yuwen, P. Vallurupalli, and L. E. Kay, _Enhancing the Sensitivity of CPMG
    Relaxation Dispersion to Conformational Exchange Processes by
    Multiple-Quantum Spectroscopy_, _Angew. Chem. Int. Ed._,
    55(38):11490-11494, 2016. https://doi.org/10.1002/anie.201605843
