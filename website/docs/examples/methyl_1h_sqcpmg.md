---
sidebar_position: 1
sidebar_label: "Methyl ¹H SQ-CPMG"
---

# Methyl ¹H SQ-CPMG experiment

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_1H_SQ).

:::

This example demonstrates the application of methyl ¹H SQ-CPMG experiment
([`cpmg_ch3_1h_sq`](experiments/cpmg/cpmg_ch3_1h_sq.md)) for studying
<sup>13</sup>CH<sub>3</sub>-labeled methyl groups. It allows measuring ¹H Δϖ while still preserving
the methyl-TROSY effects.[^1] The sensitivity is much higher than older
experiments for similar purposes.[^2] This experiment can be combined with
MQ-CPMG experiment ([`cpmg_ch3_mq`](experiments/cpmg/cpmg_ch3_mq.md)) such that
Δϖ for both ¹³C and ¹H can be obtained simultaneously.

In this experiment signals from both I = 3/2 and I = 1/2 manifolds are
collected, to get rid of potential artifacts due to off-resonance effects CPMG
pulses are applied with XY-4 phase cycle. Besides, each CPMG pulse can be
applied as 90x240y90x composite variety to achieve even better off-resonance
performance, which is indicated by `comp180_flg`. To get flat relaxation
dispersion profiles in the absence of a chemical exchange process, relaxation
compensation pulses should be applied and `ncyc_max` should be set properly.

In this example, a selective set of residues are used for obtaining global
parameters (p<sub>B</sub> and k<sub>ex</sub>) first, and then residue-specific
parameters for each residue are obtained with global parameters fixed.

[^1]:
    T. Yuwen, R. Huang, P. Vallurupalli, and L. E. Kay. A Methyl-TROSY-Based ¹H
    Relaxation Dispersion Experiment for Studies of Conformational Exchange in
    High Molecular Weight Proteins _Angew. Chem. Int. Ed._ **58**, 6250-6254
    (2019). https://doi.org/10.1002/anie.201900241

[^2]:
    V. Tugarinov2007, and L. E. Kay. Separating Degenerate ¹H Transitions in
    Methyl Group Probes for Single-Quantum 1H-CPMG Relaxation Dispersion NMR
    Spectroscopy _J. Am. Chem. Soc._ **129**, 9514-9521 (2007).
    https://doi.org/10.1021/ja0726456
