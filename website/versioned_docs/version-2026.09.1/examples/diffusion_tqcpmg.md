---
sidebar_position: 3
sidebar_label: "Measuring diffusion constants of invisible protein conformers"
---

# Diffusion rate measurement for excited state with methyl ¹H TQ-CPMG

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_1H_TQ_DIFF).

:::

This example demonstrates the use of CPMG experiment to measure diffusion rate
constant of the excited state. In order to achieve such purpose an additional
pair of gradients should be placed around the CPMG element, such idea has been
successfully applied in the modified version of methyl ¹H TQ-CPMG
experiment.[^1] In order to take into account the additional pair of gradients,
a separate module
([`cpmg_ch3_1h_tq_diff`](experiments/cpmg/cpmg_ch3_1h_tq_diff.md)) is required
to fit datasets from this experiment.

In this example methyl ¹H TQ-CPMG datasets with different strengths of
additional gradients are measured to obtain diffusion rate constant of the
excited state. Note that those parameters about gradients around the CPMG
element should be specified in experiment files, such as the strength and
duration.

[^1]:
    T. Yuwen, A. Sekhar, A. Baldwin, P. Vallurupalli, and L. E. Kay. Measuring
    Diffusion Constants of Invisible Protein Conformers by Triple-Quantum ¹H
    CPMG Relaxation Dispersion. _Angew. Chem. Int. Ed._ **57**, 16777-16780
    (2018). https://doi.org/10.1002/anie.201810868
