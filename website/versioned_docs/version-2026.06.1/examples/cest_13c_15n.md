---
sidebar_position: 5
sidebar_label: "CEST experiments for ¹³C, ¹⁵N-labeled samples"
---

# CEST experiments for ¹³C, ¹⁵N-labeled samples

:::note Files

You can get the example files from the GitHub ChemEx page:
[¹⁵N CEST](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N_LABEL_CN)
and
[¹³C CEST](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_13C_LABEL_CN).

:::

This example demonstrates the analysis of CEST datasets measured for uniformly
¹³C, ¹⁵N-labeled samples, which include both ¹³C
([`cest_13c`](experiments/cest/cest_13c.md)) and ¹⁵N
([`cest_15n`](experiments/cest/cest_15n.md)) CEST experiments. In uniformly ¹³C,
¹⁵N-labeled samples, due to extensive ¹J<sub>CC</sub> and ¹J<sub>CN</sub>
coupling network, the data analysis becomes more complicated compared with the
case when such couplings are absent. If ¹J<sub>CC</sub> or ¹J<sub>CN</sub> are
not taken into account properly it is difficult to obtain proper fitting
results, especially when the scalar couplings have a relatively large size.

ChemEx can take into account ¹J<sub>CC</sub> and ¹J<sub>CN</sub> couplings for
analyzing ¹³C or ¹⁵N CEST datasets. Since the number and size of ¹J<sub>CC</sub>
and ¹J<sub>CN</sub> for each specific ¹³C site depends on both the residue and
atom name, such information should be provided in the name of each data profile
to be taken into account properly, which is especially important for side-chain
study with ¹³C CEST.[^1]

:::info

1. `label` key should be set properly in experiment files to indicate whether
   the sample is uniformly ¹³C-labeled or not.

2. For both experiments it is implicitly assumed that the sample under study is
   always ¹⁵N-labeled, therefore it is not necessary to include `"15N"` in the
   list associated to the `label` key.

:::

[^1]:
    P. Vallurupalli, G. Bouvignies, L. E. Kay. A Computational Study of the
    Effects of ¹³C–¹³C Scalar Couplings on ¹³C CEST NMR Spectra: Towards Studies
    on a Uniformly ¹³C-Labeled Protein. _ChemBioChem_ **14**, 1709-1713 (2014).
    https://doi.org/10.1002/cbic.201300230
