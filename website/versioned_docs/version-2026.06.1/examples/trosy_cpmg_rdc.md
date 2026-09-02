---
sidebar_position: 4
sidebar_label: "Measuring Excited state RDCs with CPMG"
---

# RDC measurement for excited state with ¹⁵N TROSY CPMG

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/N15_NH_RDC).

:::

This example demonstrates the use of CPMG experiment to measure RDC parameters
for the excited state. In ¹⁵N TROSY CPMG experiment
([`cpmg_15n_tr`](../experiments/cpmg/cpmg_15n_tr)), by analyzing datasets measured
for TROSY and anti-TROSY components, the information about ¹⁵N Δϖ for both
components can be obtained, therefore RDC parameters for excited state can be
derived based on RDC parameters for ground state together with Δϖ of these two
components [^1].

Note that `antitrosy` key should be set properly in experiment files to indicate
whether the datasets are measured for TROSY or anti-TROSY component. This
example also includes datasets measured from pure in-phase CPMG experiment
([`cpmg_15n_ip`](../experiments/cpmg/cpmg_15n_ip)), which is optional and may help to
obtain even more accurate results. Note that ¹⁵N chemical shifts provided in
parameter files should correspond to the "actual" value in the absence of 1JHN,
therefore ¹⁵N chemical shifts from pure in-phase experiment should be used
instead of the TROSY version.

[^1]:
    P. Vallurupalli, D. F. Hansen, E. Stollar, E. Meirovitch, and L. E. Kay.
    Measurement of Bond Vector Orientations in Invisible Excited States of
    Proteins _Proc. Natl. Acad. Sci. USA_ **104**, 18473-18477 (2007).
    https://doi.org/10.1073/pnas.0708296104
