---
sidebar_position: 6
sidebar_label: "¹⁵N TROSY CEST experiment"
---

# ¹⁵N TROSY CEST experiment

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N_TR).

:::

This example demonstrates the analysis of datasets from ¹⁵N TROSY CEST
experiment ([`cest_15n_tr`](../experiments/cest/cest_15n_tr)). Compared with ¹⁵N pure in-phase CEST
experiment ([`cest_15n`](../experiments/cest/cest_15n)), the TROSY version benefits from TROSY
effects therefore the relaxation property is more favorable. From this
experiment it is possible to measure certain parameters about excited state that
are difficult to obtain from other experiments, such as solvent exchange
rate.[^1]

Note that `antitrosy` key should be set properly in experiment files to indicate
whether the datasets are measured for TROSY or anti-TROSY component. The current
example includes both TROSY and anti-TROSY datasets, while it is also possible
to analyze TROSY datasets alone. Also note that ¹⁵N chemical shifts provided in
parameter files should correspond to the "actual" value in the absence of
¹J<sub>HN</sub>, therefore a separate non-TROSY HSQC might be required to obtain
the correct value, or an alternative strategy is to correct chemical shifts
obtained from TROSY HSQC with the size of ¹J<sub>HN</sub>.

[^1]:
    D. Long, G. Bouvignies, and L. E. Kay. Measuring Hydrogen Exchange Rates in
    Invisible Protein Excited States. _Proc. Natl. Acad. Sci. USA_ **111**,
    8820-8825 (2014). https://doi.org/10.1073/pnas.1405011111
