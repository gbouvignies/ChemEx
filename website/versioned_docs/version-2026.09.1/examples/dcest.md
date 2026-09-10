---
sidebar_position: 7
sidebar_label: "¹⁵N D-CEST experiment"
---

# ¹⁵N D-CEST experiment

:::note Files

You can get the example files from the GitHub ChemEx page:
[regular](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N)
and
[H/D solvent exchange](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N_HD_EXCH)
studies.

:::

The two examples demonstrate the application of the ¹⁵N D-CEST experiment
([`dcest_15n`](experiments/dcest/dcest_15n.md)). D-CEST experiment accelerates
data acquisition by using the DANTE scheme to perform multifrequency irradiation
at a desired effective B<sub>1</sub> field. In D-CEST experiment, the saturation
element is replaced with DANTE pulses to perform multifrequency irradiation,
which results in dramatic savings in measurement time.[^1]

An important application of D-CEST experiment is to measure H/D solvent exchange
rates, the results should be more accurate than traditional CLEANEX experiment,
especially for studying highly flexible biomolecules such as intrinsically
disordered proteins.[^2]

:::info

Use the kinetic model `2st_hd` for H/D solvent exchange study. Besides, you
should provide the amount of D<sub>2</sub>O in the sample with the `d2o` key in
the section `[conditions]` of the **experiment.toml** file. The `d2o` key is
included in experiment files to differentiate datasets recorded with different
amounts of D<sub>2</sub>O. Note that the amount of D<sub>2</sub>O is a parameter
that can be fitted to account for potential errors during sample preparation.

:::

[^1]:
    T. Yuwen, L. E. Kay, and G. Bouvignies. Dramatic Decrease in CEST
    Measurement Times Using Multi-Site Excitation. _ChemPhysChem_ **19**,
    1707-1710 (2018). https://doi.org/10.1002/cphc.201800249

[^2]:
    T. Yuwen, A. Bah, J. P. Brady, F. Ferrage, G. Bouvignies, and L. E. Kay.
    Measuring Solvent Hydrogen Exchange Rates by Multifrequency Excitation ¹⁵N
    CEST: Application to Protein Phase Separation. _J. Phys. Chem. B_ **122**,
    11206-11217 (2018). https://doi.org/10.1021/acs.jpcb.8b06820
