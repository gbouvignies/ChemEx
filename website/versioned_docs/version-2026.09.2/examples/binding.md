---
sidebar_position: 9
sidebar_label: "Protein-ligand binding study"
---

# Protein-ligand binding study

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/2stBinding).

:::

This example demonstrates the study of protein-ligand binding with CPMG and CEST
experiments. In this example, CPMG and CEST experiments have been carried out on
samples with different protein/ligand molar ratios. Since the amount of ligand
is much less than protein in each sample, the protein-ligand bound state is only
sparsely populated, therefore CPMG and CEST can be used to probe the bound state
even if the resonances are not visible in the spectra.

Here, it is assumed that the same binding constant K<sub>d</sub> is shared among
different samples, which are distinguished by different `p_total` and `l_total`
keys in experiment files. The `2st_binding` kinetic model is used to fit all the
datasets with a two-step fitting scheme. In the first step only selective
residues near the protein-ligand interaction sites are fitted, and the goal is
to get a good estimate of k<sub>off</sub>; in the second step all residues are
included while k<sub>off</sub> is fixed to the estimated value in the previous
step, the major purpose is to obtain residue-specific parameters.

:::info

Further analysis of the fitting results can be carried out with the aid of
additional functions in ChemEx (refer to
[`Plotting best-fit parameters`](user_guide/additional_modules.mdx#plotting-best-fit-parameters)
subsection). During the whole fitting process K<sub>d</sub> is fixed to the
value obtained from ITC experiments, more details about this example can be
found in the reference.[^1]

:::

[^1]:
    C. Charlier, G. Bouvignies, P. Pelupessy, A. Walrant, R. Marquant, M.
    Kozlov, P. De Ioannes, N. Bolik-Coulon, S. Sagan, P. Cortes, A. K. Aggarwal,
    L. Carlier, and F. Ferrage. Structure and Dynamics of an Intrinsically
    Disordered Protein Region That Partially Folds upon Binding by
    Chemical-Exchange NMR _J. Am. Chem. Soc._, _139_, 12219-12227 (2017).
    https://doi.org/10.1021/jacs.7b05823
