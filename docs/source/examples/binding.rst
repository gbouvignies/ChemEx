.. _example_binding:

============================
Protein-ligand binding study
============================

This example demonstrates the study of protein-ligand binding with CPMG and
CEST experiments, see :file:`2stBinding/` under :file:`examples/Combinations/`
for this example. In this example several different samples with different
protein/ligand molar ratios are prepared, both CPMG and CEST experiments have
been carried out. Since the amount of ligand is much less than protein in
each sample, the protein-ligand bound state is only sparsely populated,
therefore CPMG and CEST can be used to probe the bound state even if the
resonances are not visible in the spectra. It is assumed that the same
binding constant Kd is shared among different samples, which are
distinguished by different :confval:`p_total` and :confval:`l_total` keys
in experiment files. The :confval:`2st_binding` kinetic model is used to
fit all the datasets with a two-step fitting scheme. In the first step
only selective residues near the protein-ligand interaction sites are
fitted, and the goal is to get a good estimate of koff; in the second
step all residues are included while koff is fixed to the estimated value
in the previous step, the major purpose is to obtain residue-specific
parameters. Further analysis of the fitting results can be carried out
with the aid of additional functions in ChemEx (refer to
:ref:`additional_visualize` subsection).  During the whole fitting
process Kd is fixed to the value obtained from ITC experiments, more
details about this example can be found in the reference [Charlier2017]_.


.. [Charlier2017] C. Charlier, et al., *Structure and Dynamics of an
   Intrinsically Disordered Protein Region That Partially Folds upon
   Binding by Chemical-Exchange NMR*, J Am Chem Soc, 139(35):12219-12227,
   2017. https://doi.org/10.1021/jacs.7b05823
