.. _example_methyl_sqcpmg:

============================
Methyl 1H SQ-CPMG experiment
============================

This example demonstrates the application of methyl 1H SQ-CPMG 
experiment (:ref:`cpmg_ch3_1h_sq`) for studying 13CH3-labeled methyl 
groups, which measures 1H Δϖ while still preserving the methyl-TROSY 
effects [Yuwen2019]_, the sensitivity is much higher than 
old experiments for similar purposes [Tugarinov2007]_. This experiment 
can be combined with MQ-CPMG experiment (:ref:`cpmg_ch3_mq`) such that 
Δϖ for both 13C and 1H can be obtained simultaneously.

In this experiment signals from both I = 3/2 and I = 1/2 manifolds
are collected, in order to get rid of potential artifacts due to 
off-resonance effects CPMG pulses are applied with XY-4 phase cycle, 
besides, each CPMG pulse can be applied as 90x240y90x composite variety
to achieve even better off-resonance performance, which is indicated by 
:confval:`comp180_flg`. In order to obtain flat relaxation dispersion 
profiles in the absence of chemical exchange process, relaxation 
compensation pulses should be applied and :confval:`ncyc_max` 
should be set properly.

See :file:`CPMG_CH3_1H_SQ/` under :file:`examples/Experiments/` for
this example. In this example a selective set of residues are used
for obtaining global parameters (p\ :sub:`b` and k\ :sub:`ex`) first, 
and then residue-specific parameters for each residue are obtained 
with global parameters fixed.


.. [Yuwen2019] T. Yuwen, R. Huang, P. Vallurupalli, and L. E. Kay, 
   *A Methyl-TROSY-Based 1H Relaxation Dispersion Experiment for
   Studies of Conformational Exchange in High Molecular Weight
   Proteins*, Angew Chem Int Ed, 58(19):6250-6254, 2019. 
   https://doi.org/10.1002/anie.201900241

.. [Tugarinov2007] V. Tugarinov2007, and L. E. Kay, *Separating 
   Degenerate 1H Transitions in Methyl Group Probes for 
   Single-Quantum 1H-CPMG Relaxation Dispersion NMR Spectroscopy*,
   J Am Chem Soc, 129(30):9514-9521, 2007. 
   https://doi.org/10.1021/ja0726456
