.. _example_diffusion_tqcpmg:

===================================================================
Diffusion rate measurement for excited state with methyl 1H TQ-CPMG
===================================================================

This example demonstrates the use of CPMG experiment to measure
diffusion rate constant of the excited state. In order to achieve
such purpose an additional pair of gradients should be placed around
the CPMG element, such idea has been successfully applied in the
modified version of methyl 1H TQ-CPMG experiment [Yuwen2018]_. In
order to take into account the additional pair of gradients, a
separate module (:ref:`cpmg_ch3_1h_tq_diff`) is required to fit
datasets from this experiment.

See :file:`CPMG_CH3_1H_TQ_DIFF/` under :file:`examples/Experiments/`
for such an example. In this example methyl 1H TQ-CPMG datasets
with different strengths of additional gradients are measured in
order to obtain diffusion rate constant of the excited state. Note
that those parameters about gradients around the CPMG element should
be specified in experiment files, such as the strength and duration.


.. [Yuwen2018] T. Yuwen, A. Sekhar, A. Baldwin, P. Vallurupalli,
   and L. E. Kay, *Measuring Diffusion Constants of Invisible
   Protein Conformers by Triple-Quantum 1H CPMG Relaxation Dispersion*,
   Angew Chem Int Ed, 57(51):16777-16780, 2018.
   https://doi.org/10.1002/anie.201810868
