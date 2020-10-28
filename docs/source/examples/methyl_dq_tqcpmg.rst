.. _example_methyl_dq_tqcpmg:

=========================================
Methyl 1H DQ-CPMG and TQ-CPMG experiments
=========================================

This example demonstrates the application of methyl 1H DQ-CPMG
(:ref:`cpmg_ch3_1h_dq`) and TQ-CPMG (:ref:`cpmg_ch3_1h_tq`) experiments
for studying 13CH3-labeled methyl groups. Both experiments are
carried out on the same sample at single magnetic field.
The major advantage of combining these two types of experiments
is that accurate exchange parameters can already be obtained at single
magnetic field. Typically it is necessary to carry out CPMG experiments
at multiple static fields to obtain accurate chemical exchange
parameters (p\ :sub:`b` and k\ :sub:`ex`), while in DQ- and TQ-CPMG
experiments the effective field strength behaves as two and three
times of the static field, respectively. If datasets from 1H SQ-CPMG
experiment (:ref:`cpmg_ch3_1h_sq`) are available it is also
possible to perform joint fits for all of these three types
of experiments.

Both DQ-CPMG [Gopalan2018]_ and TQ-CPMG [Yuwen2016]_ experiments
collect signals from I = 3/2 manifold, the selection of different
terms of coherences is mainly achieved by different types of phase
cycling. In both experiments CPMG pulses are applied with XY-16 phase
cycle to help overcome off-resonance effects, besides, CPMG pulses
have the option to be applied as 90x180y90x composite variety, which is
indicated by :confval:`comp180_flg`.

See :file:`CPMG_CH3_1H_DQ_TQ/` under :file:`examples/Combinations/` for
this example.


.. [Gopalan2018] A. Gopalan, T. Yuwen, L. E. Kay, and P. Vallurupalli,
   *A Methyl 1H Double Quantum CPMG Experiment to Study Protein
   Conformational Exchange*, J Biomol NMR, 72:79-91, 2018.
   https://doi.org/10.1007/s10858-018-0208-z

.. [Yuwen2016] T. Yuwen, P. Vallurupalli, and L. E. Kay, *Enhancing the
   Sensitivity of CPMG Relaxation Dispersion to Conformational Exchange
   Processes by Multiple-Quantum Spectroscopy*, Angew Chem Int Ed,
   55(38):11490-11494, 2016.
   https://doi.org/10.1002/anie.201605843
