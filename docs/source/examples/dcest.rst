.. _example_dcest:

=====================
15N D-CEST experiment
=====================

This example demonstrates the application of 15N D-CEST experiment
(:ref:`dcest_15n`), which obtains similar information about excited
state as regular CEST experiment (:ref:`cest_15n`) while achieves
much time-saving. In D-CEST experiment the saturation element is
replaced with DANTE pulses, which excite multiple sites
simultaneously therefore speed up the searching process for
excited states [Yuwen2018b]_. An important application of D-CEST
experiment is to measure H/D solvent exchange rates, the results should
be more accurate than traditional CLEANEX experiment, especially for
studying certain biomolecules such as disordered proteins [Yuwen2018c]_.

See :file:`DCEST_15N/` and :file:`DCEST_15N_HD_EXCH/`
under :file:`examples/Experiments/` as example for regular and H/D
solvent exchange study with 15N D-CEST, respectively.  Note that
the same experiment module is used in both examples, while
:confval:`2st_hd` kinetic model should be used for H/D solvent
exchange study, besides, certain additional information should be
provided such as :confval:`d2o` key which indicates the amount
of D\ :sub:`2`\ O in the sample. Note that :confval:`d2o` key
is included in experiment files such that samples with different
amount of D\ :sub:`2`\ O can be analyzed at the same time, while
the amount of D\ :sub:`2`\ O is allowed to be fitted to account
for potential errors during sample preparation.


.. [Yuwen2018b] T. Yuwen, L. E. Kay, G. Bouvignies, *Dramatic Decrease in
   CEST Measurement Times Using Multi-Site Excitation*, ChemPhysChem,
   19(14):1707-1710, 2018. https://doi.org/10.1002/cphc.201800249

.. [Yuwen2018c] T. Yuwen, A. Bah, J. P. Brady, F. Ferrage, G. Bouvignies,
   and L. E. Kay, *Measuring Solvent Hydrogen Exchange Rates by
   Multifrequency Excitation 15N CEST: Application to Protein Phase
   Separation*, Journal of Physical Chemistry B, 122(49):11206-11217,
   2018. https://doi.org/10.1021/acs.jpcb.8b06820
