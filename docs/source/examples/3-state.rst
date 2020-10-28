.. _example_3_state:

==============================
3-state exchange kinetic model
==============================

This example demonstrates the study of systems with 3-state
exchange using ChemEx, see :file:`DCEST_15N_3States/` under
:file:`examples/Experiments/` for this example. For systems with
3-state exchange, two separate minor dips may show up in CEST
profiles, while in general it is difficult to identify the correct
exchange model without additional information. In this example
DRD-CEST experiments have been carried out, with the aid of
additional DRD-CEST datasets it is possible to distinguish
between the two possible exchange models. The :confval:`3st`
kinetic model should be used to fit the datasets for 3-state
exchange. See :file:`README.txt` in this example for explanation
about different fitting schemes, more details about this example
can be found in the reference [Vallurupalli2019]_.


.. [Vallurupalli2019] P. Vallurupalli, V. K. Tiwari, S. Ghosh,
   *A Double-Resonance CEST Experiment To Study Multistate Protein
   Conformational Exchange: An Application to Protein Folding*,
   J Phys Chem Lett, 10(11):3051-3056, 2019.
   https://doi.org/10.1021/acs.jpclett.9b00985
