.. _example_rdc_measurement:

=====================================================
RDC measurement for excited state with 15N TROSY CPMG
=====================================================

This example demonstrates the use of CPMG experiment to measure
RDC parameters for the excited state. In 15N TROSY CPMG
experiment (:ref:`cpmg_15n_tr`), by analyzing datasets
measured for TROSY and anti-TROSY components, the information
about 15N Δϖ for both components can be obtained,
therefore RDC parameters for excited state can be derived
based on RDC parameters for ground state together with
Δϖ of these two components [Vallurupalli2007]_.

See :file:`N15_NH_RDC/` under :file:`examples/Combinations/`
for this example. Note that :confval:`antitrosy` key should
be set properly in experiment files to indicate whether the datasets
are measured for TROSY or anti-TROSY component. This example also
includes datasets measured from pure in-phase CPMG experiment
(:ref:`cpmg_15n_ip`), which is optional and may help to obtain
even more accurate results. Note that 15N chemical shifts provided
in parameter files should correspond to the "actual" value in the
absence of 1JHN, therefore 15N chemical shifts from pure in-phase
experiment should be used instead of the TROSY version.


.. [Vallurupalli2007] P. Vallurupalli, D. F. Hansen, E. Stollar,
   E. Meirovitch, and L. E. Kay, *Measurement of Bond Vector Orientations
   in Invisible Excited States of Proteins*, Proc Natl Acad Sci USA,
   104(47):18473-18477, 2007. https://doi.org/10.1073/pnas.0708296104
