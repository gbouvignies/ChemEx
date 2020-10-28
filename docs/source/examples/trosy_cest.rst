.. _example_trosy_cest:

=========================
15N TROSY CEST experiment
=========================

This example demonstrates the analysis of datasets from 15N TROSY
CEST experiment (:ref:`cest_15n_tr`). Compared with 15N pure in-phase
CEST experiment (:ref:`cest_15n`), the TROSY version benefits from
TROSY effects therefore the relaxation property is more favorable.
From this experiment it is possible to measure certain parameters
about excited state that are difficult to obtain from other experiments,
such as solvent exchange rate [Long2014]_.

See :file:`CEST_15N_TR/` under :file:`examples/Experiments/` for
this example. Note that :confval:`antitrosy` key should
be set properly in experiment files to indicate whether the datasets
are measured for TROSY or anti-TROSY component. The current example
includes both TROSY and anti-TROSY datasets, while it is also
possible to analyze TROSY datasets alone. Also note that 15N
chemical shifts provided in parameter files should correspond
to the "actual" value in the absence of 1JHN, therefore a
separate non-TROSY HSQC might be required to obtain the correct
value, or an alternative stratety is to correct chemical shifts
obtained from TROSY HSQC with the size of 1JHN.


.. [Long2014] D. Long, G. Bouvignies, L. E. Kay, *Measuring Hydrogen
   Exchange Rates in Invisible Protein Excited States*,
   Proc Natl Acad Sci USA, 111(24):8820-8825, 2014.
   https://doi.org/10.1073/pnas.1405011111
