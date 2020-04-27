.. _example_ip_ap_cest:

=======================================
In-phase/anti-phase 1H CEST experiments
=======================================

This example demonstrates the analysis of datasets from 
in-phase/anti-phase 1H CEST experiments, such scheme has been
implemented for studying both backbone (:ref:`cest_1hn_ip_ap`) and
methyl side-chain (:ref:`cest_ch3_1h_ip_ap`) sites, which has the
advantage of fully suppressing artifacts due to NOE effects. 
The in-phase/anti-phase CEST scheme requires measuring two 
sets of CEST profiles for two separate spin states and then take
subtraction [Yuwen2017]_, [Yuwen2017a]_, an alternative scheme 
is to start from in-phase magnetization and detect anti-phase 
magnetization at the end of CEST period [Yuwen2018a]_. These two 
schemes are equivalent while the latter scheme achieves higher
sensitivity and much time-saving, therefore should be preferred.

See :file:`CEST_1HN_IP_AP/` and :file:`CEST_CH3_1H_IP_AP/` 
under :file:`examples/Experiments/` as example for backbone and
methyl side-chain study, respectively. Note that 1H chemical shifts 
provided in parameter files should correspond to the "actual" value 
in the absence of 1JHN or 1JHC, therefore a separate experiment 
that fully suppresses 1JHN or 1JHC during final detection might 
be required to obtain the correct value.


.. [Yuwen2017] T. Yuwen, A. Sekhar, L. E. Kay, *Separating Dipolar 
   and Chemical Exchange Magnetization Transfer Processes in 1H‐CEST*,
   Angew Chem Int Ed, 56(22):6122-6125, 2017. 
   https://doi.org/10.1002/anie.201610759

.. [Yuwen2017a] T. Yuwen, R. Huang, and L. E. Kay, *Probing Slow 
   Timescale Dynamics in Proteins Using Methyl 1H CEST*, J Biomol NMR,
   68:215-224, 2017. https://doi.org/10.1007/s10858-017-0121-x

.. [Yuwen2018a] T. Yuwen, and L. E. Kay, *A New Class of CEST Experiment 
   Based on Selecting Different Magnetization Components at the Start 
   and End of the CEST Relaxation Element: An Application to 1H CEST*,
   J Biomol NMR, 70:93-102, 2018. https://doi.org/10.1007/s10858-017-0161-2
