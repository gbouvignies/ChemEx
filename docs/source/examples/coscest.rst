.. _example_coscest:

=======================
1HN COS-CEST experiment
=======================

This example demonstrates the analysis of datasets from 1HN COS-CEST
experiment (:ref:`coscest_1hn_ip_ap`), which is based on
in-phase/anti-phase scheme while the excitation element is replaced
with cos-shaped pulses. Similar to D-CEST scheme, COS-CEST also
performs multiple excitation therefore achieves much time-saving
compared with regular CEST. Since in COS-CEST the excitation 
region is highly restricted compared with D-CEST, it may lead
to much enhanced sensitivity due to relaxation effects in
many cases, an important application is 1HN COS-CEST 
for studying fully protonated proteins [Yuwen2018d]_.

See :file:`COSCEST_1HN_IP_AP/` under :file:`examples/Experiments/`
for this example. Note that since COS-CEST requires numerical
integration for cos-shaped pulses with multiple time steps, the
fitting speed could be much slower than other experment modules.
The :confval:`cos_res` key in experiment files should not be set 
too large, and typically should be much smaller than actually used 
experimentally so long as reasonable fitting results can be achieved.
Similar to the requirement in other in-phase/anti-phase CEST 
experiments, 1H chemical shifts provided in parameter files should 
correspond to the "actual" value in the absence of 1JHN.


.. [Yuwen2018d] T. Yuwen, G. Bouvignies, L. E. Kay, *Exploring Methods
   to Expedite the Recording of CEST Datasets Using Selective Pulse 
   Excitation*, J Magn Reson, 292:1-7, 2018. 
   https://doi.org/10.1016/j.jmr.2018.04.013
