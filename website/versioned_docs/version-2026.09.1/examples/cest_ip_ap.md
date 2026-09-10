---
sidebar_position: 7
sidebar_label: "¹H CEST experiments"
---

# In-phase/anti-phase ¹H CEST experiments

:::note Files

You can get the example files from the GitHub ChemEx page:
[backbone](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_1HN_IP_AP)
and
[methyl side-chain](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_CH3_1H_IP_AP)
sites.

:::

This example demonstrates the analysis of datasets from in-phase/anti-phase ¹H
CEST experiments, such scheme has been implemented for studying both backbone
([`cest_1hn_ip_ap`](../experiments/cest/cest_1hn_ip_ap)) and methyl side-chain
([`cest_ch3_1h_ip_ap`](../experiments/cest/cest_ch3_1h_ip_ap)) sites, which has the
advantage of fully suppressing artifacts due to NOE effects. The
in-phase/anti-phase CEST scheme requires measuring two sets of CEST profiles for
two separate spin states and then take subtraction,[^1]<sup>,</sup>[^2] an
alternative scheme is to start from in-phase magnetization and detect anti-phase
magnetization at the end of CEST period.[^3] These two schemes are equivalent
while the latter scheme achieves higher sensitivity and much time-saving,
therefore should be preferred.

:::info

Note that ¹H chemical shifts provided in parameter files should correspond to
the "actual" value in the absence of ¹J<sub>HN</sub> or ¹J<sub>HC</sub>,
therefore a separate experiment that fully suppresses ¹J<sub>HN</sub> or
¹J<sub>HC</sub> during final detection might be required to obtain the correct
value.

:::

[^1]:
    T. Yuwen, A. Sekhar, and L. E. Kay. Separating Dipolar and Chemical Exchange
    Magnetization Transfer Processes in 1H‐CEST. _Angew. Chem. Int. Ed._ **56**,
    6122-6125 (2017). https://doi.org/10.1002/anie.201610759

[^2]:
    T. Yuwen, R. Huang, and L. E. Kay, Probing Slow Timescale Dynamics
    in Proteins Using Methyl ¹H CEST. _J. Biomol. NMR_ **68**, 215-224 (2017).
    https://doi.org/10.1007/s10858-017-0121-x

[^3]:
    T. Yuwen, and L. E. Kay, A New Class of CEST Experiment Based on Selecting
    Different Magnetization Components at the Start and End of the CEST
    Relaxation Element: An Application to ¹H CEST. _J. Biomol. NMR_ **70**,
    93-102 (2018). https://doi.org/10.1007/s10858-017-0161-2
