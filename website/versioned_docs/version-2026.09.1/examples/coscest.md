---
sidebar_position: 8
sidebar_label: "¹Hᴺ COS-CEST experiment"
---

# ¹Hᴺ COS-CEST experiment

:::note Files

You can get
[the example files from the GitHub ChemEx page](https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/COSCEST_1HN_IP_AP).

:::

This example demonstrates the analysis of datasets from ¹Hᴺ COS-CEST experiment
([`coscest_1hn_ip_ap`](experiments/dcest/coscest_1hn_ip_ap.md), which is based
on in-phase/anti-phase scheme while the excitation element is replaced with
cos-modulated excitation fields. Similar to D-CEST scheme, COS-CEST performs
multiple excitations, leading to significant time-savings compared to regular
CEST. Cos-modulated excitation fields offer some advantages over the DANTE
sequence in systems where the probe spins in question are dipolar or scalar
coupled to like spins that can therefore be simultaneously perturbed by the
broad-banded nature of the DANTE element. An important application is ¹Hᴺ
COS-CEST for the study of slow chemical exchange processes in fully protonated
proteins.[^1]

:::info

Since COS-CEST requires numerical integration for cos-shaped pulses with
multiple time steps, the fitting speed could be much slower than other
experiment modules. The `cos_res` key in experiment files should not be set too
large, and typically should be much smaller than actually used experimentally so
long as reasonable fitting results can be achieved. Similar to the requirement
in other in-phase/anti-phase CEST experiments, ¹H chemical shifts provided in
parameter files should correspond to the "actual" value in the absence of
¹J<sub>HN</sub>.

:::

[^1]:
    T. Yuwen, G. Bouvignies, and L. E. Kay. Exploring Methods to Expedite the
    Recording of CEST Datasets Using Selective Pulse Excitation _J. Magn.
    Reson._ **292**, 1-7 (2018). https://doi.org/10.1016/j.jmr.2018.04.013
