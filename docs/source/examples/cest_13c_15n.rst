.. _example_cest_13c_15n:

=============================================
CEST experiments for 13C, 15N-labeled samples
=============================================

This example demonstrates the analysis of CEST datasets measured for 
uniformly 13C, 15N-labeled samples, which include both 13C 
(:ref:`cest_13c`) and 15N (:ref:`cest_15n`) CEST experiments. In 
uniformly 13C, 15N-labeled samples, due to extensive 1JCC and 1JCN
coupling network, the data analysis becomes more complicated compared
with the case when such couplings are absent. If 1JCC or 1JCN are
not taken into account properly it is difficult to obtain proper
fitting results, especially when the scalar couplings have
relatively large size.

ChemEx allows taking into account 1JCC and 1JCN couplings for analyzing
13C or 15N CEST datasets. Since the number and size of 1JCC and 1JCN 
for each specific 13C site depends on both the residue type and atom 
name, such information should be provided in the name of each data 
profile in order to be taken into account properly, which is especially
important for side-chain study with 13C CEST [Vallurupalli2014]_.

See :file:`CEST_13C_LABEL_CN/` and :file:`CEST_15N_LABEL_CN/` 
under :file:`examples/Experiments/` as example for 13C and 15N CEST
for studying uniformly 13C, 15N-labeled samples, respectively. 
Note that :confval:`label` key should be set properly 
in experiment files to indicate whether the sample is uniformly
13C-labeled or not. Also note that for both experiments it is 
implicitly assumed that the sample under study is always 
15N-labeled, therefore it is not necessary to include ``"15N"`` 
in :confval:`label` key.


.. [Vallurupalli2014] P. Vallurupalli, G. Bouvignies, L. E. Kay, 
   *A Computational Study of the Effects of 13C–13C Scalar Couplings 
   on 13C CEST NMR Spectra: Towards Studies on a Uniformly 13C-Labeled 
   Protein*, ChemBioChem, 14(14):1709-1713, 2014. 
   https://doi.org/10.1002/cbic.201300230
