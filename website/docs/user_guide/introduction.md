---
sidebar_position: 1
---

# Introduction

## Context

Biological macromolecules such as proteins and nucleic acids are inherently
flexible molecules, whose function (or malfunction) critically depends on
interconversions between different conformational states. Understanding how
biomolecules work therefore requires a quantitative characterization of their
thermally accessible conformational landscape. The task is challenging, however,
because many of the populated conformers are only transiently formed and
marginally populated, so that they remain invisible to most standard biophysical
techniques. Over the past two decades, nuclear magnetic resonance (NMR) has
emerged as an extremely powerful tool to study these elusive states at atomic
resolution. Central to this success has been the development of chemical
exchange-based approaches, such as (CPMG, R$_{1ρ}$) relaxation dispersion and
chemical exchange saturation transfer (CEST) experiments, with the numerical
tools needed to extract the kinetic (rates) and thermodynamic (populations)
parameters associated with the exchange process and the structural information
(chemical shifts) of the sparsely populated, "invisible" excited states.


![exchange_cest_cpmg_figure](/img/exchange_cest_cpmg_figure.png)

## About ChemEx

ChemEx is an open-source [Python](https://www.python.org) application for the
analysis of NMR chemical exchange data, whose general idea is to integrate the
evolution matrix of the spin-system of interest over the pulse sequence and
extract the best-fit parameters by least-squares optimization. As ChemEx does
not rely on analytical equations to fit the data sets, any type of experiments
(e.g., D-CEST/COS-CEST) or kinetic models (e.g., 3-state exchange model) can be
simulated, and most experimental details (e.g., finite pulse width,
off-resonance effects, etc.) can be taken into account. ChemEx offers a wide
range of pulse sequences and kinetic models to choose from. Some of its main
features include multi-step fits, joint analysis of datasets from different
experiments, error analysis, grid search, etc. This documentation provides an
overview of these different features, along with a description of different
examples illustrating their use.

ChemEx is a pure [Python](https://www.python.org) package that builds upon well
established open-source packages, including [NumPy](https://numpy.org),
[SciPy](https://scipy.org), [Matplotlib](https://matplotlib.org),
[Rich](https://rich.readthedocs.io/en/stable/), and
[pydantic](https://pydantic-docs.helpmanual.io). The minimization procedure is
carried out using the [LMFIT](https://lmfit.github.io/lmfit-py/) module, which
supports many of the optimization algorithms available in
[SciPy](https://scipy.org/).
