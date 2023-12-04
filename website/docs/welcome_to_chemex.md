---
sidebar_position: 1
---

# Welcome to ChemEx

## Quick Overview

ChemEx is a comprehensive tool for analyzing NMR experimental data to characterize chemical exchange processes. This guide will help you get started with ChemEx, covering everything from installation to initial usage. Primarily, ChemEx facilitates analysis of experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

Before installing ChemEx, ensure you have Python installed. For a hassle-free setup, especially for beginners, we recommend the [Anaconda Distribution](https://www.anaconda.com/distribution/). It comes bundled with Python, Numpy, and a suite of scientific computing tools.

## Installation Methods{#installation}

ChemEx offers flexible installation options: through **conda**, **pip**, or directly [from the source](https://github.com/gbouvignies/ChemEx/). Choose the method that best suits your setup.

### Using conda

For an optimal experience, install ChemEx using **conda**. This method is ideal as it handles dependencies efficiently and is available in all versions of [Anaconda](https://www.anaconda.com/distribution/) and [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

```shell
# Create and activate a new environment for ChemEx
conda create -n chemex
conda activate chemex

# Ensure Python 3.11 is installed
conda install python=3.11

# Add conda-forge channel and install ChemEx
conda config --env --add channels conda-forge
conda install chemex
```

### Using pip

Alternatively, you can install ChemEx using pip, available via the [Python Package Index](https://pypi.org/project/chemex/).

```shell
pip install chemex
```

:::tip

Utilize a virtual environment when installing with pip. [Here's a helpful guide](https://dev.to/bowmanjd/python-tools-for-managing-virtual-environments-3bko#howto).

:::

### From source

For those interested in the latest features and updates, install the development version of ChemEx directly from [GitHub](https://github.com/gbouvignies/ChemEx/).

```shell
pip install git+https://github.com/gbouvignies/ChemEx.git
```


## Performance Optimization

ChemEx performance is significantly enhanced by using numerically optimized modules. We recommend installing [NumPy](https://numpy.org) and [SciPy](https://scipy.org) with Intel® Math Kernel Library (Intel® MKL) for the best performance. These optimized modules are available via [Anaconda](https://www.anaconda.com/distribution/) or the [Intel® Distribution for Python](https://software.intel.com/en-us/distribution-for-python). Note that the standard pip installations of these modules might not offer the same level of optimization.
