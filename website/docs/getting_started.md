---
sidebar_position: 1
lastmod: 2022-04-19T20:32:18.485Z
---

# Getting started

## Overview

ChemEx is a tool for analyzing datasets from different types of NMR experiments
to study chemical exchange processes. The most commonly used experiments include
Carr–Purcell–Meiboom–Gill (CPMG) relaxation dispersion and Chemical Exchange
Saturation Transfer (CEST).

## Installation

The only prerequisite for installing ChemEx is Python itself. If you don’t have
Python yet and want the simplest way to get started, we recommend you use the
[Anaconda](https://www.anaconda.com/distribution/) Distribution -- it includes
Python, Numpy, and many other commonly used packages for scientific computing
and data science.

ChemEx can be installed with **conda**, with **pip**, or
[from source](https://github.com/gbouvignies/ChemEx/).

import Tabs from '@theme/Tabs'; import TabItem from '@theme/TabItem';

<Tabs>
<TabItem value="conda" label="Using conda" default>

The recommended way to install ChemEx is via **conda** (which is included in all
versions of [Anaconda](https://www.anaconda.com/distribution/) and
[Miniconda](https://docs.conda.io/en/latest/miniconda.html)):

```shell
# Best practice, use an environment rather than install in the base env
conda create -n chemex
conda activate chemex

# If the default version of Python is less than 3.9
conda install python=3.9

# ChemEx is available from the conda-forge channel
conda config --env --add channels conda-forge

# The actual install command
conda install chemex
```

</TabItem>
<TabItem value="pip" label="Using pip">

ChemEx is also available via the
[Python Package Index](https://pypi.org/project/chemex/) using pip:

```shell
pip install chemex
```

Also when using pip, it's good practice to use a virtual environment -- see
[this guide](https://dev.to/bowmanjd/python-tools-for-managing-virtual-environments-3bko#howto)
for details on using virtual environments.

:::tip

Due to the extensive numerical calculations during the fitting process, the
performance of ChemEx highly depends on the efficiency of numerical calculation
modules. It is highly recommended to install
[NumPy](https://numpy.org)/[SciPy](https://scipy.org) modules compiled with
Intel® Math Kernel Library (Intel® MKL), which can be obtained from
[Anaconda](https://www.anaconda.com/distribution/) or
[Intel® Distribution for Python](https://software.intel.com/en-us/distribution-for-python).
The [NumPy](https://numpy.org)/[SciPy](https://scipy.org) modules installed with
**pip** are usually not optimized and may show lower performance.

:::

</TabItem>
<TabItem value="source" label="From source">

The development version can be installed directly from
[Github](https://github.com/gbouvignies/ChemEx/) via **pip**:

```shell
pip install git+https://github.com/gbouvignies/ChemEx.git
```

</TabItem>
</Tabs>
