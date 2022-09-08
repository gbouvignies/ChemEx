# ChemEx

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

## Overview

ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

## Installation

The easiest way to install `chemex` is via [conda](http://conda.pydata.org):

```bash
conda install -c conda-forge chemex
```

If your version of python is less than 3.9, you can also install `chemex` in a separate conda environment enforcing the use of python 3.9+:

```bash
conda create -c conda-forge -n chemex python=3.10 chemex
conda activate chemex
```

`chemex` is also available via the [Python package index](https://pypi.python.org/pypi/chemex) using `pip`:

```bash
pip install chemex
```

The development version can be installed directly from github via `pip`:

```bash
pip install git+https://github.com/gbouvignies/chemex.git
```
