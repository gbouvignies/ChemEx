ChemEx
======

[![Build Status](https://travis-ci.org/gbouvignies/ChemEx.svg?branch=develop)](https://travis-ci.org/gbouvignies/ChemEx)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/da6be0c1863c4655b1ee15006bc90f36)](https://app.codacy.com/app/gbouvignies/chemex?utm_source=github.com&utm_medium=referral&utm_content=gbouvignies/chemex&utm_campaign=Badge_Grade_Dashboard)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)


Overview
---------

ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

Installation
------------

The easiest way to install `chemex` is via [conda](http://conda.pydata.org):
```bash
conda install -c conda-forge chemex
```
If your version of python is less than 3.5, you can also install `chemex` in a separate conda environment enforcing the use of python 3.7:
```bash
conda create -c conda-forge -n chemex python=3.7 chemex
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

Dependencies
------------

  * [Python>=3.5](https://www.python.org/downloads/)
  * [SciPy>=1.0](https://www.scipy.org/install.html)
  * [NumPy>=1.0](https://www.scipy.org/scipylib/download.html)
  * [Matplotlib>=2.0](http://matplotlib.org/users/installing.html)
  * [LmFit>=0.9.11](https://lmfit.github.io/lmfit-py/)
  * [ASTEVAL>=0.9.11](https://github.com/newville/asteval)
