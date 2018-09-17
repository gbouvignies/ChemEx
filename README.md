ChemEx
======

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/da6be0c1863c4655b1ee15006bc90f36)](https://app.codacy.com/app/gbouvignies/chemex?utm_source=github.com&utm_medium=referral&utm_content=gbouvignies/chemex&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/gbouvignies/chemex.svg?branch=develop)](https://travis-ci.org/gbouvignies/chemex)


Overview
---------

ChemEx is an analysis program for chemical exchange detected by NMR.

It is designed to take almost any kind of NMR data to aid the analysis,
but the principle techniques are CPMG relaxation dispersion and Chemical
Exchange Saturation Transfer.

Quick install
-------------

In a clean directory run:

    git clone https://github.com/gbouvignies/chemex

Navigate to the chemex directory and run:

    python setup.py install --user

Quick update
------------

To get the latest code using git:

    git pull

Then simply re-run the install command above.


Prerequisites
-------------

You should have installed on your system:

  * [Python>=3.5](https://www.python.org/downloads/)
  * [SciPy>=1.0](https://www.scipy.org/install.html)
  * [NumPy>=1.0](https://www.scipy.org/scipylib/download.html)
  * [Matplotlib>=2.0](http://matplotlib.org/users/installing.html)
  * [LmFit>=0.9.11](https://lmfit.github.io/lmfit-py/)
  * [ASTEVAL>=0.9.11](https://github.com/newville/asteval)

