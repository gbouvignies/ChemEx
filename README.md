ChemEx
======

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

  * [Python 3.5.x](https://www.python.org/downloads/)
  * [SciPy 0.15](https://www.scipy.org/install.html)
  * [NumPy 0.15.1](https://www.scipy.org/scipylib/download.html)
  * [Matplotlib 2.0.x](http://matplotlib.org/users/installing.html)
  * [lmfit 0.9.6](https://lmfit.github.io/lmfit-py/)

Newer versions of the above should be fine, as will be some older versions.
Please report any compatibility issues that you run into.
