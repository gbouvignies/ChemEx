.. _chemex_start:

===============
Getting Started
===============

What is ChemEx
--------------

ChemEx a Python module for analyzing datasets from different types 
of NMR experiments for chemical exchange process study, the most 
commonly used experiments include Carr–Purcell–Meiboom–Gill
relaxation dispersion (CPMG) and Chemical Exchange Saturation
Transfer (CEST).


Where to get ChemEx
-------------------

The source code and documentation of ChemEx for all platforms 
(Linux, OS X and Windows) are available for 
`download <https://github.com/gbouvignies/ChemEx/releases>`_ from
`github <https://github.com>`_.


Requirements
------------

ChemEx requires `NumPy <https://numpy.org>`_, 
`SciPy <https://www.scipy.org>`_, `matplotlib <https://matplotlib.org/>`_
and `LMfit <https://lmfit.github.io/lmfit-py/>`_ to be installed.  
An easy way of obtaining these Python modules is to use a
Python distribution that provides these packages, such as 
`Anaconda <https://www.anaconda.com/distribution/>`_ or 
`Intel® Distribution for Python <https://software.intel.com/en-us/distribution-for-python>`_.

.. attention::
   Due to the highly extensive numerical calculation feature during 
   the fitting process, the performance of ChemEx highly depends on 
   the efficiency of numerical calculation modules.  It is highly 
   recommended to install `SciPy`_/`NumPy`_ modules compiled with 
   Intel® Math Kernel Library (Intel® MKL), which can be obtained 
   from `Anaconda`_ or `Intel® Distribution for Python`_. 
   The `SciPy`_/`NumPy`_ modules installed with :command:`pip` are 
   usually not optimized and may show much lower performance in
   many cases.


Installation
------------

.. highlight:: console

The easiest way to install ChemEx is via 
`conda <https://conda.io/en/latest/>`_ (which comes with 
`Anaconda`_) from the command line (indicated with ``$``)::

   $ conda install -c conda-forge chemex

Note that there is minimum version of Python required. If the version
of Python is less than 3.7, ChemEx can be installed in a separate conda
environment which enforces the use of Python 3.7::

   $ conda create -c conda-forge -n chemex python=3.7 chemex
   $ conda activate chemex

ChemEx is also available via the Python package index using :command:`pip`::

   $ pip install chemex

The development version can be installed directly from 
`github`_ via :command:`pip`::

   $ pip install git+https://github.com/gbouvignies/chemex.git

.. highlight:: default

.. _ChemEx: https://github.com/gbouvignies/ChemEx/
