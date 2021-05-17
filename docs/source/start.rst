.. _chemex_start:

===============
Getting Started
===============

What is ChemEx
--------------

ChemEx a Python module for analyzing datasets from different types
of NMR experiments to study chemical exchange processes. The most
commonly used experiments include Carr–Purcell–Meiboom–Gill (CPMG)
relaxation dispersion and Chemical Exchange Saturation
Transfer (CEST).


Where to get ChemEx
-------------------

The source code and documentation of ChemEx for all platforms
(Linux, OS X and Windows) is available for
`download <https://github.com/gbouvignies/ChemEx/releases>`_ from
`GitHub <https://github.com>`_.


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
   Due to the extensive numerical calculations during the fitting
   process, the performance of ChemEx highly depends on the efficiency
   of numerical calculation modules. It is highly recommended to
   install `SciPy`_/`NumPy`_ modules compiled with Intel® Math Kernel
   Library (Intel® MKL), which can be obtained from `Anaconda`_ or
   `Intel® Distribution for Python`_. The `SciPy`_/`NumPy`_ modules
   installed with :command:`pip` are usually not optimized and may show
   much lower performance in many cases.


Installation
------------

.. highlight:: console

The easiest way to install ChemEx is via
`conda <https://conda.io/en/latest/>`_ (which comes with
`Anaconda`_) from the command line (indicated with ``$``)::

   $ conda install -c conda-forge chemex

Note that the minimum required version of Python is 3.8. If the version
of Python is less than 3.8, ChemEx can be installed in a separate conda
environment, which enforces the use of Python 3.8::

   $ conda create -c conda-forge -n chemex python=3.8 chemex
   $ conda activate chemex

ChemEx is also available via the Python package index using :command:`pip`::

   $ pip install chemex

The development version can be installed directly from
`GitHub`_ via :command:`pip`::

   $ pip install git+https://github.com/gbouvignies/chemex.git

.. highlight:: default

.. _ChemEx: https://github.com/gbouvignies/ChemEx/
