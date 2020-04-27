.. _chemex_develop:

==================
Developement Guide
==================

This guide provides instruction for setting up the environment for 
developing ChemEx, an overview of the project layout, the summary 
of all modules and the contribution process.


Requirements
------------

To create an environment for developing ChemEx the following additional
programs should be installed aside form the required Python modules:

* `git <https://git-scm.com>`_

* `Sphinx <https://www.sphinx-doc.org/>`_


Source code
-----------

ChemEx uses `github <https://github.com>`_ for source code hosting.  
For access to the source code, see the 
`ChemEx github site <https://github.com/gbouvignies/chemex>`_. 
To check out the latest version of ChemEx use :command:`git`:

.. code-block:: console
    
   $ git clone git://github.com/gbouvignies/chemex.git

ChemEx is a pure Python module, the root directory can be included in 
``PYTHONPATH`` directly, or a symbolic link can be added into the 
``site-packages/`` directory of the installed Python program. In this
way any modification to the ChemEx source tree will be picked up when
running ChemEx.


Project layout
--------------

The directory layout of the ChemEx project is as follows:

* ``chemex/``: Contain source code for the project.

* ``docs/``: Contain the setup file and source code for building the
  ChemEx documentation using `Sphinx`_.

* ``examples/``: Contain numerous examples to demonstrate the application
  of each experiment module in ChemEx.

* ``sandbox/`` (optional): Suggested location to store code, data and 
  other stuff which are not ready to be included in ChemEx. This directory
  is not required and will be ignored by :command:`git` using the default
  :file:`.gitignore` file, it is mainly created for developement purposes.


Chemex modules
--------------

.. toctree::
   :maxdepth: 1

   modules/chemex/chemex
   modules/chemex/cli
   modules/chemex/fitting
   modules/chemex/helper


Containers modules
------------------

.. toctree::
   :maxdepth: 1

   modules/containers/cest
   modules/containers/conditions
   modules/containers/cpmg
   modules/containers/experiment
   modules/containers/helper
   modules/containers/noise
   modules/containers/plot
   modules/containers/relaxation
   modules/containers/shift


Experiments modules
-------------------

.. toctree::
   :maxdepth: 1

   modules/experiments/init
   modules/experiments/helper

Refer to :ref:`chemex_experiments` section for more details about each 
experiment module. The ``config/`` subdirectory contains sample config 
files for each experiment.


Nmr modules
-----------

.. toctree::
   :maxdepth: 1

   modules/nmr/constants
   modules/nmr/liouvillian
   modules/nmr/propagator
   modules/nmr/rates
   modules/nmr/spin_system


Parameters modules
------------------

.. toctree::
   :maxdepth: 1

   modules/parameters/helper
   modules/parameters/kinetics
   modules/parameters/liouvillian
   modules/parameters/name
   modules/parameters/settings


Tools modules
-------------

.. toctree::
   :maxdepth: 1

   modules/tools/pick_cest
   modules/tools/plot_param


Suggestions
-----------

When working with the ChemEx source code please consider the following when
making any updates.

* Coding style: When adding a new module for a new experiment please try
  to follow the coding style of other existing modules. The easiest way is 
  to create a copy from an existing experiment module that is similar to 
  the new experiment, and then make necessary modifications.

* Documentation: All public functions and classes should have docstring that
  follows the `numpydoc docstring standard 
  <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_.
  Private functions and classes may have shorter dostrings. The ChemEx 
  documentation is built using `Sphinx`_, which translates 
  `reStructuredText <https://docutils.sourceforge.io/rst.html>`_ formatted 
  documents (including docstring) into html/pdf/epub. When adding new 
  functions, classes or parameters to ChemEx please update the docstring and 
  make any necessary changes to the Sphinx files in the ``docs/``
  directory.
  
* Examples: Numerous example showing the real world use of ChemEx are 
  provided in the ``examples/`` directory. Contributions of additional
  examples are welcome and appreciated.  


Reporting bugs
--------------

The preferred location for submitting feature requests and reporting bugs
is the `github issue tracker <https://github.com/gbouvignies/chemex/issues>`_.
Reports are also welcomed on the 
`ChemEx mailing list <https://groups.io/g/chemex>`_ or by
contacting `Guillaume Bouvignies <gbouvignies@gmail.com>`_ directly.


Contributions
-------------

Contribution of source code or examples to ChemEx is welcomed provided that
the contents can be distributed under the 
`New BSD License <https://opensource.org/licenses/BSD-3-Clause>`_. The 
preferred method for contributing is to create a feature branch on a 
github fork of ChemEx and submit a pull request, although patches 
are also accepted. Refer to the 
`NumPy git workflow <https://docs.scipy.org/doc/numpy/dev/gitwash/index.html>`_
for more details on how to submit a pull request or prepare a patch.
