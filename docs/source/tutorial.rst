.. _chemex_tutorial:

========
Tutorial
========

Introduction
------------

This tutorial provides an overview of some major features of ChemEx.
ChemEx achieves the fitting purpose by minimizing a predefined *χ*\ :sup:`2`
target function, which typically has the following form:

.. math::
   \chi^2 = \sum_{i}{\left(\frac{I^{expt}_i-I^{calc}_i}{\sigma^{expt}_i}\right)^2}

For any given sets of parameters, ChemEx performs numerical simulation and
calculates the corresponding *χ*\ :sup:`2`. The aim of the fitting
process is to locate the sets of parameters that generate *χ*\ :sup:`2` minimum,
which is carried out with Levenberg–Marquardt non-linear optimization
using `LMfit <https://lmfit.github.io/lmfit-py/>`_ module. Unlike most other
programs for chemical exchange data analysis, ChemEx does not fit datasets
with analytical equations therefore most experimental details (e.g. finite
pulse width, off-resonance effects, etc.) can be taken into account.
Besides, for certain types of experiments (e.g. D-CEST/COS-CEST) or
kinetic models (e.g. 3-state exchange model) it is difficult to fit
with analytical equations, therefore numerical fitting
is the only way to obtain exchange parameters in such cases.


Running ChemEx
--------------

ChemEx is designed to use with the command line (indicated with ``$``), a
typical command for running ChemEx for data fitting purpose is like this:

.. code-block:: console

   $ chemex fit -e <FILE> \
                -p <FILE> \
                -m <FILE> \
                -d <MODEL> \
                -o <DIR>

Such command is usually saved in a shell script (typically called
:file:`run.sh` in most examples) to save some typing efforts. The
first argument determines the major purpose of running ChemEx, and
should be one of the following:

================   ========================================================
``info``           Show experiments that can be fit
``config``         Show sample configuration files of experiment modules
``fit``            Start a fit
``simulate``       Start a simulation
``pick_cest``      Plot CEST profiles for interactive dip picking
``plot_param``     Plot selected parameters based on fitting output results
================   ========================================================

The second argument and afterwards can be one of the following:

======================  ===============================================================
``-h`` or ``--help``    Show help message
``-e`` ``<FILE>``       Input files containing experimental setup and data location
``-p`` ``<FILE>``       Input files containing the initial values of fitting parameters
``-m`` ``<FILE>``       Input file containing the fitting method
``-d`` ``<MODEL>``      Exchange model used to fit the datasets
``-o`` ``<DIR>``        Directory for output files
``--plot``              Plotting level (``nothing``, ``normal`` (default), ``all``)
``--include`` ``<ID>``  Residue(s) to include in the fit
``--exclude`` ``<ID>``  Residue(s) to exclude from the fit
``--mc`` ``<N>``        Number of Monte-Carlo simulations
``--bs`` ``<N>``        Number of Bootstrap simulations
======================  ===============================================================

Usually multiple arguments should be specified in order to provide
complete information for running ChemEx.

.. tip::
   For those arguments requiring input files (e.g. ``-p``), it is
   possible to provide multiple files using command line wildcard
   characters (e.g. ``*``), instead of specifying the name of each
   file individually.


File formats
------------

The input and output files of ChemEx use
`TOML <https://en.wikipedia.org/wiki/TOML>`_ file format, the syntax
of which mainly consists of :confval:`key = value` pairs,
:confval:`[section names]` and :confval:`# comments`. Detailed
description about this file format can be found in the
`TOML github site <https://github.com/toml-lang/toml>`_.


Experiment files
----------------

The experiment files (indicated with ``-e``) contain information about
the experiments that have been carried out, such as experiment name and
Larmor frequency. An example experiment file looks like this:

.. literalinclude:: _static/experiment.toml
   :language: toml

Several most commonly used keys among different types of experiments are
summarized as below:

+-------------------------+--------------------------+-------------------------------------------+
| Section Name            | Key                      |     Description                           |
+=========================+==========================+===========================================+
| :confval:`[experiment]` |  :confval:`name`         |     Experiment name                       |
+                         +--------------------------+-------------------------------------------+
|                         |  :confval:`carrier`      |  RF carrier of studied nuclei in ppm      |
+-------------------------+--------------------------+-------------------------------------------+
| :confval:`[conditions]` | :confval:`h_larmor_frq`  |   Magnetic field strength in MHz          |
+                         +--------------------------+-------------------------------------------+
|                         | :confval:`label`         |   Labeling scheme of the sample           |
+-------------------------+--------------------------+-------------------------------------------+
| :confval:`[data]`       | :confval:`path`          |  The directory containing the data files  |
+                         +--------------------------+-------------------------------------------+
|                         | :confval:`error`         |  The method for error estimation          |
+                         +--------------------------+-------------------------------------------+
|                         | :confval:`profiles`      |  The name of each data file               |
+-------------------------+--------------------------+-------------------------------------------+

Each experiment may also contain certain specific keys, the complete list
of keys can be found in the sample config file generated with:

.. code-block:: console

   $ chemex config <NAME>

Detailed descriptions about the meaning of each key are included as
comments in sample config files.

.. important::
   Under :confval:`profiles` section, the name of each dataset should be
   properly chosen according to the spin system for each experiment, e.g.
   :confval:`G2N` is suitable for experiments based on single-spin
   system, while :confval:`G2N-HN` is suitable for experiments based
   on two-spin system, etc. Note that there are some exceptions in
   certain cases though.

.. note::
   It is possible to combine results from multiple types of experiments
   to fit together, see :file:`2stBinding/` under
   :file:`examples/Experiments/` for such
   :ref:`an example <example_binding>`.


Data files
----------

The location of data files is specified within experiment files. Data
files typically have three columns containing the following information:

+---------------------------+--------------+---------------+-----------------+
| Experiment type           | First column | Second column |  Third column   |
+===========================+==============+===============+=================+
|    CPMG                   |   ncyc_cp    |   Intensity   |   Uncertainty   |
+---------------------------+--------------+               +                 +
|    CEST/D-CEST/COS-CEST   | Offset (Hz)  |               |                 |
+---------------------------+--------------+               +                 +
|    Relaxation             |   Time (s)   |               |                 |
+---------------------------+--------------+---------------+-----------------+

An example data file looks like this:

.. literalinclude:: _static/data.out
   :language: python

Data files are just peak intensity tables that can be created from
the output of many peak fitting programs, such as the :command:`autoFit`
subroutine in `NMRPipe <https://www.ibbr.umd.edu/nmrpipe/index.html>`_
program.

.. note::
   The input data files for :ref:`shift_experiments` have special format,
   refer to the :file:`Shifts/` example under :file:`examples/Combinations/`
   to learn how to create data files for :ref:`shift_experiments`.


Parameter files
---------------

The parameter files (indicated with ``-p``) contain initial estimates of
parameters to be used during the fitting process, which typically look
like this:

.. literalinclude:: _static/parameters.toml
   :language: toml

The parameter values under :confval:`[GLOBAL]` section apply to all residues,
while residue-specific parameters can be specified under other sections with
the parameter name, such as :confval:`[CS_A]`. Multiple parameter files can
be provided if necessary.

.. note::
   If certain required parameter is not specified in parameter files, a
   default value will be assigned as the initial value,  which is
   determined by each specific experiment module.

.. attention::
   Due to the multidimensional searching feature of *χ*\ :sup:`2`
   minimization process, it is essential to set suitable initial value
   for each parameter to avoid getting trapped in a local minimum.

.. tip::
   Setting model-free parameters (e.g. :confval:`TAUC_A`) is a simple way
   to obtain initial estimates of relaxation parameters
   (e.g. :confval:`R1_A`, :confval:`R2_A`, etc.). For every 2.6 kDa molecular
   weight, the overall tumbling time is approximately 1 ns at T = 300 K for
   biomolecules in H\ :sub:`2`\ O. Assuming similar molecular structure at
   different conditions, the overall tumbling time is proportional to η/T,
   where η is solution viscosity and T is temperature in Kelvin.

.. tip::
   In parameter files, each parameter not only can have the initial value
   specified, but also the upper and lower bounds during the fitting process.
   Such boundary helps to prevent over-fitting that certain parameters
   are enforced to unreasonable values in order to minimize *χ*\ :sup:`2`.
   With the use of boundary, each parameter key should look like
   [:confval:`value`, :confval:`min`, :confval:`max`] or
   [:confval:`value`, :confval:`min`, :confval:`max`, :confval:`brute_step`].
   Note that :confval:`brute_step` corresponds to the step size between
   :confval:`min` and :confval:`max`, it is only required for brute-force
   grid search, which also includes :ref:`χ2 surface plot generation
   <additional_chi2>`. See the parameter file for :file:`RELAXATION_NZ/`
   under :file:`Examples/Experiments/` as an actual example to learn how
   to set parameter boundary.


Method files
------------

The method file (indicated with ``-m``) contains the fitting methods to
be used during the fitting process, which typically looks like this:

.. literalinclude:: _static/method.toml
   :language: toml

The method file defines the behavior of each parameter during the
fitting process and may contain multiple fitting steps. Each parameter
typically can be fixed (i.e. ``"fix"``) or fitted (i.e. ``"fit"``) during
each fitting step. If the method file is not provided or certain parameter
not specified, a default behavior will be assigned according to each
experiment module.

.. note::
   If the fitting method contains multiple fitting steps, the behavior of
   each parameter always inherits from the previous fitting step if not
   set in the current step.

.. note::
   The :confval:`INCLUDE` key in method file allows selecting a subset
   of residues for analysis during each fitting step. The residue name
   should match those in experiment files, ``"all"`` is the default value
   which indicates all residues to be included in the current fitting
   step.

.. tip::
   The set of residues to be included in global fits should be chosen
   carefully.  A commonly used multi-step fitting strategy is to select
   a subset of residues with relatively large CPMG dispersion or good
   quality CEST minor dips to obtain global parameters
   (p\ :sub:`b`, k\ :sub:`ex`) first, and then carry out single-residue
   fits with (p\ :sub:`b`, k\ :sub:`ex`) fixed to obtain residue-specific
   parameters (e.g. Δϖ) in the next step. In CPMG experiments, in order
   to get reasonable initial estimates of Δϖ for each residue, an
   additional single-residue fitting step can be carried out at the
   very beginning, see the method file for :file:`CPMG_CH3_1H_SQ/` under
   :file:`Examples/Experiments/` for such an example.

.. tip::
   In method file each parameter not only can be set to ``"fix"``
   or ``"fit"``, but also as constraints which can be any kind of
   relationship with other parameters. See the method file for
   :file:`DCEST_15N_HD_EXCH/` under :file:`Examples/Experiments/`
   for such an example.


Kinetic models
--------------

The kinetic model (indicated with ``-d``) indicates the type of exchange
model to be used for the data analysis, which can be one of the following:

====================  =============================================================================
``2st``               2-state exchange model (default)
``3st``               3-state exchange model
``4st``               4-state exchange model
``2st_rs``            2-state exchange model for residue-specific study
``2st_hd``            2-state exchange model for H/D solvent exchange study
``2st_eyring``        2-state exchange model for temperature-dependent study
``3st_eyring``        3-state exchange model for temperature-dependent study
``2st_binding``       2-state exchange model for ligand binding study
``4st_hd``            4-state exchange model for simultaneous normal and H/D solvent exchange study
====================  =============================================================================

Multiple states involved in exchange process are distinguished with
different parameter suffix :confval:`A`, :confval:`B`, :confval:`C`,
:confval:`D`, etc. For example, :confval:`R1_A` represents R\ :sub:`1`
rate of the major state (i.e.  ground state), :confval:`R2_B`
represents R\ :sub:`2` rate of the first minor state, and so on.

.. note::
   For each kinetic model, it is possible to add ``.mf`` suffix to create
   a new kinetic model (e.g. ``2st.mf``). In such kinetic models,
   model-free parameters (e.g. :confval:`TAUC_A`, :confval:`S2_A`, etc.)
   are directly fitted instead of each individual relaxation parameter
   (e.g. :confval:`R1_A`, :confval:`R2_A`, etc.). See
   :file:`CEST_15N_TR/` under :file:`Examples/Experiments/`
   as an example.


Output files
------------

The fitting output typically contains the following files and directories:

* ``statistics.toml``: Contain all statsitical information about the
  fitting results, such as *χ*\ :sup:`2`.

* ``Data/``: Contain all data files used for the fitting, which are
  identical to the input data files.

* ``Parameters/``: Contain fitting results as three separate files
  ``fitted.toml``, ``fixed.toml`` and ``constrained.toml``, which contain
  output parameters that are fitted, fixed and constrained during the
  fitting process, respectively.

* ``Plot/``: Contain fitting results as plots (in ``.pdf`` format) and also
  the raw datasets (both the original input and fitted data points) for
  creating the plots. Example fitting results for CPMG and CEST
  experiments are shown below:

   .. figure:: _static/cpmg_800mhz_fit.*
      :scale: 60
      :align: center
      :alt: CPMG fitting results
      :figclass: align-center

      Example of CPMG fitting results

   .. figure:: _static/cest_26hz_fit.*
      :scale: 60
      :align: center
      :alt: CEST fitting results
      :figclass: align-center

      Example of CEST fitting results

.. note::
   In plots of CEST/D-CEST/COS-CEST fitting results, the positions for
   ground and excited states are indicated by solid and dashed vertical
   lines respectively, besides, data points that are filtered out
   from the fitting are shown in light color. In plots of
   D-CEST/COS-CEST fitting results, the "folded" positions for ground
   and excited states are indicated by "*" at the vertical lines.

If the fittng method has multiple fitting steps, each step will create
its own output subdirectory.  During any fitting step, if all global
parameters (p\ :sub:`b`, k\ :sub:`ex`, etc.) are fixed or
residue-specific fitting model is used (e.g. ``2st_rs``), then
two separate subdirectories ``All/`` and ``Cluster/`` will be created,
which contain fitting results for all residues and each individual
residue, respectively.
