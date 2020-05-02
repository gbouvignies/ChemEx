.. _chemex_additional:

====================
Additional Functions
====================

Aside from fitting input datasets to obtain chemical exchange parameters, 
ChemEx also has many other useful functions, several most important ones 
are summarized as below.


Simulate CPMG or CEST profiles based on a given sets of input parameters
------------------------------------------------------------------------

ChemEx allows simulating CPMG or CEST profiles based on a given sets of 
input parameters, which is useful to learn about the effects of each 
individual parameter on the final results. In order to make use of this 
function the first argument should be set to ``simulate``. Since all 
parameters are fixed to the initial input value, it is not necessary to 
provide the method file.  A typical command for simulation purpose with 
ChemEx is like this: 

.. code-block:: console

   $ chemex simulate -e <FILE> \
                     -p <FILE> \
                     -d <MODEL> \
                     -o <DIR>

Example simulation results for CPMG and CEST experiments are shown below:
                      
   .. figure:: _static/cpmg_800mhz_simu.*
      :scale: 60
      :align: center 
      :alt: CPMG simulation results
      :figclass: align-center 

      Example of CPMG simulation results

   .. figure:: _static/cest_26hz_simu.*
      :scale: 60
      :align: center 
      :alt: CEST simulation results
      :figclass: align-center 

      Example of CEST simulation results


Obtain initial estimates of Δϖ for CEST experiments
---------------------------------------------------------------------

In CEST experiments, in order to avoid getting trapped in a local
minimum, it is necessary to choose suitable initial value of Δϖ.
ChemEx comes with a module :ref:`tools_pick_cest` for manually 
picking the major and minor dips of CEST profiles, which correspond
to the ground and excited states, respectively. A typical command
for such purpose is like this:

.. code-block:: console

   $ chemex pick_cest -e <FILE> -o <DIR>

After typing this command, a window showing all CEST profiles will 
come out. For each profile first click on the major dip and then on the 
minor dip. Note that in certain profiles only one dip could be visible, 
which indicates the minor dip is overlapped with the major dip, therefore 
the major dip should be clicked twice. When done with any profile, click 
the :guilabel:`Next` or :guilabel:`Previous` button to proceed to the 
next or previous profile. The :guilabel:`Clear` button allows cleaning
the selection in the current profile. After all profiles are finished 
click the :guilabel:`Quit` button, then two separate files will be 
created: :file:`cs_a.toml` and :file:`dw_ab.toml` that contain chemical 
shifts of the major state and chemical shift difference between the major 
and minor states, respectively. Try to run the :file:`pick_cest.sh` script 
under :file:`CEST_15N/` example to learn how to make use of this function.


.. _additional_visualize:

Visualize fitting results interactively
---------------------------------------

ChemEx comes with a module :ref:`tools_plot_param` that allows 
visualizing the fitting results interactively, a typical command 
for such purpose is like this:

.. code-block:: console

   $ chemex plot_param -p <FILE> -n <NAME>

See :file:`2stBinding/` example to learn how to make use of this function.
After finish running :file:`run.sh`, the chemical shift differences between
the free and bound states can be displayed with:

.. code-block:: console

   $ chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n DW_AB

and the transverse relaxation rates of both states can be compared with:

.. code-block:: console

   $ chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n R2

These two commands are saved in the :file:`plot_param.sh` script in 
:ref:`this example <example_binding>`. From these two observables, 
the core region of the interaction site can be clearly located. Aside 
from the core region, there is also a tail with increased R\ :sub:`2` 
rates located at C-terminal end of the interaction site and with very 
little chemical shift perturbation. This region is likely involved 
in the transient interactions with the binding partner, which 
causes certain degree of steric restriction to this region.


.. _additional_chi2:

Create *χ*\ :sup:`2` surface plots for CPMG or CEST experiments
----------------------------------------------------------------

*χ*\ :sup:`2` surface plot is commonly used for showing the 
dependence of *χ*\ :sup:`2` on each parameter. In order to calculate
*χ*\ :sup:`2` surface map, a grid set of parameters should be chosen.
A commonly used scheme is to calculate the dependence of 
*χ*\ :sup:`2` on p\ :sub:`b` and  k\ :sub:`ex`, besides, 
one-dimensional *χ*\ :sup:`2` surface plot can be created
based on the dependence on each individual parameter. With the 
:ref:`tools_chi2_surface` module in ChemEx, *χ*\ :sup:`2` surface 
plot can be easily created. A typical command for such purpose is 
like this:

.. code-block:: console

   $ chemex chi2_surface -e <FILE> \
                         -p <FILE> \
                         -m <FILE> \
                         -d <MODEL> \
                         -o <DIR>

Try to run the :file:`chi2_surface.sh` script in :file:`CPMG_15N_IP/` 
example to learn how to create *χ*\ :sup:`2` surface plots, which include 
both one- and two-dimensional examples.

