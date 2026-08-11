output_partial
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_partial = 1

**Meaning**
   Controls whether pySED writes element- and Cartesian-direction-resolved SED
   files. Enable it when you want to plot contributions such as one element
   in the ``x``, ``y``, or ``z`` direction.

**Default**
   ``0``.

**Outputs**
   pySED writes partial files under ``<out_files_name>_partial_SED/`` with names
   such as ``<out_files_name>.SED_O_x``. Their unit is ``eV/THz``.

**Notes**
   Enable this in compute mode. Plot selected partial components later with
   ``plot_partial_SED``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_partial_SED <plot_partial_SED>`
- :doc:`out_files_name <out_files_name>`
- :doc:`basis_lattice_file <basis_lattice_file>`
