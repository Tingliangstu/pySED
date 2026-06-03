output_partial
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_partial = 1

**Meaning**
   Controls whether pySED writes atom-type and Cartesian-direction resolved SED
   files. Enable it when you want to plot contributions such as one atom type
   in the ``x``, ``y``, or ``z`` direction.

**Default**
   ``0``.

**Outputs**
   pySED writes partial files under ``<out_files_name>_partial_SED/`` with names
   such as ``<out_files_name>.SED_type1_x``.

**Notes**
   Enable this in compute mode. Plot selected partial components later with
   ``plot_partial_SED``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`plot_partial_SED <plot_partial_SED>`
- :doc:`plot_SED <plot_SED>`
- :doc:`out_files_name <out_files_name>`
- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`file_format <file_format>`
- :doc:`dump_xyz_file <dump_xyz_file>`
