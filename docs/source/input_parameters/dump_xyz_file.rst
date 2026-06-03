dump_xyz_file
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   dump_xyz_file = '../gpumd_run/dump.xyz'

**Meaning**
   Path to the GPUMD extended XYZ trajectory file, usually ``dump.xyz``. This
   file must contain the positions and velocities saved from the MD production
   run.

**Default**
   ``dump.xyz``.

**Notes**
   This is used only when ``file_format = 'gpumd'``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`out_files_name <out_files_name>`
- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`file_format <file_format>`
- :doc:`pos_file <pos_file>`
- :doc:`vels_file <vels_file>`
- :doc:`lammps_unit <lammps_unit>`
