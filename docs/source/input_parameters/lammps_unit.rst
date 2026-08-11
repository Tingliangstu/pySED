lammps_unit
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lammps_unit = 'metal'

**Meaning**
   LAMMPS unit style used to convert velocities into the units expected by
   pySED. Set this to match the ``units`` command in the LAMMPS simulation.

**Default**
   ``metal``.

**Allowed values**
   ``metal`` or ``real``.

**Notes**
   ``metal`` velocities are interpreted as Angstrom/ps. ``real`` velocities are
   interpreted as Angstrom/fs.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`out_files_name <out_files_name>`
- :doc:`basis_lattice_file <basis_lattice_file>`
