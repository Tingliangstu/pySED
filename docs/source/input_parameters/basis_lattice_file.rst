basis_lattice_file
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   basis_lattice_file = '../structure/basis.in'

**Meaning**
   Path to ``basis.in``, the file that maps atoms in the MD supercell to the
   primitive-cell basis used by pySED. pySED needs this mapping to connect the
   trajectory atoms with commensurate q-points.

**Default**
   ``basis.in``.

**Notes**
   The file maps each atom id to a unit-cell index, basis index, and mass. It is
   required for q-point construction and SED projection.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`out_files_name <out_files_name>`
- :doc:`file_format <file_format>`
- :doc:`dump_xyz_file <dump_xyz_file>`
- :doc:`pos_file <pos_file>`
- :doc:`vels_file <vels_file>`
- :doc:`lammps_unit <lammps_unit>`
