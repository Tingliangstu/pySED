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

- :doc:`trajectory_file <trajectory_file>`
- :doc:`out_files_name <out_files_name>`
- :doc:`lammps_unit <lammps_unit>`
