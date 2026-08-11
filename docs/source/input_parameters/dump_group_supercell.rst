dump_group_supercell
~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   dump_group_supercell = 80.712 0 0  40.356 69.899 0  0 0 12.356

**Meaning**
   The local cell represented by ``basis.in`` when the trajectory contains
   atoms dumped from one group of a larger simulation.

**Default**
   ``None``. pySED uses the simulation cell stored in the trajectory.

When to use it
--------------

Use this parameter for local SED after dumping one group from a larger model.
For example, GPUMD can dump a pure-material region from a heterostructure with
``dump_netcdf <grouping_method> <group_id> ...``. The trajectory then contains
only that group's atoms, but its NetCDF cell still describes the complete
simulation box.

The atoms are correct; only the cell used to choose the local Q points needs
to be changed. Set ``dump_group_supercell`` to the cell of the local region
described by the matching ``basis.in``. pySED uses it to make the Q-point
grid and leaves atom positions and velocities unchanged.

Do not use this parameter for a whole-system trajectory whose ``basis.in``
already describes the full simulation cell.

**Notes**
   Write the three local lattice vectors as nine numbers in Angstrom. They
   must use the same Cartesian coordinate frame as the trajectory and
   ``prim_unitcell``.

   Set ``rescale_prim = 0`` and provide the matching ``prim_unitcell``. The
   program checks that the local cell and ``basis.in`` describe the same number
   of primitive cells. If they do not match, use the ``basis.in`` and local
   cell generated for the same dump group.

Example: a local diamond group
------------------------------

.. code-block:: text

   trajectory_file       = trajectory_group_09.nc
   basis_lattice_file     = group_09/basis.in
   prim_unitcell          = 2.52225 0 0  1.26112 2.18433 0  -1.26112 -0.72811 2.05941
   dump_group_supercell   = 80.711996 0 0  40.355998 69.898639 0  0 0 12.356450
   rescale_prim           = 0

Each local group needs its own matching ``basis.in``. Groups with a different
material or a different thickness can require a different
``dump_group_supercell``.

:doc:`Back to Parameter Index <../input_parameters>`

Related parameters
------------------

- :doc:`basis_lattice_file <basis_lattice_file>`
- :doc:`prim_unitcell <prim_unitcell>`
- :doc:`supercell_dim <supercell_dim>`
- :doc:`rescale_prim <rescale_prim>`
- :doc:`q_path <q_path>`
