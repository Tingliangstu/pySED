SED Workflow and Trajectories
=============================

Workflow
--------

pySED calculates an eigenvector-free, kinetic-energy-weighted SED:

.. code-block:: text

   trajectory + basis.in + primitive cell + Q path
                         |
                         v
               block-averaged SED
                         |
              +----------+-----------+
              |                      |
              v                      v
      dispersion and slices     partial SED
              |
              v
      peak detection and line-shape fitting
              |
              v
       frequency and lifetime outputs

Recommended order:

1. Validate sampling metadata and the atom-to-basis mapping.
2. Compute with ``action = compute``.
3. Plot the dispersion and representative single-Q spectra.
4. Tune peak detection and fitting on selected Q points.
5. Enable all-Q fitting only after representative fits are satisfactory.

Trajectory formats
------------------

``trajectory_file`` accepts:

- GPUMD extended XYZ with lattice, positions, and velocities;
- a combined LAMMPS custom dump with box bounds and
  ``id type x y z vx vy vz``;
- compatible GPUMD, LAMMPS, or pySED NetCDF.

For LAMMPS text output, use one sorted dump:

.. code-block:: text

   variable dump_stride equal 20
   dump pysed all custom ${dump_stride} trajectory.dump id type x y z vx vy vz
   dump_modify pysed sort id

Set ``lammps_unit = metal`` for velocities in Angstrom/ps or ``real`` for
Angstrom/fs. The unit style is part of cache identity, so changing it invalidates
a previously converted cache.

For GPUMD text output, request velocity fields:

.. code-block:: text

   dump_xyz -1 1 10 dump.xyz velocity

The reader follows the extended-XYZ ``Properties=`` schema rather than
assuming a fixed field order.

Cache, direct reading, and NetCDF
---------------------------------

For text trajectories, the default

.. code-block:: text

   trajectory_read_mode = cache

streams the source once into a reusable ``.pysed.nc`` cache next to it.
The cache fingerprint includes source identity and conversion-affecting
settings. ``trajectory_read_mode = direct`` consumes text sequentially
without an intermediate file. Native NetCDF input is always read directly.

``trajectory_prefetch = 1`` prepares at most one following block in a
background thread. Disable it when the extra raw block would exceed memory.
``netcdf_batch_size`` controls text-conversion batches and
``netcdf_compression_level`` selects lossless zlib compression from 0 to 9.

Sampling consistency
--------------------

The requested frame count is

.. math::

   N_\mathrm{frames}
   =
   \frac{\mathtt{total\_num\_steps}}
        {\mathtt{output\_data\_stride}}.

It must be at least ``num_blocks`` and divisible by ``num_blocks``.
``time_step`` is the MD integration step in fs; the actual sample spacing is
``time_step * output_data_stride``.

Basis mapping and local groups
------------------------------

``basis.in`` maps every trajectory atom to a repeated cell, basis atom, and
mass:

.. code-block:: text

   atoms_ids unitcell_index basis_index mass_types
   1  1  1  28.0855
   2  2  1  28.0855

Atom IDs and count must match the trajectory. Each repeated cell must contain
the expected basis indices. ``supercell_dim`` describes repetitions used for
commensurate Q-point construction.

When the trajectory contains a periodic local group rather than the full MD
cell, set ``dump_group_supercell`` to that group's three lattice vectors.
It changes the cell used for reciprocal-space construction; it does not alter
coordinates or the ``basis.in`` mapping.

Backends and partial SED
------------------------

``backend = numpy`` uses one process at ``max_cores = 1`` and persistent
worker processes otherwise. ``backend = cupy`` uses one CUDA GPU for the SED
kernel; plotting and fitting remain CPU operations.

Set ``output_partial = 1`` during computation to preserve element- and
Cartesian-resolved components. Use ``plot_partial_SED = Si`` for an element
sum or ``plot_partial_SED = Si x`` for one direction.
