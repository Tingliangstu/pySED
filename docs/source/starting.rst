Quick Start
===========

Check the installation
----------------------

.. code-block:: bash

   pySED -h

The lowercase aliases ``pysed`` and ``PYSED`` are also installed, but the
manual uses ``pySED`` consistently.

Create an input file
--------------------

Save the following as ``input.in`` next to a GPUMD extended-XYZ trajectory.
Adapt the atom count, sampling, structure, and Q path to the simulation:

.. code-block:: text

   action  = compute
   method  = sed
   backend = numpy

   trajectory_file       = ../gpumd_run/dump.xyz
   trajectory_read_mode  = cache
   trajectory_prefetch   = 1
   out_files_name        = silicon
   lammps_unit           = metal

   num_atoms             = 16000
   total_num_steps       = 100000
   time_step             = 1.0
   output_data_stride    = 20
   num_blocks            = 5
   max_cores             = 4

   basis_lattice_file    = ../structure/basis.in
   prim_unitcell         = 3.86697 0 0  1.93349 3.34890 0  1.93349 1.11630 3.15737
   supercell_dim         = 20 20 20
   rescale_prim          = 1

   num_qpaths            = 5
   q_path_name           = GXUKGL
   q_path                = 0 0 0  1/2 0 1/2  5/8 1/4 5/8  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2

   plot_cutoff_freq      = 20
   plot_interval         = 5
   plot_slice            = 1
   qpoint_slice_index    = 2
   if_show_figures       = 0

Run the stages
--------------

.. code-block:: bash

   pySED input.in

Use ``action = compute`` to calculate the numerical SED, ``action = plot``
to regenerate plots, and ``action = fit`` to detect and fit spectral peaks.
``action = thinking`` reuses outputs already present and selects the next
missing stage.

A failed or inapplicable stage exits nonzero. In particular, plot and fit do
not report success when their required numerical outputs are absent.

Prepare the trajectory
----------------------

pySED accepts GPUMD extended XYZ, one combined LAMMPS custom dump, and
compatible NetCDF trajectories. A LAMMPS text trajectory must contain atom IDs,
types, positions, velocities, and box information:

.. code-block:: text

   dump pysed all custom ${dump_stride} trajectory.dump id type x y z vx vy vz
   dump_modify pysed sort id

Ordinary LAMMPS XYZ output and separate position/velocity dumps are not pySED
3.0 inputs. See :doc:`sed_workflow` for trajectory and ``basis.in``
requirements, and :doc:`example` for complete cases.
