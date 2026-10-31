Quick Start
===========

The pySED workflow is a post-processing workflow. First run molecular dynamics
to produce a trajectory with coordinates and velocities, then run pySED to
prepare the trajectory, construct commensurate q-points, calculate SED, plot
the result, and optionally fit Lorentzian peaks.

Recommended Workflow
--------------------

1. **Prepare the primitive structure.**
   Start from a POSCAR or compatible structure file for the primitive cell.

2. **Generate the MD supercell and ``basis.in``.**
   Use ``pySED.structure.generate_data.structure_maker``. The generated
   ``basis.in`` maps every atom in the MD supercell to its unit-cell index,
   basis index, and mass. pySED needs this mapping for reciprocal-space SED.

3. **Run the MD simulation.**
   Use GPUMD or LAMMPS. Equilibrate first, then run an NVE production trajectory
   for SED.

4. **Edit ``input.in``.**
   Set the trajectory path, number of atoms, MD step information, primitive
   cell, supercell dimensions, and q-path. Choose
   ``backend = numpy`` or ``cupy``. The complete control file is provided
   in each example's ``SED/input.in``.

5. **Compute SED.**

   .. code-block:: bash

      cd example/MoS2_gpumd/SED
      pysed input.in

   For the first run, use ``action = compute``. pySED writes ``.SED``, ``.Qpts``,
   ``.THz``, and ``.Q_distances_and_labels`` files.

6. **Plot SED.**
   Set ``action = plot`` and run pySED again. Tune ``plot_cutoff_freq``,
   ``plot_interval``, ``plot_color``, ``colorbar_min``, and ``colorbar_max`` for
   a clean figure.

7. **Fit lifetimes if needed.**
   First set ``plot_slice = 1`` and choose ``qpoint_slice_index`` to inspect a
   single q-point with ``action = plot``. Then use ``action = fit`` and
   tune ``peak_min_significance``. For Lorentzian fits, set
   ``fitting_function = lorentz``. After a good single-q fit, set
   ``lorentz_fit_all_qpoint = 1``.

GPUMD Trajectory
----------------

For GPUMD users, the production run should write an extended XYZ trajectory
with positions and velocities. A minimal production block is:

.. code-block:: text

   ensemble       nve
   dump_exyz      10     1
   run            500000

The first number after ``dump_exyz`` is the output stride in MD steps. Use the
same value for ``output_data_stride`` in ``input.in``. For details, see the
GPUMD manual page for
`dump_exyz <https://gpumd.org/gpumd/input_parameters/dump_exyz.html>`_.

In ``input.in``, use:

.. code-block:: text

   trajectory_file = '../gpumd_run/dump.xyz'

LAMMPS Trajectory
-----------------

LAMMPS writes positions and velocities together in one custom dump for
pySED 3.0. The required dump format is sorted by atom id:

.. code-block:: text

   dump            pysed  all  custom  ${dt_dump}  trajectory.dump  id type x y z vx vy vz
   dump_modify     pysed  sort  id

   run             2097152

In ``input.in``, use:

.. code-block:: text

   trajectory_file = '../lammps_run/trajectory.dump'
   lammps_unit = 'metal'

Use ``lammps_unit = 'metal'`` for velocities in Angstrom/ps and
``lammps_unit = 'real'`` for velocities in Angstrom/fs.

What to Check Before Running pySED
----------------------------------

- ``num_atoms`` equals the number of atoms in the trajectory and the maximum
  atom id in ``basis.in``.
- ``total_num_steps``, ``time_step``, and ``output_data_stride`` match the MD
  production run.
- ``prim_unitcell`` and ``supercell_dim`` reconstruct the MD supercell.
- ``basis_lattice_file`` points to the correct ``basis.in`` file.
- ``q_path_name`` has ``num_qpaths + 1`` labels.
- ``q_path`` contains ``num_qpaths + 1`` reduced-coordinate triples.
- The trajectory contains coordinates and velocities for every saved frame.

Recommended First Examples
--------------------------

Start with an existing example before applying pySED to a new system:

- 1D: `CNT <https://github.com/AI4MD-lab/pySED/tree/codex/pysed-3.0-architecture/example/CNT>`_
- 2D: `Graphene <https://github.com/AI4MD-lab/pySED/tree/codex/pysed-3.0-architecture/example/In_plane_graphene_gpumd>`_
  or `MoS₂ <https://github.com/AI4MD-lab/pySED/tree/codex/pysed-3.0-architecture/example/MoS2_gpumd>`_
- 3D: `Silicon <https://github.com/AI4MD-lab/pySED/tree/codex/pysed-3.0-architecture/example/Silicon_primitive_gpumd>`_

The illustrated 2.4 results show the expected dispersion features. Recompute
the data with 3.0 before comparing intensities or fitting lifetimes; the
normalization and fitting controls changed. A small executable check is
available at :doc:`sed_workflow`.
