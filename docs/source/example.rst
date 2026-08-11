Examples
========

All SED example control files use ``input.in`` and the pySED 3.0 action
workflow. Generated numerical files and figures are intentionally not stored
beside active examples; run ``compute``, then ``plot`` or ``fit``.

Recommended examples
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30

   * - Folder
     - Trajectory
     - System
     - Main purpose
   * - ``example/deterministic_sed``
     - generated LAMMPS custom dump
     - 8-atom synthetic crystal
     - Fast parser, cache/direct, action, and known-frequency validation.
   * - ``example/CNT``
     - GPUMD extxyz
     - carbon nanotube
     - One-dimensional G-A dispersion and all-Q fitting.
   * - ``example/In_plane_graphene_gpumd``
     - GPUMD extxyz
     - graphene
     - Exact fractional G-M-K-G path and contour plotting.
   * - ``example/MoS2_gpumd``
     - GPUMD extxyz
     - bulk MoS2
     - Low-frequency out-of-plane modes and fitting.
   * - ``example/Silicon_primitive_gpumd``
     - GPUMD extxyz
     - bulk silicon
     - Multi-segment three-dimensional Q path.

The four cases under ``example/For_old_version_example`` have also been
migrated. They demonstrate combined LAMMPS custom dumps for Aluminum FCC,
high-entropy oxide, PTCDI-C8, and conventional-cell silicon. Their directory
name records their origin; their active ``input.in`` files are v3 inputs.

Run the deterministic case
--------------------------

.. code-block:: bash

   cd example/deterministic_sed
   python generate_trajectory.py
   pySED input.in
   pySED input-plot.in
   pySED input-fit.in

The known X-point signal is centered at 6.25 THz. ``expected.json`` records
the frame count, Q point, frequency resolution, and cache/direct comparison
tolerance. Use ``input-direct.in`` for the direct-reader calculation.

Run a GPUMD case
----------------

Each GPUMD material folder has ``structure/``, ``gpumd_run/``, and ``SED/``
subdirectories:

.. code-block:: bash

   cd example/CNT/structure
   python generate_gpumd_xyz.py

   cd ../gpumd_run
   gpumd

   cd ../SED
   pySED input.in

Verify that the GPUMD output stride equals ``output_data_stride``. The text
trajectory must include lattice, positions, and velocities.

Run a LAMMPS case
-----------------

The migrated inputs use one custom dump:

.. code-block:: text

   dump pysed all custom ${dump_stride} trajectory.dump id type x y z vx vy vz
   dump_modify pysed sort id

Generate or copy the matching structure and ``basis.in``, run
``lammps_run/in.lammps``, then run ``pySED input.in`` from the ``SED``
directory. Separate ``pos.dat`` and ``vels.dat`` are not v3 inputs.

Historical outputs
------------------

Exactly 81 generated pySED 2.4 example files are preserved at
``example/legacy/pySED-2.4/outputs``. The archive retains original relative
paths and includes ``MANIFEST.tsv`` and ``SHA256SUMS``. These files use
the old interface and normalization and are not expected pySED 3.0 results.

The images under ``Ref_Phonon_dispersion_from_phonopy`` are harmonic
reference calculations. They must not be confused with finite-temperature SED
output.
