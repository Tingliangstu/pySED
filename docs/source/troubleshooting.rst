Troubleshooting
===============

If you have any questions, please email `liangting.zj@gmail.com <mailto:liangting.zj@gmail.com>`_ or `wwjiang96@163.com <mailto:wwjiang96@163.com>`_.

Input is rejected
-----------------

Run ``pySED input.in`` and read the first reported parameter error. pySED
3.0 intentionally rejects v2.4 names; use :doc:`migration_v3` rather than
mixing old and new keys.

Trajectory cannot be read
-------------------------

For LAMMPS text, verify that every frame has box bounds and
``id type x y z vx vy vz``, and that IDs are stable. For GPUMD extxyz,
verify that ``Properties=`` declares position and velocity fields. Confirm
that ``num_atoms`` agrees with both the trajectory and ``basis.in``.

Frame/block mismatch
--------------------

The requested frame count is
``total_num_steps / output_data_stride``. It must be divisible by
``num_blocks``. Correct the metadata or select another block count; do not
silently discard frames.

Plot or fit reports missing files
---------------------------------

Run ``action = compute`` first and keep ``out_files_name`` unchanged.
Plot and fit now fail nonzero rather than reporting success without their
numerical prerequisites.

Cache appears stale
-------------------

pySED fingerprints the trajectory and conversion settings, including
``lammps_unit``. If a cache was copied or manually renamed, remove the
``.pysed.nc`` file and let pySED recreate it.

GPU backend fails
-----------------

Confirm that the installed CuPy wheel matches the CUDA Toolkit and driver.
Switch to ``device = cpu`` to distinguish an input/data problem from a
CUDA installation problem.

Further help
------------

Report reproducible issues at
`GitHub <https://github.com/AI4MD-lab/pySED/issues>`_ and include the
input file, pySED version, Python version, traceback, and a minimal trajectory
when sharing is permitted.
