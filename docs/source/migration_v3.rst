Migrating from pySED 2.4
==========================

pySED 3.0 is a breaking upgrade. Rename ``input_SED.in`` to ``input.in``,
use one trajectory file, and select work through ``action``.

Parameter replacements
----------------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Removed v2.4 parameter
     - v3 replacement
   * - ``plot_SED``, ``lorentz``
     - ``action = compute | plot | fit | thinking``
   * - ``file_format``
     - Automatic detection from ``trajectory_file``.
   * - ``dump_xyz_file``, ``pos_file``, ``vels_file``
     - One ``trajectory_file``; LAMMPS must combine positions and velocities.
   * - ``num_splits``
     - ``num_blocks``.
   * - ``compress``, ``output_hdf5``
     - ``trajectory_read_mode`` and an automatic ``.pysed.nc`` cache.
   * - ``use_parallel``
     - ``max_cores = 1`` for serial or a larger value for multiprocessing.
   * - ``peak_height``, ``peak_prominence``
     - Dimensionless ``peak_min_significance``.
   * - ``lorentz_fit_cutoff``
     - ``lorentz_fit_freq_min`` and ``lorentz_fit_freq_max``.
   * - ``with_eigs``
     - Removed without replacement.

Old parameter names are rejected with migration guidance. SED output is
eV/THz by default. The temporary :doc:`legacy_unit <input_parameters/legacy_unit>`
option restores only the historical intensity scale for comparisons;
it does not restore old input names or fitting behavior. Set it to ``True``
when reading 2.x ``.SED`` files without scale metadata. This option will be
removed in a future major release.

Q-path construction
-------------------

The current GCD/Bézout solver keeps only supercell-compatible points and uses
consistent rational endpoints for selection and output. It also checks
zero-length paths and rejects non-integer cell transformations using an
absolute tolerance. The input syntax is unchanged; explicit fractions such
as ``1/3`` and ``7/15`` retain their exact values. See
:doc:`input_parameters/q_path` for decimal normalization and :doc:`theory`
for the periodicity condition.

LAMMPS migration
----------------

Replace two dumps such as ``pos.dat`` and ``vels.dat`` with one:

.. code-block:: text

   dump pysed all custom ${dump_stride} trajectory.dump id type x y z vx vy vz
   dump_modify pysed sort id

Then set ``trajectory_file = trajectory.dump`` and the appropriate
``lammps_unit``.
