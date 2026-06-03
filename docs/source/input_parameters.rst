Input Parameters
================

pySED reads a plain-text control file, usually named ``input_SED.in``. Each
active line has the form

.. code-block:: text

   parameter_name = value

Blank lines and text after ``#`` are ignored. If no input file is passed on the
command line, pySED reads ``input_SED.in`` in the current directory:

.. code-block:: bash

   pysed
   pysed input_SED.in
   pySED input_SED.in

This page documents the parameters accepted in ``input_SED.in`` by the current
parser. Use the table below as a quick index; clicking a parameter name jumps to
the detailed explanation.

Parameter Index
---------------

.. list-table:: ``input_SED.in`` parameters
   :header-rows: 1
   :widths: 24 24 52

   * - Parameter
     - Group
     - Purpose
   * - `num_atoms`_
     - MD simulation
     - Number of atoms in the MD trajectory.
   * - `total_num_steps`_
     - MD simulation
     - Number of MD steps used for the production trajectory.
   * - `time_step`_
     - MD simulation
     - MD time step in fs.
   * - `output_data_stride`_
     - MD simulation
     - Trajectory output interval in MD steps.
   * - `num_splits`_
     - Averaging
     - Number of trajectory blocks used for averaging.
   * - `compress`_
     - Averaging
     - Compress trajectory data to HDF5 before SED.
   * - `out_files_name`_
     - Files
     - Prefix for pySED output files.
   * - `basis_lattice_file`_
     - Files
     - Path to ``basis.in``.
   * - `file_format`_
     - Files
     - Select GPUMD or LAMMPS trajectory reader.
   * - `dump_xyz_file`_
     - Files
     - GPUMD extended XYZ trajectory path.
   * - `pos_file`_
     - Files
     - LAMMPS position trajectory path.
   * - `vels_file`_
     - Files
     - LAMMPS velocity trajectory path.
   * - `lammps_unit`_
     - Files
     - LAMMPS unit style for velocity conversion.
   * - `output_hdf5`_
     - Files
     - Compressed trajectory database name.
   * - `output_partial`_
     - Files
     - Write atom-type and direction-resolved SED files.
   * - `use_parallel`_
     - Parallelization
     - Enable multiprocessing over q-points.
   * - `max_cores`_
     - Parallelization
     - Maximum number of worker processes.
   * - `prim_unitcell`_
     - Crystal structure
     - Primitive lattice vectors.
   * - `prim_axis`_
     - Crystal structure
     - Primitive-axis transformation.
   * - `supercell_dim`_
     - Crystal structure
     - Supercell repeats; controls q-point resolution.
   * - `rescale_prim`_
     - Crystal structure
     - Reconstruct the primitive cell from the trajectory cell.
   * - `num_qpaths`_
     - Q-path
     - Number of high-symmetry path segments.
   * - `q_path_name`_
     - Q-path
     - Labels for high-symmetry points.
   * - `q_path`_
     - Q-path
     - High-symmetry q-points in reduced coordinates.
   * - `plot_SED`_
     - Plotting
     - Select compute mode or plot/fitting mode.
   * - `plot_cutoff_freq`_
     - Plotting
     - Maximum plotted frequency in THz.
   * - `plot_interval`_
     - Plotting
     - Frequency tick interval.
   * - `plot_color`_
     - Plotting
     - Matplotlib colormap.
   * - `colorbar_min`_
     - Plotting
     - Lower log-scale colorbar limit.
   * - `colorbar_max`_
     - Plotting
     - Upper log-scale colorbar limit.
   * - `use_contourf`_
     - Plotting
     - Use contour-style SED plotting.
   * - `qpoint_slice_index`_
     - Plotting
     - Zero-based q-point index for slice plots and fitting.
   * - `plot_slice`_
     - Plotting
     - Plot one q-point spectrum.
   * - `plot_partial_SED`_
     - Plotting
     - Plot one atom-type or direction-resolved partial SED.
   * - `if_show_figures`_
     - Plotting
     - Show figures interactively.
   * - `lorentz`_
     - Lorentzian fitting
     - Enable Lorentzian peak fitting.
   * - `peak_height`_
     - Lorentzian fitting
     - Minimum detected peak height.
   * - `peak_prominence`_
     - Lorentzian fitting
     - Minimum detected peak prominence.
   * - `initial_guess_hwhm`_
     - Lorentzian fitting
     - Initial HWHM guess.
   * - `peak_max_hwhm`_
     - Lorentzian fitting
     - Maximum fitted HWHM.
   * - `lorentz_fit_cutoff`_
     - Lorentzian fitting
     - Maximum fitting frequency in THz.
   * - `modulate_factor`_
     - Lorentzian fitting
     - Shrink fitting windows around detected peaks.
   * - `lorentz_fit_all_qpoint`_
     - Lorentzian fitting
     - Fit peaks at all q-points.
   * - `re_output_total_freq_lifetime`_
     - Lorentzian fitting
     - Rebuild total frequency-lifetime output after refitting.
   * - `with_eigs`_
     - Reserved
     - Reserved for eigenvector-related development.

MD Simulation Parameters
------------------------

num_atoms
~~~~~~~~~

**Syntax**

.. code-block:: text

   num_atoms = 13824

**Meaning**
   Number of atoms in the MD supercell.

**Default**
   ``0``. This must be set for real calculations.

**Example**
   In the MoS2 GPUMD example, ``num_atoms = 13824``.

**Notes**
   The value must match the maximum atom id in ``basis.in``. pySED stops if
   these two values are inconsistent.

total_num_steps
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   total_num_steps = 500000

**Meaning**
   Total number of MD steps in the production trajectory used for SED.

**Default**
   ``0``. This must be set.

**Notes**
   The number of trajectory frames used by pySED is
   ``total_num_steps / output_data_stride``.

time_step
~~~~~~~~~

**Syntax**

.. code-block:: text

   time_step = 1

**Meaning**
   MD time step in femtoseconds.

**Default**
   ``0``. This must be set.

**Notes**
   Use the time step from the production MD run, not from a separate relaxation
   run.

output_data_stride
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_data_stride = 50

**Meaning**
   Number of MD steps between saved trajectory frames.

**Default**
   ``0``. This must be set.

**Example**
   If GPUMD uses ``dump_exyz 50 1``, set ``output_data_stride = 50``.

**Notes**
   This value determines the frequency grid together with ``time_step`` and
   ``total_num_steps``.

Data Splitting and Compression
------------------------------

num_splits
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_splits = 5

**Meaning**
   Number of trajectory blocks used for block averaging.

**Default**
   ``1``.

**Notes**
   Values such as ``5`` or ``10`` often give smoother SED maps. Each split uses
   fewer frames, so very large values reduce frequency resolution.

compress
~~~~~~~~

**Syntax**

.. code-block:: text

   compress = 1

**Meaning**
   Compress trajectory coordinates and velocities into an HDF5 file before SED
   computation.

**Default**
   ``1``.

**Notes**
   Keep this enabled for normal workflows. If the HDF5 file already exists,
   pySED reuses it.

Input and Output Files
----------------------

out_files_name
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   out_files_name = 'bulk_MoS2'

**Meaning**
   Prefix used for output files.

**Default**
   No parser default. Set it explicitly.

**Outputs**
   ``bulk_MoS2.SED``, ``bulk_MoS2.Qpts``, ``bulk_MoS2.THz``,
   ``bulk_MoS2.Q_distances_and_labels``, and ``bulk_MoS2-SED.png``.

basis_lattice_file
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   basis_lattice_file = '../structure/basis.in'

**Meaning**
   Path to the basis mapping file generated by pySED structure utilities.

**Default**
   ``basis.in``.

**Notes**
   The file maps each atom id to a unit-cell index, basis index, and mass. It is
   required for q-point construction and SED projection.

file_format
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   file_format = 'gpumd'

**Meaning**
   Selects the trajectory reader.

**Default**
   ``gpumd``.

**Allowed values**
   ``gpumd`` or ``lammps``.

dump_xyz_file
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   dump_xyz_file = '../gpumd_run/dump.xyz'

**Meaning**
   GPUMD extended XYZ trajectory file containing positions and velocities.

**Default**
   ``dump.xyz``.

**Notes**
   This is used only when ``file_format = 'gpumd'``.

pos_file
~~~~~~~~

**Syntax**

.. code-block:: text

   pos_file = 'pos.dat'

**Meaning**
   LAMMPS position trajectory.

**Default**
   ``pos.dat``.

**Notes**
   This is used only when ``file_format = 'lammps'``.

vels_file
~~~~~~~~~

**Syntax**

.. code-block:: text

   vels_file = 'vels.dat'

**Meaning**
   LAMMPS velocity trajectory.

**Default**
   ``vels.dat``.

**Notes**
   This is used only when ``file_format = 'lammps'``.

lammps_unit
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lammps_unit = 'metal'

**Meaning**
   LAMMPS unit style used for velocity conversion.

**Default**
   ``metal``.

**Allowed values**
   ``metal`` or ``real``.

**Notes**
   ``metal`` velocities are interpreted as Angstrom/ps. ``real`` velocities are
   interpreted as Angstrom/fs.

output_hdf5
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_hdf5 = 'vel_pos_compress.hdf5'

**Meaning**
   Name of the compressed trajectory database.

**Default**
   ``vel_pos_compress.hdf5``.

output_partial
~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   output_partial = 1

**Meaning**
   Output partial SED components for atom types and Cartesian directions.

**Default**
   ``0``.

**Outputs**
   pySED writes partial files under ``<out_files_name>_partial_SED/`` with names
   such as ``<out_files_name>.SED_type1_x``.

**Notes**
   Enable this in compute mode. Plot selected partial components later with
   ``plot_partial_SED``.

Parallelization
---------------

use_parallel
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   use_parallel = 1

**Meaning**
   Enable multiprocessing over q-points.

**Default**
   ``1``.

**Notes**
   On Windows, or when memory is limited, set ``use_parallel = 0`` or use a small
   ``max_cores``.

max_cores
~~~~~~~~~

**Syntax**

.. code-block:: text

   max_cores = 4

**Meaning**
   Maximum number of worker processes used in parallel mode.

**Default**
   ``4``.

**Notes**
   More cores can reduce compute time but increase memory use because each worker
   handles trajectory data.

Crystal Structure
-----------------

prim_unitcell
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   prim_unitcell = 3.11905 0 0 -1.55952 2.70117 0 0 0 12.6227

**Meaning**
   Primitive-cell lattice vectors in Angstrom, written as nine numbers and
   reshaped into a 3 by 3 matrix.

**Default**
   No useful default. Set it explicitly.

**Notes**
   The rows are lattice vectors. Use the primitive cell consistent with
   ``basis.in`` and ``supercell_dim``.

prim_axis
~~~~~~~~~

**Syntax**

.. code-block:: text

   prim_axis = 0.0 0.5 0.5 0.5 0.0 0.5 0.5 0.5 0.0

**Meaning**
   Transformation from the input cell to a primitive cell, following the same
   idea as phonopy primitive axes.

**Default**
   ``None``.

**Notes**
   Use only when the input structure cell differs from the primitive cell needed
   for q-point construction. For this workflow, keep the MD cell consistent and
   verify the transformed primitive cell printed by pySED.

supercell_dim
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   supercell_dim = 12 12 16

**Meaning**
   Number of primitive-cell repeats used to build the MD supercell.

**Default**
   ``1 1 1``.

**Notes**
   This is one of the most important parameters in ``input_SED.in`` because it
   determines the repetition matrix used for commensurate q-point selection.
   pySED does not use a separate ``num_qpoints`` parameter; the available
   q-points are determined by ``supercell_dim``, the primitive cell, the
   trajectory cell, and ``q_path``.

   For a simple diagonal supercell with repetition matrix
   :math:`P=\mathrm{diag}(N_x,N_y,N_z)`, a path from Gamma to the zone boundary
   along one reduced reciprocal direction, for example
   ``q_path = 0 0 0  0 0 0.5``, gives allowed points
   :math:`0, 1/N_z, 2/N_z, \ldots, \lfloor N_z/2 \rfloor/N_z`. The number of
   q-points is therefore :math:`\lfloor N_z/2 \rfloor + 1`, or
   :math:`N_z/2 + 1` when :math:`N_z` is even.

   A non-orthogonal real-space cell does not by itself change this counting if
   the supercell is still a diagonal repetition of the primitive lattice in
   reduced coordinates. What matters is the integer repetition matrix
   :math:`P` satisfying ``supercell = P @ primitive``.

   For a non-diagonal or transformed supercell, the number of q-points depends
   on the selected path direction. pySED keeps only fractional positions
   :math:`f` on the line
   :math:`q(f)=q_{\mathrm{start}}+f(q_{\mathrm{end}}-q_{\mathrm{start}})`
   that satisfy :math:`q(f)P^T \in \mathbb{Z}^3`. In this case, inspect the
   printed ``Number of q-points generated`` message or the output ``.Qpts``
   file. Increase the supercell repeats in the lattice directions that project
   onto the desired q-path to improve resolution.

rescale_prim
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   rescale_prim = 1

**Meaning**
   Reconstruct the primitive cell from the actual trajectory cell when the
   trajectory cell differs from the ideal ``prim_unitcell * supercell_dim``.

**Default**
   ``1``.

**Notes**
   This is useful after NPT relaxation. For strict NVT workflows, the trajectory
   cell should normally match the expected cell.

Q-Point Path
------------

num_qpaths
~~~~~~~~~~

**Syntax**

.. code-block:: text

   num_qpaths = 3

**Meaning**
   Number of path segments.

**Default**
   ``None``. Set it explicitly.

**Notes**
   A path with labels ``GMKG`` has three segments: Gamma to M, M to K, and K to
   Gamma.

q_path_name
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   q_path_name = 'GMKG'

**Meaning**
   Labels for high-symmetry points used on the SED plot.

**Default**
   ``GA``.

**Notes**
   Use ``G`` for Gamma. The number of labels must be ``num_qpaths + 1``.

q_path
~~~~~~

**Syntax**

.. code-block:: text

   q_path = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333333 0.3333333 0.0  0.0 0.0 0.0

**Meaning**
   High-symmetry q-points in reduced reciprocal coordinates.

**Default**
   No useful default. Set it explicitly.

**Notes**
   Provide ``num_qpaths + 1`` triples. Fractions such as ``1/3`` are accepted.
   pySED keeps only q-points commensurate with the supercell.

Plotting Options
----------------

plot_SED
~~~~~~~~

**Syntax**

.. code-block:: text

   plot_SED = 0

**Meaning**
   Select compute mode or plot/fitting mode.

**Default**
   ``0``.

**Modes**
   ``0`` computes SED data from the trajectory. ``1`` reads existing SED output
   files and plots or fits them.

plot_cutoff_freq
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_cutoff_freq = 20

**Meaning**
   Maximum frequency shown in the SED map, in THz.

**Default**
   ``None``.

**Notes**
   If omitted, pySED plots the available FFT frequency range.

plot_interval
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_interval = 5

**Meaning**
   Frequency tick interval on the y-axis, in THz.

**Default**
   ``5``.

plot_color
~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_color = 'RdBu_r'

**Meaning**
   Matplotlib colormap for the SED intensity.

**Default**
   ``RdBu_r``.

**Examples**
   ``jet``, ``Spectral``, ``inferno``.

colorbar_min
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   colorbar_min = -20

**Meaning**
   Lower colorbar limit after log scaling.

**Default**
   ``None``.

**Notes**
   Run an initial plot first, then tune this value for contrast.

colorbar_max
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   colorbar_max = -4

**Meaning**
   Upper colorbar limit after log scaling.

**Default**
   ``None``.

use_contourf
~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   use_contourf = 1

**Meaning**
   Force contour-style plotting instead of image-style plotting.

**Default**
   ``0``.

**Notes**
   pySED uses contour plotting automatically for multi-segment q-paths. For a
   single path, set this to ``1`` when the contour style gives a cleaner figure.

qpoint_slice_index
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   qpoint_slice_index = 2

**Meaning**
   Zero-based q-point index used for single-q-point slice plotting and fitting.

**Default**
   ``0``.

plot_slice
~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_slice = 1

**Meaning**
   Plot the SED spectrum at one q-point.

**Default**
   ``0``.

**Notes**
   Use this before all-q-point Lorentz fitting to tune peak detection settings.

plot_partial_SED
~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   plot_partial_SED = 3
   plot_partial_SED = 3 x

**Meaning**
   Plot the partial SED for a one-based atom type index. With only the type
   index, pySED sums x, y, and z directions. With a direction, pySED plots only
   that component.

**Default**
   ``0``.

**Allowed directions**
   ``x``, ``y``, or ``z``.

**Notes**
   Requires partial files generated earlier with ``output_partial = 1``.

if_show_figures
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   if_show_figures = 0

**Meaning**
   Show figures interactively.

**Default**
   ``0``.

**Notes**
   Use ``0`` on clusters or headless environments.

Lorentzian Fitting
------------------

lorentz
~~~~~~~

**Syntax**

.. code-block:: text

   lorentz = 1

**Meaning**
   Enable Lorentzian fitting for the selected q-point or for all q-points.

**Default**
   ``0``.

peak_height
~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_height = 8.0e-7

**Meaning**
   Minimum peak height passed to ``scipy.signal.find_peaks``.

**Default**
   No parser default.

**Notes**
   If too high, real peaks are missed. If too low, noise may be fitted.

peak_prominence
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_prominence = 6.0e-7

**Meaning**
   Minimum peak prominence passed to ``scipy.signal.find_peaks``.

**Default**
   No parser default.

**Notes**
   Prominence measures how strongly a peak stands above the surrounding
   baseline.

initial_guess_hwhm
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   initial_guess_hwhm = 0.001

**Meaning**
   Initial guess for the Lorentzian HWHM.

**Default**
   ``0.001``.

peak_max_hwhm
~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   peak_max_hwhm = 0.1

**Meaning**
   Upper bound for the fitted HWHM.

**Default**
   ``1e6``.

**Notes**
   Set a smaller value when broad or noisy peaks cause unreasonable fits.

lorentz_fit_cutoff
~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lorentz_fit_cutoff = 20

**Meaning**
   Maximum frequency included in Lorentzian fitting, in THz.

**Default**
   ``None``.

modulate_factor
~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   modulate_factor = 2

**Meaning**
   Shrinks the fitting range around each detected peak by shifting the left and
   right fitting boundaries inward.

**Default**
   ``0``.

**Notes**
   This can help when neighboring peaks or noisy baselines distort the fit.

lorentz_fit_all_qpoint
~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   lorentz_fit_all_qpoint = 1

**Meaning**
   Fit all q-points and collect frequency-lifetime data.

**Default**
   ``0``.

**Notes**
   Set ``lorentz = 1`` at the same time. Tune fitting parameters on selected
   q-points first.

re_output_total_freq_lifetime
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Syntax**

.. code-block:: text

   re_output_total_freq_lifetime = 1

**Meaning**
   Rebuild ``TOTAL-LORENTZ-Qpoints.Fre_lifetime`` after re-fitting selected
   q-points.

**Default**
   ``0``.

**Notes**
   Use this with ``lorentz_fit_all_qpoint = 0`` after improving a single
   q-point fit.

Reserved Parameter
------------------

with_eigs
~~~~~~~~~

**Syntax**

.. code-block:: text

   with_eigs = 1

**Meaning**
   Reserved for eigenvector-related development.

**Default**
   ``None``.

**Notes**
   This key is accepted by the parser but is not part of the current recommended
   user workflow. pySED currently uses the eigenvector-free SED expression.

Minimal ``input_SED.in`` Example
--------------------------------

The following compact example is adapted from the
`MoS2 GPUMD input_SED.in <https://github.com/Tingliangstu/pySED/blob/main/example/MoS2_gpumd/SED/input_SED.in>`_.
It is a good starting point for a GPUMD trajectory. First run it with
``plot_SED = 0`` to compute the SED data. After the ``.SED``, ``.Qpts``, and
``.THz`` files exist, change ``plot_SED = 1`` to plot or fit the result.

.. code-block:: text

   # MD trajectory information
   num_atoms = 13824
   total_num_steps = 500000
   time_step = 1
   output_data_stride = 50

   # Averaging and compression
   num_splits = 5
   compress = 1

   # Files
   out_files_name = 'bulk_MoS2'
   basis_lattice_file = '../structure/basis.in'
   dump_xyz_file = '../gpumd_run/dump.xyz'
   output_hdf5 = 'vel_pos_compress.hdf5'
   file_format = 'gpumd'

   # Parallelization
   use_parallel = 1
   max_cores = 4

   # Crystal structure
   prim_unitcell = 3.11905 0 0 -1.55952 2.70117 0 0 0 12.6227
   supercell_dim = 12 12 16
   rescale_prim = 1

   # Q-path: Gamma to A
   num_qpaths = 1
   q_path_name = 'GA'
   q_path = 0.0 0.0 0.0  0.0 0.0 0.5

   # Compute first. Change to plot_SED = 1 after SED files exist.
   plot_SED = 0
   plot_cutoff_freq = 2.0
   plot_interval = 0.5
   qpoint_slice_index = 0
   plot_slice = 1
   if_show_figures = 0

   # Fitting is off in the first run.
   lorentz = 0
   lorentz_fit_cutoff = 2
   peak_height = 1.0e-6
   peak_prominence = 1.0e-7
   initial_guess_hwhm = 0.0005
   lorentz_fit_all_qpoint = 0
