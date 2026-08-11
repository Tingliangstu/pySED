Input Parameters
================

pySED reads a plain-text control file. The preferred filename is ``input.in``;
another path can be passed on the command line. Each active line has the form

.. code-block:: text

   parameter_name = value

Blank lines and text following ``#`` are ignored. With no argument, pySED
searches for ``input.in`` first and then the legacy filename
``input_SED.in``:

.. code-block:: bash

   pySED input.in
   pySED

The table contains all 45 public parameters in the pySED 3.0 SED workflow.
Every parameter name opens a dedicated page with syntax, defaults, constraints,
examples, and related settings.

.. toctree::
   :hidden:
   :maxdepth: 1

   input_parameters/lammps_unit
   input_parameters/num_atoms
   input_parameters/total_num_steps
   input_parameters/time_step
   input_parameters/output_data_stride
   input_parameters/max_cores
   input_parameters/out_files_name
   input_parameters/if_show_figures
   input_parameters/basis_lattice_file
   input_parameters/prim_unitcell
   input_parameters/prim_axis
   input_parameters/supercell_dim
   input_parameters/rescale_prim
   input_parameters/num_qpaths
   input_parameters/q_path_name
   input_parameters/q_path
   input_parameters/output_partial
   input_parameters/plot_partial_SED
   input_parameters/plot_cutoff_freq
   input_parameters/plot_interval
   input_parameters/plot_color
   input_parameters/colorbar_min
   input_parameters/colorbar_max
   input_parameters/use_contourf
   input_parameters/plot_slice
   input_parameters/qpoint_slice_index
   input_parameters/lorentz_fit_all_qpoint
   input_parameters/initial_guess_hwhm
   input_parameters/peak_max_hwhm
   input_parameters/modulate_factor
   input_parameters/re_output_total_freq_lifetime

Complete parameter index
------------------------

.. list-table::
   :class: parameter-index
   :header-rows: 1
   :widths: 28 17 55

   * - Parameter
     - Default
     - Purpose
   * - :doc:`action <input_parameters/action>`
     - ``thinking``
     - Select automatic dispatch, compute, plot, or fit.
   * - :doc:`method <input_parameters/method>`
     - ``sed``
     - Select the SED workflow; ``sed`` is the only accepted value.
   * - :doc:`backend <input_parameters/backend>`
     - ``numpy``
     - Select NumPy or the optional CuPy backend.
   * - :doc:`trajectory_file <input_parameters/trajectory_file>`
     - ``dump.xyz``
     - Name a GPUMD extxyz, combined LAMMPS dump, or NetCDF trajectory.
   * - :doc:`trajectory_read_mode <input_parameters/trajectory_read_mode>`
     - ``cache``
     - Choose converted NetCDF caching or direct text reading.
   * - :doc:`trajectory_prefetch <input_parameters/trajectory_prefetch>`
     - ``1``
     - Prepare one trajectory block ahead in a background thread.
   * - :doc:`netcdf_batch_size <input_parameters/netcdf_batch_size>`
     - ``64``
     - Set frames parsed per text-to-NetCDF batch.
   * - :doc:`netcdf_compression_level <input_parameters/netcdf_compression_level>`
     - ``1``
     - Set zlib compression for generated NetCDF caches.
   * - :doc:`lammps_unit <input_parameters/lammps_unit>`
     - ``metal``
     - Select LAMMPS velocity-unit conversion.
   * - :doc:`num_atoms <input_parameters/num_atoms>`
     - ``0``
     - Set the number of trajectory atoms.
   * - :doc:`total_num_steps <input_parameters/total_num_steps>`
     - ``0``
     - Set MD production steps represented by the requested frames.
   * - :doc:`time_step <input_parameters/time_step>`
     - ``0.0``
     - Set the MD integration step in fs.
   * - :doc:`output_data_stride <input_parameters/output_data_stride>`
     - ``0``
     - Set MD steps between saved trajectory frames.
   * - :doc:`num_blocks <input_parameters/num_blocks>`
     - ``5``
     - Set equal trajectory blocks used for spectral averaging.
   * - :doc:`max_cores <input_parameters/max_cores>`
     - ``4``
     - Limit NumPy worker processes.
   * - :doc:`out_files_name <input_parameters/out_files_name>`
     - ``pySED``
     - Set the output prefix.
   * - :doc:`if_show_figures <input_parameters/if_show_figures>`
     - ``0``
     - Display figures in addition to saving them.
   * - :doc:`basis_lattice_file <input_parameters/basis_lattice_file>`
     - ``basis.in``
     - Map trajectory atoms to repeated cells, basis atoms, and masses.
   * - :doc:`prim_unitcell <input_parameters/prim_unitcell>`
     - ``None``
     - Set primitive-cell lattice vectors.
   * - :doc:`prim_axis <input_parameters/prim_axis>`
     - ``None``
     - Apply an optional primitive-axis transformation.
   * - :doc:`supercell_dim <input_parameters/supercell_dim>`
     - ``1 1 1``
     - Set supercell repetitions and commensurate Q resolution.
   * - :doc:`dump_group_supercell <input_parameters/dump_group_supercell>`
     - ``None``
     - Set the periodic cell represented by a locally dumped group.
   * - :doc:`rescale_prim <input_parameters/rescale_prim>`
     - ``1``
     - Reconstruct the primitive cell from the trajectory cell.
   * - :doc:`num_qpaths <input_parameters/num_qpaths>`
     - ``1``
     - Set the number of connected Q-path segments.
   * - :doc:`q_path_name <input_parameters/q_path_name>`
     - ``GA``
     - Label the Q-path vertices.
   * - :doc:`q_path <input_parameters/q_path>`
     - ``None``
     - Define Q-path vertices in reduced coordinates.
   * - :doc:`output_partial <input_parameters/output_partial>`
     - ``0``
     - Save element- and Cartesian-resolved SED components.
   * - :doc:`plot_partial_SED <input_parameters/plot_partial_SED>`
     - disabled
     - Select an element and optional Cartesian component to plot.
   * - :doc:`plot_cutoff_freq <input_parameters/plot_cutoff_freq>`
     - ``None``
     - Set maximum plotted frequency in THz.
   * - :doc:`plot_interval <input_parameters/plot_interval>`
     - ``5.0``
     - Set frequency-axis tick spacing.
   * - :doc:`plot_color <input_parameters/plot_color>`
     - ``RdBu_r``
     - Select the dispersion colormap.
   * - :doc:`colorbar_min <input_parameters/colorbar_min>`
     - ``None``
     - Set an optional lower logarithmic color limit.
   * - :doc:`colorbar_max <input_parameters/colorbar_max>`
     - ``None``
     - Set an optional upper logarithmic color limit.
   * - :doc:`use_contourf <input_parameters/use_contourf>`
     - ``0``
     - Draw the SED heatmap with filled contours.
   * - :doc:`plot_slice <input_parameters/plot_slice>`
     - ``0``
     - Plot one Q-point spectrum.
   * - :doc:`qpoint_slice_index <input_parameters/qpoint_slice_index>`
     - ``0``
     - Select the zero-based Q point for plotting or fitting.
   * - :doc:`lorentz_fit_all_qpoint <input_parameters/lorentz_fit_all_qpoint>`
     - ``0``
     - Fit every Q point or only the selected Q point.
   * - :doc:`lorentz_fit_freq_min <input_parameters/lorentz_fit_freq_min>`
     - ``None``
     - Set minimum fitted frequency in THz.
   * - :doc:`lorentz_fit_freq_max <input_parameters/lorentz_fit_freq_max>`
     - ``None``
     - Set maximum fitted frequency in THz.
   * - :doc:`fitting_function <input_parameters/fitting_function>`
     - ``auto``
     - Select Lorentz, velocity-DHO, or automatic AICc choice.
   * - :doc:`peak_min_significance <input_parameters/peak_min_significance>`
     - ``4.0``
     - Set the robust local peak-significance threshold.
   * - :doc:`initial_guess_hwhm <input_parameters/initial_guess_hwhm>`
     - ``0.001``
     - Set optimizer initial HWHM in THz.
   * - :doc:`peak_max_hwhm <input_parameters/peak_max_hwhm>`
     - ``1e6``
     - Set fitted HWHM upper bound in THz.
   * - :doc:`modulate_factor <input_parameters/modulate_factor>`
     - ``0``
     - Trim samples from both ends of each local fit range.
   * - :doc:`re_output_total_freq_lifetime <input_parameters/re_output_total_freq_lifetime>`
     - ``0``
     - Rebuild combined lifetime output after a single-Q refit.
