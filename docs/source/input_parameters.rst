Input Parameters
================

pySED reads a plain-text control file. The preferred file name is
``input.in``, but another file name can be passed on the command line. Each
active line has the form

.. code-block:: text

   parameter_name = value

Blank lines and text after ``#`` are ignored. If no input file is passed on the
command line, pySED searches for ``input.in`` and then the legacy filename
``input_SED.in`` in the current directory:

.. code-block:: bash

   pysed
   pysed input.in
   pySED input.in

Boolean parameters accept ``True``/``False`` (case-insensitive) or ``1``/``0``.

The table below lists the parameters accepted by the pySED input file. Each
parameter has its own page with the parameter name as the title, plus syntax,
meaning, defaults, examples, practical notes, and related parameters.

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Input file parameters

   input_parameters/num_atoms
   input_parameters/total_num_steps
   input_parameters/time_step
   input_parameters/output_data_stride
   input_parameters/out_files_name
   input_parameters/legacy_unit
   input_parameters/basis_lattice_file
   input_parameters/lammps_unit
   input_parameters/output_partial
   input_parameters/max_cores
   input_parameters/prim_unitcell
   input_parameters/prim_axis
   input_parameters/supercell_dim
   input_parameters/rescale_prim
   input_parameters/num_qpaths
   input_parameters/q_path_name
   input_parameters/q_path
   input_parameters/plot_cutoff_freq
   input_parameters/plot_interval
   input_parameters/plot_color
   input_parameters/colorbar_min
   input_parameters/colorbar_max
   input_parameters/use_contourf
   input_parameters/qpoint_slice_index
   input_parameters/plot_slice
   input_parameters/plot_partial_SED
   input_parameters/if_show_figures
   input_parameters/initial_guess_hwhm
   input_parameters/peak_max_hwhm
   input_parameters/modulate_factor
   input_parameters/lorentz_fit_all_qpoint
   input_parameters/re_output_total_freq_lifetime

.. raw:: html

   <script>
   (function () {
     var pages = {
       "action": "input_parameters/action.html",
       "method": "input_parameters/method.html",
       "num-atoms": "input_parameters/num_atoms.html",
       "total-num-steps": "input_parameters/total_num_steps.html",
       "time-step": "input_parameters/time_step.html",
       "output-data-stride": "input_parameters/output_data_stride.html",
       "out-files-name": "input_parameters/out_files_name.html",
       "legacy-unit": "input_parameters/legacy_unit.html",
       "basis-lattice-file": "input_parameters/basis_lattice_file.html",
       "lammps-unit": "input_parameters/lammps_unit.html",
       "output-partial": "input_parameters/output_partial.html",
       "trajectory-file": "input_parameters/trajectory_file.html",
       "trajectory-read-mode": "input_parameters/trajectory_read_mode.html",
       "trajectory-prefetch": "input_parameters/trajectory_prefetch.html",
       "netcdf-batch-size": "input_parameters/netcdf_batch_size.html",
       "netcdf-compression-level": "input_parameters/netcdf_compression_level.html",
       "max-cores": "input_parameters/max_cores.html",
       "backend": "input_parameters/backend.html",
       "prim-unitcell": "input_parameters/prim_unitcell.html",
       "prim-axis": "input_parameters/prim_axis.html",
       "supercell-dim": "input_parameters/supercell_dim.html",
       "rescale-prim": "input_parameters/rescale_prim.html",
       "dump-group-supercell": "input_parameters/dump_group_supercell.html",
       "num-qpaths": "input_parameters/num_qpaths.html",
       "q-path-name": "input_parameters/q_path_name.html",
       "q-path": "input_parameters/q_path.html",
       "plot-cutoff-freq": "input_parameters/plot_cutoff_freq.html",
       "plot-interval": "input_parameters/plot_interval.html",
       "plot-color": "input_parameters/plot_color.html",
       "colorbar-min": "input_parameters/colorbar_min.html",
       "colorbar-max": "input_parameters/colorbar_max.html",
       "use-contourf": "input_parameters/use_contourf.html",
       "qpoint-slice-index": "input_parameters/qpoint_slice_index.html",
       "plot-slice": "input_parameters/plot_slice.html",
       "plot-partial-sed": "input_parameters/plot_partial_SED.html",
       "if-show-figures": "input_parameters/if_show_figures.html",
       "initial-guess-hwhm": "input_parameters/initial_guess_hwhm.html",
       "peak-max-hwhm": "input_parameters/peak_max_hwhm.html",
       "modulate-factor": "input_parameters/modulate_factor.html",
       "lorentz-fit-all-qpoint": "input_parameters/lorentz_fit_all_qpoint.html",
       "re-output-total-freq-lifetime": "input_parameters/re_output_total_freq_lifetime.html",
       "lorentz-fit-freq-min": "input_parameters/lorentz_fit_freq_min.html",
       "lorentz-fit-freq-max": "input_parameters/lorentz_fit_freq_max.html",
       "fitting-function": "input_parameters/fitting_function.html",
       "peak-min-significance": "input_parameters/peak_min_significance.html",
       "num-blocks": "input_parameters/num_blocks.html",
       "num-splits": "input_parameters/num_blocks.html",
       "compress": "input_parameters/trajectory_read_mode.html",
       "output-hdf5": "input_parameters/trajectory_read_mode.html",
       "file-format": "input_parameters/trajectory_file.html",
       "dump-xyz-file": "input_parameters/trajectory_file.html",
       "pos-file": "input_parameters/trajectory_file.html",
       "vels-file": "input_parameters/trajectory_file.html",
       "use-parallel": "input_parameters/backend.html",
       "plot-sed": "input_parameters/action.html",
       "lorentz": "input_parameters/fitting_function.html",
       "peak-height": "input_parameters/peak_min_significance.html",
       "peak-prominence": "input_parameters/peak_min_significance.html",
       "lorentz-fit-cutoff": "input_parameters/lorentz_fit_freq_max.html",
       "with-eigs": "migration_v3.html"
     };
     function redirectHash() {
       var key = window.location.hash ? window.location.hash.slice(1) : "";
       if (pages[key]) {
         window.location.replace(pages[key]);
       }
     }
     redirectHash();
     window.addEventListener("hashchange", redirectHash);
   }());
   </script>

.. rubric:: Parameter pages

.. rst-class:: parameter-index

.. list-table::
   :header-rows: 1
   :widths: 24 24 52

   * - Parameter
     - Group
     - Purpose
   * - `action <input_parameters/action.html>`_
     - Workflow
     - Select automatic dispatch, compute, plot, or fit.
   * - `method <input_parameters/method.html>`_
     - Workflow
     - Select the SED workflow; ``sed`` is the only accepted value.
   * - `num_atoms <input_parameters/num_atoms.html>`_
     - MD simulation
     - Number of atoms in the MD trajectory.
   * - `total_num_steps <input_parameters/total_num_steps.html>`_
     - MD simulation
     - Number of MD steps used for the production trajectory.
   * - `time_step <input_parameters/time_step.html>`_
     - MD simulation
     - MD time step in fs.
   * - `output_data_stride <input_parameters/output_data_stride.html>`_
     - MD simulation
     - Trajectory output interval; controls the maximum resolvable frequency.
   * - `out_files_name <input_parameters/out_files_name.html>`_
     - Files
     - Prefix for pySED output files.
   * - `legacy_unit <input_parameters/legacy_unit.html>`_
     - Output scale
     - Temporarily use the pySED 2.x intensity scale for historical comparisons.
   * - `basis_lattice_file <input_parameters/basis_lattice_file.html>`_
     - Files
     - Path to ``basis.in``.
   * - `lammps_unit <input_parameters/lammps_unit.html>`_
     - Files
     - LAMMPS unit style for velocity conversion.
   * - `output_partial <input_parameters/output_partial.html>`_
     - Files
     - Write element and direction-resolved SED files.
   * - `trajectory_file <input_parameters/trajectory_file.html>`_
     - Files
     - Name a GPUMD extxyz, combined LAMMPS dump, or NetCDF trajectory.
   * - `trajectory_read_mode <input_parameters/trajectory_read_mode.html>`_
     - Files
     - Choose converted NetCDF caching or direct text reading.
   * - `trajectory_prefetch <input_parameters/trajectory_prefetch.html>`_
     - Files
     - Prepare one trajectory block ahead in a background thread.
   * - `netcdf_batch_size <input_parameters/netcdf_batch_size.html>`_
     - Files
     - Set frames parsed per text-to-NetCDF batch.
   * - `netcdf_compression_level <input_parameters/netcdf_compression_level.html>`_
     - Files
     - Set zlib compression for generated NetCDF caches.
   * - `max_cores <input_parameters/max_cores.html>`_
     - Parallelization
     - Maximum number of worker processes.
   * - `backend <input_parameters/backend.html>`_
     - Parallelization
     - Select NumPy or the optional CuPy backend.
   * - `prim_unitcell <input_parameters/prim_unitcell.html>`_
     - Crystal structure
     - Primitive lattice vectors.
   * - `prim_axis <input_parameters/prim_axis.html>`_
     - Crystal structure
     - Primitive-axis transformation.
   * - `supercell_dim <input_parameters/supercell_dim.html>`_
     - Crystal structure
     - Supercell repeats; controls q-point resolution.
   * - `rescale_prim <input_parameters/rescale_prim.html>`_
     - Crystal structure
     - Reconstruct the primitive cell from the trajectory cell.
   * - `dump_group_supercell <input_parameters/dump_group_supercell.html>`_
     - Crystal structure
     - Set the periodic cell represented by a locally dumped group.
   * - `num_qpaths <input_parameters/num_qpaths.html>`_
     - Q-path
     - Number of high-symmetry path segments.
   * - `q_path_name <input_parameters/q_path_name.html>`_
     - Q-path
     - Labels for high-symmetry points.
   * - `q_path <input_parameters/q_path.html>`_
     - Q-path
     - High-symmetry q-points in reduced coordinates.
   * - `plot_cutoff_freq <input_parameters/plot_cutoff_freq.html>`_
     - Plotting
     - Maximum plotted frequency in THz.
   * - `plot_interval <input_parameters/plot_interval.html>`_
     - Plotting
     - Frequency tick interval.
   * - `plot_color <input_parameters/plot_color.html>`_
     - Plotting
     - Matplotlib colormap.
   * - `colorbar_min <input_parameters/colorbar_min.html>`_
     - Plotting
     - Lower log-scale colorbar limit.
   * - `colorbar_max <input_parameters/colorbar_max.html>`_
     - Plotting
     - Upper log-scale colorbar limit.
   * - `use_contourf <input_parameters/use_contourf.html>`_
     - Plotting
     - Use contour-style SED plotting.
   * - `qpoint_slice_index <input_parameters/qpoint_slice_index.html>`_
     - Plotting
     - Zero-based q-point index for slice plots and fitting.
   * - `plot_slice <input_parameters/plot_slice.html>`_
     - Plotting
     - Plot one q-point spectrum.
   * - `plot_partial_SED <input_parameters/plot_partial_SED.html>`_
     - Plotting
     - Plot one element or direction-resolved partial SED.
   * - `if_show_figures <input_parameters/if_show_figures.html>`_
     - Plotting
     - Show figures interactively.
   * - `initial_guess_hwhm <input_parameters/initial_guess_hwhm.html>`_
     - Lorentzian fitting
     - Initial HWHM guess.
   * - `peak_max_hwhm <input_parameters/peak_max_hwhm.html>`_
     - Lorentzian fitting
     - Maximum fitted HWHM.
   * - `modulate_factor <input_parameters/modulate_factor.html>`_
     - Lorentzian fitting
     - Shrink fitting windows around detected peaks.
   * - `lorentz_fit_all_qpoint <input_parameters/lorentz_fit_all_qpoint.html>`_
     - Lorentzian fitting
     - Fit peaks at all q-points.
   * - `re_output_total_freq_lifetime <input_parameters/re_output_total_freq_lifetime.html>`_
     - Lorentzian fitting
     - Rebuild total frequency-lifetime output after refitting.
   * - `lorentz_fit_freq_min <input_parameters/lorentz_fit_freq_min.html>`_
     - Peak fitting
     - Set minimum fitted frequency in THz.
   * - `lorentz_fit_freq_max <input_parameters/lorentz_fit_freq_max.html>`_
     - Peak fitting
     - Set maximum fitted frequency in THz.
   * - `fitting_function <input_parameters/fitting_function.html>`_
     - Peak fitting
     - Select Lorentz, velocity-DHO, or automatic AICc choice.
   * - `peak_min_significance <input_parameters/peak_min_significance.html>`_
     - Peak fitting
     - Set the robust local peak-significance threshold.
   * - `num_blocks <input_parameters/num_blocks.html>`_
     - Averaging
     - Set equal trajectory blocks used for spectral averaging.
