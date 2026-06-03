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

Use the table below as the quick index for ``input_SED.in``. Click a parameter
name to open its detailed syntax, meaning, defaults, examples, and practical
notes.

Parameter Index
---------------

.. rst-class:: parameter-index

.. list-table:: ``input_SED.in`` parameters
   :header-rows: 1
   :widths: 24 24 52

   * - Parameter
     - Group
     - Purpose
   * - `num_atoms <input_parameter_details.html#num-atoms>`_
     - MD simulation
     - Number of atoms in the MD trajectory.
   * - `total_num_steps <input_parameter_details.html#total-num-steps>`_
     - MD simulation
     - Number of MD steps used for the production trajectory.
   * - `time_step <input_parameter_details.html#time-step>`_
     - MD simulation
     - MD time step in fs.
   * - `output_data_stride <input_parameter_details.html#output-data-stride>`_
     - MD simulation
     - Trajectory output interval; controls the maximum resolvable frequency.
   * - `num_splits <input_parameter_details.html#num-splits>`_
     - Averaging
     - Number of trajectory blocks used for averaging.
   * - `compress <input_parameter_details.html#compress>`_
     - Averaging
     - Compress trajectory data to HDF5 before SED.
   * - `out_files_name <input_parameter_details.html#out-files-name>`_
     - Files
     - Prefix for pySED output files.
   * - `basis_lattice_file <input_parameter_details.html#basis-lattice-file>`_
     - Files
     - Path to ``basis.in``.
   * - `file_format <input_parameter_details.html#file-format>`_
     - Files
     - Select GPUMD or LAMMPS trajectory reader.
   * - `dump_xyz_file <input_parameter_details.html#dump-xyz-file>`_
     - Files
     - GPUMD extended XYZ trajectory path.
   * - `pos_file <input_parameter_details.html#pos-file>`_
     - Files
     - LAMMPS position trajectory path.
   * - `vels_file <input_parameter_details.html#vels-file>`_
     - Files
     - LAMMPS velocity trajectory path.
   * - `lammps_unit <input_parameter_details.html#lammps-unit>`_
     - Files
     - LAMMPS unit style for velocity conversion.
   * - `output_hdf5 <input_parameter_details.html#output-hdf5>`_
     - Files
     - Compressed trajectory database name.
   * - `output_partial <input_parameter_details.html#output-partial>`_
     - Files
     - Write atom-type and direction-resolved SED files.
   * - `use_parallel <input_parameter_details.html#use-parallel>`_
     - Parallelization
     - Enable multiprocessing over q-points.
   * - `max_cores <input_parameter_details.html#max-cores>`_
     - Parallelization
     - Maximum number of worker processes.
   * - `prim_unitcell <input_parameter_details.html#prim-unitcell>`_
     - Crystal structure
     - Primitive lattice vectors.
   * - `prim_axis <input_parameter_details.html#prim-axis>`_
     - Crystal structure
     - Primitive-axis transformation.
   * - `supercell_dim <input_parameter_details.html#supercell-dim>`_
     - Crystal structure
     - Supercell repeats; controls q-point resolution.
   * - `rescale_prim <input_parameter_details.html#rescale-prim>`_
     - Crystal structure
     - Reconstruct the primitive cell from the trajectory cell.
   * - `num_qpaths <input_parameter_details.html#num-qpaths>`_
     - Q-path
     - Number of high-symmetry path segments.
   * - `q_path_name <input_parameter_details.html#q-path-name>`_
     - Q-path
     - Labels for high-symmetry points.
   * - `q_path <input_parameter_details.html#q-path>`_
     - Q-path
     - High-symmetry q-points in reduced coordinates.
   * - `plot_SED <input_parameter_details.html#plot-sed>`_
     - Plotting
     - Select compute mode or plot/fitting mode.
   * - `plot_cutoff_freq <input_parameter_details.html#plot-cutoff-freq>`_
     - Plotting
     - Maximum plotted frequency in THz.
   * - `plot_interval <input_parameter_details.html#plot-interval>`_
     - Plotting
     - Frequency tick interval.
   * - `plot_color <input_parameter_details.html#plot-color>`_
     - Plotting
     - Matplotlib colormap.
   * - `colorbar_min <input_parameter_details.html#colorbar-min>`_
     - Plotting
     - Lower log-scale colorbar limit.
   * - `colorbar_max <input_parameter_details.html#colorbar-max>`_
     - Plotting
     - Upper log-scale colorbar limit.
   * - `use_contourf <input_parameter_details.html#use-contourf>`_
     - Plotting
     - Use contour-style SED plotting.
   * - `qpoint_slice_index <input_parameter_details.html#qpoint-slice-index>`_
     - Plotting
     - Zero-based q-point index for slice plots and fitting.
   * - `plot_slice <input_parameter_details.html#plot-slice>`_
     - Plotting
     - Plot one q-point spectrum.
   * - `plot_partial_SED <input_parameter_details.html#plot-partial-sed>`_
     - Plotting
     - Plot one atom-type or direction-resolved partial SED.
   * - `if_show_figures <input_parameter_details.html#if-show-figures>`_
     - Plotting
     - Show figures interactively.
   * - `lorentz <input_parameter_details.html#lorentz>`_
     - Lorentzian fitting
     - Enable Lorentzian peak fitting.
   * - `peak_height <input_parameter_details.html#peak-height>`_
     - Lorentzian fitting
     - Minimum detected peak height.
   * - `peak_prominence <input_parameter_details.html#peak-prominence>`_
     - Lorentzian fitting
     - Minimum detected peak prominence.
   * - `initial_guess_hwhm <input_parameter_details.html#initial-guess-hwhm>`_
     - Lorentzian fitting
     - Initial HWHM guess.
   * - `peak_max_hwhm <input_parameter_details.html#peak-max-hwhm>`_
     - Lorentzian fitting
     - Maximum fitted HWHM.
   * - `lorentz_fit_cutoff <input_parameter_details.html#lorentz-fit-cutoff>`_
     - Lorentzian fitting
     - Maximum fitting frequency in THz.
   * - `modulate_factor <input_parameter_details.html#modulate-factor>`_
     - Lorentzian fitting
     - Shrink fitting windows around detected peaks.
   * - `lorentz_fit_all_qpoint <input_parameter_details.html#lorentz-fit-all-qpoint>`_
     - Lorentzian fitting
     - Fit peaks at all q-points.
   * - `re_output_total_freq_lifetime <input_parameter_details.html#re-output-total-freq-lifetime>`_
     - Lorentzian fitting
     - Rebuild total frequency-lifetime output after refitting.
   * - `with_eigs <input_parameter_details.html#with-eigs>`_
     - Reserved
     - Reserved for eigenvector-related development.
