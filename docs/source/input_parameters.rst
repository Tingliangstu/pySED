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

The table below is the left-sidebar index for ``input_SED.in`` input
parameters. Click a parameter name to open a dedicated page with its syntax,
meaning, defaults, examples, practical notes, and related parameters.

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: input_SED.in parameters

   input_parameters/num_atoms
   input_parameters/total_num_steps
   input_parameters/time_step
   input_parameters/output_data_stride
   input_parameters/num_splits
   input_parameters/compress
   input_parameters/out_files_name
   input_parameters/basis_lattice_file
   input_parameters/file_format
   input_parameters/dump_xyz_file
   input_parameters/pos_file
   input_parameters/vels_file
   input_parameters/lammps_unit
   input_parameters/output_hdf5
   input_parameters/output_partial
   input_parameters/use_parallel
   input_parameters/max_cores
   input_parameters/prim_unitcell
   input_parameters/prim_axis
   input_parameters/supercell_dim
   input_parameters/rescale_prim
   input_parameters/num_qpaths
   input_parameters/q_path_name
   input_parameters/q_path
   input_parameters/plot_SED
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
   input_parameters/lorentz
   input_parameters/peak_height
   input_parameters/peak_prominence
   input_parameters/initial_guess_hwhm
   input_parameters/peak_max_hwhm
   input_parameters/lorentz_fit_cutoff
   input_parameters/modulate_factor
   input_parameters/lorentz_fit_all_qpoint
   input_parameters/re_output_total_freq_lifetime
   input_parameters/with_eigs


Parameter Index
---------------

.. raw:: html

   <span id="num-atoms"></span>
   <span id="total-num-steps"></span>
   <span id="time-step"></span>
   <span id="output-data-stride"></span>
   <span id="num-splits"></span>
   <span id="compress"></span>
   <span id="out-files-name"></span>
   <span id="basis-lattice-file"></span>
   <span id="file-format"></span>
   <span id="dump-xyz-file"></span>
   <span id="pos-file"></span>
   <span id="vels-file"></span>
   <span id="lammps-unit"></span>
   <span id="output-hdf5"></span>
   <span id="output-partial"></span>
   <span id="use-parallel"></span>
   <span id="max-cores"></span>
   <span id="prim-unitcell"></span>
   <span id="prim-axis"></span>
   <span id="supercell-dim"></span>
   <span id="rescale-prim"></span>
   <span id="num-qpaths"></span>
   <span id="q-path-name"></span>
   <span id="q-path"></span>
   <span id="plot-sed"></span>
   <span id="plot-cutoff-freq"></span>
   <span id="plot-interval"></span>
   <span id="plot-color"></span>
   <span id="colorbar-min"></span>
   <span id="colorbar-max"></span>
   <span id="use-contourf"></span>
   <span id="qpoint-slice-index"></span>
   <span id="plot-slice"></span>
   <span id="plot-partial-sed"></span>
   <span id="if-show-figures"></span>
   <span id="lorentz"></span>
   <span id="peak-height"></span>
   <span id="peak-prominence"></span>
   <span id="initial-guess-hwhm"></span>
   <span id="peak-max-hwhm"></span>
   <span id="lorentz-fit-cutoff"></span>
   <span id="modulate-factor"></span>
   <span id="lorentz-fit-all-qpoint"></span>
   <span id="re-output-total-freq-lifetime"></span>
   <span id="with-eigs"></span>


.. rst-class:: parameter-index

.. list-table:: Parameter index for ``input_SED.in``
   :header-rows: 1
   :widths: 24 24 52

   * - Parameter
     - Group
     - Purpose
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
   * - `num_splits <input_parameters/num_splits.html>`_
     - Averaging
     - Number of trajectory blocks used for averaging.
   * - `compress <input_parameters/compress.html>`_
     - Averaging
     - Compress trajectory data to HDF5 before SED.
   * - `out_files_name <input_parameters/out_files_name.html>`_
     - Files
     - Prefix for pySED output files.
   * - `basis_lattice_file <input_parameters/basis_lattice_file.html>`_
     - Files
     - Path to ``basis.in``.
   * - `file_format <input_parameters/file_format.html>`_
     - Files
     - Select GPUMD or LAMMPS trajectory reader.
   * - `dump_xyz_file <input_parameters/dump_xyz_file.html>`_
     - Files
     - GPUMD extended XYZ trajectory path.
   * - `pos_file <input_parameters/pos_file.html>`_
     - Files
     - LAMMPS position trajectory path.
   * - `vels_file <input_parameters/vels_file.html>`_
     - Files
     - LAMMPS velocity trajectory path.
   * - `lammps_unit <input_parameters/lammps_unit.html>`_
     - Files
     - LAMMPS unit style for velocity conversion.
   * - `output_hdf5 <input_parameters/output_hdf5.html>`_
     - Files
     - Compressed trajectory database name.
   * - `output_partial <input_parameters/output_partial.html>`_
     - Files
     - Write atom-type and direction-resolved SED files.
   * - `use_parallel <input_parameters/use_parallel.html>`_
     - Parallelization
     - Enable multiprocessing over q-points.
   * - `max_cores <input_parameters/max_cores.html>`_
     - Parallelization
     - Maximum number of worker processes.
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
   * - `num_qpaths <input_parameters/num_qpaths.html>`_
     - Q-path
     - Number of high-symmetry path segments.
   * - `q_path_name <input_parameters/q_path_name.html>`_
     - Q-path
     - Labels for high-symmetry points.
   * - `q_path <input_parameters/q_path.html>`_
     - Q-path
     - High-symmetry q-points in reduced coordinates.
   * - `plot_SED <input_parameters/plot_SED.html>`_
     - Plotting
     - Select compute mode or plot/fitting mode.
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
     - Plot one atom-type or direction-resolved partial SED.
   * - `if_show_figures <input_parameters/if_show_figures.html>`_
     - Plotting
     - Show figures interactively.
   * - `lorentz <input_parameters/lorentz.html>`_
     - Lorentzian fitting
     - Enable Lorentzian peak fitting.
   * - `peak_height <input_parameters/peak_height.html>`_
     - Lorentzian fitting
     - Minimum detected peak height.
   * - `peak_prominence <input_parameters/peak_prominence.html>`_
     - Lorentzian fitting
     - Minimum detected peak prominence.
   * - `initial_guess_hwhm <input_parameters/initial_guess_hwhm.html>`_
     - Lorentzian fitting
     - Initial HWHM guess.
   * - `peak_max_hwhm <input_parameters/peak_max_hwhm.html>`_
     - Lorentzian fitting
     - Maximum fitted HWHM.
   * - `lorentz_fit_cutoff <input_parameters/lorentz_fit_cutoff.html>`_
     - Lorentzian fitting
     - Maximum fitting frequency in THz.
   * - `modulate_factor <input_parameters/modulate_factor.html>`_
     - Lorentzian fitting
     - Shrink fitting windows around detected peaks.
   * - `lorentz_fit_all_qpoint <input_parameters/lorentz_fit_all_qpoint.html>`_
     - Lorentzian fitting
     - Fit peaks at all q-points.
   * - `re_output_total_freq_lifetime <input_parameters/re_output_total_freq_lifetime.html>`_
     - Lorentzian fitting
     - Rebuild total frequency-lifetime output after refitting.
   * - `with_eigs <input_parameters/with_eigs.html>`_
     - Reserved
     - Reserved for eigenvector-related development.
