.. pySED documentation master file, created by
   sphinx-quickstart on Sun Sep 22 01:30:31 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pySED's documentation!
=================================

**pySED** is a Python package for calculating kinetic-energy-weighted phonon
spectral energy density (SED) from molecular dynamics trajectories. It can read
trajectories from GPUMD and LAMMPS, construct commensurate q-points, plot SED
maps, and fit Lorentzian peaks to estimate phonon lifetimes.

.. toctree::
   :maxdepth: 1
   :caption: Table of contents

   introduction
   requirements
   installation
   starting
   input_parameters
   tutorials
   example
   theory
   tips
   publications
   reference
   troubleshooting

.. toctree::
   :maxdepth: 1
   :caption: pySED 3.0 additions

   command_line
   sed_workflow
   input_parameters/action
   input_parameters/method
   input_parameters/backend
   input_parameters/trajectory_file
   input_parameters/trajectory_read_mode
   input_parameters/trajectory_prefetch
   input_parameters/netcdf_batch_size
   input_parameters/netcdf_compression_level
   input_parameters/num_blocks
   input_parameters/dump_group_supercell
   input_parameters/fitting_function
   input_parameters/peak_min_significance
   input_parameters/lorentz_fit_freq_min
   input_parameters/lorentz_fit_freq_max
   peak_fitting
   output_files
   sed_units
   migration_v3
   code_architecture
