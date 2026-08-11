pySED Documentation
===================

**pySED** calculates kinetic-energy-weighted phonon spectral energy density
from molecular-dynamics trajectories. Version 3.0 reads GPUMD extended XYZ,
combined LAMMPS custom dumps, and NetCDF; supports CPU multiprocessing and an
optional CuPy backend; and provides plotting, partial SED, robust peak
detection, line-shape fitting, and linewidth-derived lifetimes.

.. toctree::
   :maxdepth: 1
   :caption: Table of contents

   introduction
   requirements
   installation
   starting
   command_line
   sed_workflow
   input_parameters
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
   tutorials
   example
   theory
   migration_v3
   tips
   publications
   reference
   troubleshooting
