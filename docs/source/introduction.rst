Introduction
============

The spectral energy density (SED) method projects atomic velocities from
molecular-dynamics trajectories into reciprocal space. The resulting
kinetic-energy-weighted map reveals finite-temperature vibrational dispersion;
line widths can be fitted to estimate relaxation times.

pySED 3.0 provides:

- GPUMD extended-XYZ, combined LAMMPS custom-dump, and NetCDF readers;
- reusable NetCDF caches or direct sequential text reading;
- block-averaged NumPy computation with persistent worker processes;
- an optional one-GPU CuPy backend;
- exact rational Q-path input and commensurate finite-supercell Q points;
- total, element-resolved, and Cartesian-resolved SED;
- eV/THz one-sided spectra with energy-preserving normalization;
- dispersion and single-Q plotting;
- robust local peak detection, Lorentz and velocity-DHO fitting, automatic
  AICc model choice, and linewidth-derived lifetime output.

The production workflow is eigenvector-free. Fitted peaks are therefore
spectral features rather than automatic phonon-branch assignments; difficult
or overlapping features require physical inspection.
