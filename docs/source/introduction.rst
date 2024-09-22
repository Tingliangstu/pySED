.. highlight:: rst

Introduction
============

The **Spectral Energy Density (SED)** technique employed by **pySED** is a method that can directly 
predict phonon dispersion relations and their corresponding lifetimes from atomic velocities 
obtained in large-scale molecular dynamics (MD) simulations. 
SED represents the intensity of the kinetic energy distribution of a system across different phonon modes, 
effectively providing a dispersion weighted by the kinetic energy carried by each phonon.

**pySED** can read and process large-scale simulation trajectories from both `LAMMPS <https://www.lammps.org/#gsc.tab=0>`_ and `GPUMD <https://gpumd.org/>`_,
enhancing its versatility and applicability in various research scenarios.

Main Features
-------------
- **SED Visualization**: Generates high-quality SED images for detailed analysis of phonon dispersion and lifetimes.
- **Lorentzian Fitting**: Utilizes Lorentzian fitting to accurately model SED peaks.
- **Comprehensive Peak Fitting**: Capable of fitting all peaks simultaneously to output phonon lifetimes or fitting them individually as needed.
- **User-Friendly Interface**: Designed for simplicity and ease of use.
- **Qualitative Phonon Lifetimes**: Provides qualitative insights into phonon lifetimes based on SED analysis.
- **Parallel Processing**: Supports parallel execution to accelerate SED computations and obtain results quickly.
