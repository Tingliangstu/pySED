Requirements
============

**pySED** requires a few python modules to be installed in your computer in order to work properly.
In general, pySED installation will automatically download the required packages.
All this packages can be downloaded easily using PyPI or conda repositories.

Mandatory requirements
----------------------
- Python >= 3.10
- numpy
- seaborn (Visual SED images)
- netCDF4 (Trajectory reading and reusable NetCDF caches)
- scipy (For fast Lorentz fitting)
- matplotlib (Plotting)

CuPy is optional for GPU computation; see :doc:`input_parameters/device`.
The MoS2 tutorial notebook also uses Jupyter and pandas.
