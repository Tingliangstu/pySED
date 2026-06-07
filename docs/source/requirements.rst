Requirements
============

**pySED** requires a few python modules to be installed in your computer in order to work properly.
In general, pySED installation will automatically download the required packages.
All this packages can be downloaded easily using PyPI or conda repositories.

Mandatory requirements
----------------------
- python>=3.7
- numpy>=1.15.0
- seaborn (Visual SED images)
- h5py (Used to compress large trajectories files)
- scipy (For fast Lorentz fitting)

Optional GPU requirements
-------------------------
The default pySED installation is CPU-only.  The scattering kernels can use a
CuPy backend when the user explicitly requests it and installs a CUDA-specific
CuPy wheel:

- ``pip install ".[gpu-cuda12x]"`` for CUDA 12.x
- ``pip install ".[gpu-cuda11x]"`` for CUDA 11.x

GPU acceleration is opt-in because CuPy wheels depend on the local CUDA runtime.
Use ``backend="cupy"`` in the DSF/current-correlation APIs after installing the
matching CuPy package.  Use ``backend="cpu"`` or omit the backend argument for
the standard NumPy/SciPy path.


