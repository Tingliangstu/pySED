[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)
![pySED Logo](https://github.com/Tingliangstu/pySED/blob/main/docs/source/_static/logo.png)

# pySED

pySED is a Python package for calculating kinetic-energy-weighted phonon
spectral energy density (SED) from molecular dynamics trajectories. It supports
GPUMD and LAMMPS trajectory formats, plots phonon SED maps, and fits SED peaks
with Lorentzian functions to estimate phonon lifetimes.

Online manual: <https://pysed.readthedocs.io/en/latest/>

## Main Features

- Compute phonon SED from GPUMD or LAMMPS molecular dynamics trajectories.
- Construct commensurate q-points from primitive-cell and supercell information.
- Plot publication-quality SED maps.
- Fit SED peaks at a selected q-point or across all q-points.
- Output frequency-lifetime data from Lorentzian fitting.
- Compute and plot partial SED by atom type and Cartesian direction.
- Run SED calculations in parallel for faster q-point processing.

## Installation

Installing from source is recommended because pySED is actively maintained.

```bash
git clone https://github.com/Tingliangstu/pySED.git
cd pySED
pip install .
```

You can also install directly from GitHub:

```bash
pip install git+https://github.com/Tingliangstu/pySED.git
```

Verify the installation:

```bash
pysed -h
```

or

```bash
pySED -h
```

## Documentation

Start with these manual sections:

- Theory: <https://pysed.readthedocs.io/en/latest/theory.html>
- Quick Start: <https://pysed.readthedocs.io/en/latest/starting.html>
- Input Parameters: <https://pysed.readthedocs.io/en/latest/input_parameters.html>
- Examples: <https://pysed.readthedocs.io/en/latest/example.html>

The `pysed -h` command also lists the supported `input_SED.in` parameters.

## Basic Workflow

1. Generate the MD supercell and `basis.in` file.
2. Run GPUMD or LAMMPS to output coordinates and velocities.
3. Edit `input_SED.in`.
4. Run pySED in compute mode with `plot_SED = 0`.
5. Run pySED in plotting mode with `plot_SED = 1`.
6. Optionally fit Lorentzian peaks for phonon lifetimes.

## Examples

The example library is available at
<https://github.com/Tingliangstu/pySED/tree/main/example>.

### 1D Systems

- Carbon nanotube: <https://github.com/Tingliangstu/pySED/tree/main/example/CNT>

### 2D Systems

- In-plane graphene:
  <https://github.com/Tingliangstu/pySED/tree/main/example/In_plane_graphene_gpumd>
- MoS<sub>2</sub> out-of-plane modes:
  <https://github.com/Tingliangstu/pySED/tree/main/example/MoS2_gpumd>

### 3D Systems

- Bulk silicon:
  <https://github.com/Tingliangstu/pySED/tree/main/example/Silicon_primitive_gpumd>

## Citations

If you use pySED in your research, please cite:

[1] Ting Liang, Wenwu Jiang, Ke Xu, Hekai Bu, Zheyong Fan, Wengen Ouyang, and
Jianbin Xu, "PYSED: A tool for extracting kinetic-energy-weighted phonon
dispersion and lifetime from molecular dynamics simulations", *Journal of
Applied Physics* **138**, 075101 (2025).
<https://doi.org/10.1063/5.0278798>

[2] J. A. Thomas, J. E. Turney, R. M. Iutzi, C. H. Amon, and A. J. H.
McGaughey, "Predicting phonon dispersion relations and lifetimes from the
spectral energy density", *Physical Review B* **81**, 081411 (2010).
<https://doi.org/10.1103/PhysRevB.81.081411>

## Publications Using pySED

See <https://github.com/Tingliangstu/pySED/tree/main/publications>.

## Contact

Ting Liang: liangting.zj@gmail.com

Wenwu Jiang: wwjiang96@163.com

Questions and suggestions can also be raised through the pySED GitHub issues
page.
