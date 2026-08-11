[![Documentation Status](https://readthedocs.org/projects/pysed/badge/?version=latest)](https://pysed.readthedocs.io/en/latest/)
[![CI](https://github.com/Tingliangstu/pySED/actions/workflows/ci.yml/badge.svg)](https://github.com/Tingliangstu/pySED/actions/workflows/ci.yml)

# pySED

pySED calculates kinetic-energy-weighted phonon spectral energy density (SED),
phonon dispersion, linewidths, and lifetimes from molecular-dynamics
trajectories.

Version 3.0 is a breaking upgrade. It provides one validated SED workflow for
GPUMD extended XYZ, LAMMPS custom dumps, and NetCDF trajectories, with direct
streaming or a reusable `.pysed.nc` cache.

## Highlights

- NumPy execution on one or more CPU processes, with shared-memory block
  processing and optional trajectory prefetch.
- Optional CuPy execution on a single NVIDIA GPU.
- Total and element/direction-resolved SED in eV/THz.
- Exact fractional q paths, including values such as `7/15` and `-17/31`.
- Dispersion and single-q plots with linear or logarithmic controls.
- Significance-based peak detection and Lorentz, velocity-DHO, or automatic
  AICc-selected fitting.
- Per-q-point and combined frequency/lifetime outputs.
- Local-group q grids through `dump_group_supercell`.

## Installation

pySED 3.0 requires Python 3.10 or newer.

```bash
git clone https://github.com/Tingliangstu/pySED.git
cd pySED
python -m pip install .
pySED -h
```

For development:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

CuPy is optional and is deliberately not installed automatically. Install the
CuPy wheel that matches the CUDA runtime on the target machine, then set
`backend = cupy`.

## Quick start

Create `input.in` in the analysis directory. A minimal computation contains
the trajectory sampling information, structure, and q path:

```text
method = sed
action = thinking
backend = numpy

trajectory_file = trajectory.dump
trajectory_read_mode = cache
trajectory_prefetch = 1
lammps_unit = metal

num_atoms = 512
total_num_steps = 100000
time_step = 1.0
output_data_stride = 20
num_blocks = 5

basis_lattice_file = basis.in
prim_unitcell = 5.43 0 0  0 5.43 0  0 0 5.43
supercell_dim = 8 8 8
num_qpaths = 1
q_path_name = GX
q_path = 0 0 0  1/2 0 0

out_files_name = Silicon
```

Then run:

```bash
pySED input.in
```

`action = thinking` advances through compute, plot, and fit according to the
outputs that already exist. The stages can also be selected explicitly with
`action = compute`, `action = plot`, or `action = fit`.

The [online manual](https://pysed.readthedocs.io/en/latest/) contains the full
45-parameter reference, trajectory formats, units, output files, fitting
workflow, troubleshooting guide, and v2.4-to-v3 migration table. Runnable
inputs are under [`example/`](https://github.com/Tingliangstu/pySED/tree/main/example).
Outputs from pySED 2.4 are retained only as labelled historical artifacts
under `example/legacy/pySED-2.4/`.

## Citation

If pySED contributes to published work, please cite:

Ting Liang, Wenwu Jiang, Ke Xu, Hekai Bu, Zheyong Fan, Wengen Ouyang, and
Jianbin Xu, “PYSED: A tool for extracting kinetic-energy-weighted phonon
dispersion and lifetime from molecular dynamics simulations,” *Journal of
Applied Physics* **138**, 075101 (2025),
[doi:10.1063/5.0278798](https://doi.org/10.1063/5.0278798).

For the underlying SED method, see J. A. Thomas et al., *Physical Review B*
**81**, 081411 (2010),
[doi:10.1103/PhysRevB.81.081411](https://doi.org/10.1103/PhysRevB.81.081411).

Additional publications using pySED are listed in
[`publications/`](https://github.com/Tingliangstu/pySED/tree/main/publications).

## License and contact

pySED 3.0 and later are licensed under
`GPL-3.0-or-later`; see [`LICENSE`](LICENSE). Cross-project synchronization
provenance is recorded in [`SYNCING.md`](SYNCING.md).

Questions and bug reports are welcome through
[GitHub Issues](https://github.com/Tingliangstu/pySED/issues) or by email at
liangting.zj@gmail.com.
