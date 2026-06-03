# Bulk Silicon SED Example

This example calculates the phonon SED of bulk silicon with GPUMD and pySED,
then compares the result with a lattice-dynamics reference.

## Purpose

Use this example to learn how to:

- run a three-dimensional crystalline SED workflow,
- use a multi-segment high-symmetry q-path,
- check primitive-cell and supercell consistency,
- compare pySED output with NEP-driven lattice dynamics,
- perform single-q-point Lorentzian fitting.

## Folder Layout

- `structure/`
  - `POSCAR_prim`: primitive silicon structure.
  - `generate_lammps_data.py`: creates structure files and `basis.in`.
- `gpumd_run/`
  - `run.in`: GPUMD run file.
  - `nep.txt`: NEP potential.
- `SED/`
  - `input_SED.in`: pySED control file.
  - `compare_LD/`: scripts and data for lattice-dynamics comparison.

## Step 1: Generate Structure Files

```bash
cd example/Silicon_primitive_gpumd/structure
python generate_lammps_data.py
```

The script reads `POSCAR_prim`, builds a `20 x 20 x 20` supercell, and writes
the `basis.in` file used by pySED.

## Step 2: Run GPUMD

```bash
cd ../gpumd_run
gpumd
```

The production block uses:

```text
ensemble       nve
dump_exyz      10     1
run            500000
```

The matching pySED settings are:

```text
total_num_steps = 500000
time_step = 1
output_data_stride = 10
dump_xyz_file = '../gpumd_run/dump.xyz'
```

## Step 3: Compute or Plot SED

```bash
cd ../SED
pysed input_SED.in
```

Use `plot_SED = 0` for the first compute run. Use `plot_SED = 1` after
`silicon.SED`, `silicon.Qpts`, and `silicon.THz` exist.

The q-path is:

```text
num_qpaths = 5
q_path_name = 'GXUKGL'
q_path = 0.0 0.0 0.0  0.5 0.0 0.5  0.625 0.25 0.625  0.375 0.375 0.75  0.0 0.0 0.0  0.5 0.5 0.5
```

## Step 4: Compare with Lattice Dynamics

The `SED/compare_LD/` directory contains:

- `get_phonon_dispersion.py`: calculates the NEP-driven lattice-dynamics
  dispersion using calorine and phonopy.
- `plot_phonon_dis_NEP_SED.py`: overlays the reference branches on pySED output.

The expected comparison figure is `SED/compare_LD/Silicon.png`.

## Step 5: Fit a q-Point

The example includes a single-q fitting result for q-point index 2. To repeat or
adjust the fit:

```text
plot_SED = 1
plot_slice = 1
qpoint_slice_index = 2
lorentz = 1
lorentz_fit_cutoff = 20
```

Tune `peak_height`, `peak_prominence`, and `initial_guess_hwhm` until the fitted
peaks match the visible SED slice.

## Checks

- `num_atoms = 16000` must match the trajectory and `basis.in`.
- `supercell_dim = 20 20 20` must match the generated supercell.
- `prim_unitcell` must be consistent with the primitive silicon cell.
- If the generated q-point count is unexpected, check `q_path`, `q_path_name`,
  and the supercell size.
