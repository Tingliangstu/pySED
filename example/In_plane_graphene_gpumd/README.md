# In-Plane Graphene SED Example

This example calculates the in-plane phonon SED of graphene with GPUMD and
pySED, then compares the SED map with a lattice-dynamics reference.

## Purpose

Use this example to learn how to:

- build a two-dimensional graphene supercell,
- calculate SED along the `G-M-K-G` path,
- use `use_contourf` for multi-path plotting,
- compare pySED output with NEP-driven lattice dynamics.

## Folder Layout

- `structure/`
  - `POSCAR_graphene`: primitive graphene structure.
  - `generate_gpumd_xyz.py`: creates `model.xyz` and `basis.in`.
- `gpumd_run/`
  - `run.in`: GPUMD run file.
  - `nep.txt`: NEP potential.
- `SED/`
  - `input_SED.in`: pySED control file.
  - `compare_LD/`: scripts and data for lattice-dynamics comparison.

## Step 1: Generate Structure Files

```bash
cd example/In_plane_graphene_gpumd/structure
python generate_gpumd_xyz.py
```

The script reads `POSCAR_graphene`, builds a `40 x 40 x 1` supercell, and
writes:

- `model.xyz` for GPUMD,
- `basis.in` for pySED.

The GPUMD structure uses `pbc="T T F"` because graphene is periodic in-plane.

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

The corresponding pySED settings are:

```text
total_num_steps = 500000
time_step = 1
output_data_stride = 10
dump_xyz_file = '../gpumd_run/dump.xyz'
```

## Step 3: Compute SED

```bash
cd ../SED
pysed input_SED.in
```

For the first run:

```text
plot_SED = 0
```

The q-path is:

```text
num_qpaths = 3
q_path_name = 'GMKG'
q_path = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333333 0.3333333 0.0  0.0 0.0 0.0
```

## Step 4: Plot

After the SED data are written, set:

```text
plot_SED = 1
use_contourf = 1
```

Tune `plot_cutoff_freq`, `plot_interval`, `colorbar_min`, and `colorbar_max`
if the plot contrast is not clear.

## Step 5: Compare with Lattice Dynamics

The `SED/compare_LD/` directory contains:

- `get_phonon_dispersion.py`: calculates the NEP-driven lattice-dynamics
  dispersion using calorine and phonopy.
- `plot_phonon_dis_NEP_SED.py`: overlays the lattice-dynamics branches on the
  pySED SED map.

The expected comparison figure is `SED/compare_LD/Graphene.png`.

## Checks

- `num_atoms = 3200` must match the generated supercell.
- `supercell_dim = 40 40 1` must match `generate_gpumd_xyz.py`.
- `dump_xyz_file` must point to `../gpumd_run/dump.xyz`.
