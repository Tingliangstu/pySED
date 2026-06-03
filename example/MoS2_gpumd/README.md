# MoS2 SED Example

This example calculates the low-frequency out-of-plane SED of layered MoS2
with GPUMD trajectories and pySED.

## Purpose

Use this example to learn how to:

- build a layered-material supercell,
- analyze the `G-A` out-of-plane q-path,
- plot low-frequency SED branches,
- fit Lorentzian peaks in a narrow frequency range.

## Folder Layout

- `structure/`
  - `POSCAR_MoS2`: input structure.
  - `generate_lammps_data.py`: creates structure files and `basis.in`.
- `gpumd_run/`
  - `run.in`: GPUMD run file.
  - `nep.txt`: NEP potential.
- `SED/`
  - `input_SED.in`: pySED control file.
  - Existing SED and Lorentzian fitting outputs.

## Step 1: Generate Structure Files

```bash
cd example/MoS2_gpumd/structure
python generate_lammps_data.py
```

The script reads `POSCAR_MoS2`, builds a `12 x 12 x 16` supercell, and writes
the `basis.in` file needed by pySED.

## Step 2: Run GPUMD

```bash
cd ../gpumd_run
gpumd
```

The production block uses:

```text
ensemble       nve
dump_exyz      50     1
run            500000
```

The matching pySED settings are:

```text
total_num_steps = 500000
time_step = 1
output_data_stride = 50
dump_xyz_file = '../gpumd_run/dump.xyz'
```

## Step 3: Compute SED

```bash
cd ../SED
pysed input_SED.in
```

For the first calculation:

```text
plot_SED = 0
```

The example uses:

```text
num_qpaths = 1
q_path_name = 'GA'
q_path = 0.0 0.0 0.0  0.0 0.0 0.5
```

## Step 4: Plot Low-Frequency Modes

After computing SED, set:

```text
plot_SED = 1
plot_cutoff_freq = 2.0
plot_interval = 0.5
```

The low cutoff keeps the out-of-plane branches visible.

## Step 5: Lorentzian Fitting

For single-q fitting, set:

```text
plot_slice = 1
lorentz = 1
lorentz_fit_cutoff = 2
```

Tune:

- `qpoint_slice_index`
- `peak_height`
- `peak_prominence`
- `initial_guess_hwhm`

After fitting all q-points, pySED writes
`TOTAL-LORENTZ-Qpoints.Fre_lifetime`.

## Expected Result

The final SED image is `SED/bulk_MoS2-SED.png`.

## Checks

- `num_atoms = 13824` must match `basis.in` and the trajectory.
- `supercell_dim = 12 12 16` must match the generated supercell.
- Use a small frequency cutoff when inspecting low-frequency branches.
