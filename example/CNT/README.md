# CNT SED Example

This example calculates the spectral energy density of a carbon nanotube (CNT)
with GPUMD trajectories and pySED.

## Purpose

Use this example to learn how to:

- build a one-dimensional supercell for SED,
- generate the `basis.in` file required by pySED,
- run GPUMD with `dump_exyz`,
- compute and plot SED along the tube axis,
- fit SED peaks and collect phonon lifetimes.

## Folder Layout

- `structure/`
  - `POSCAR_CNT`: input structure.
  - `generate_gpumd_xyz.py`: creates `model.xyz` and `basis.in`.
- `gpumd_run/`
  - `run.in`: GPUMD run file.
  - `nep.txt`: NEP potential.
- `SED/`
  - `input_SED.in`: pySED control file.
  - Existing SED outputs and the final `CNT-SED.svg` figure.

## Step 1: Generate Structure Files

From the structure folder:

```bash
cd example/CNT/structure
python generate_gpumd_xyz.py
```

The script reads `POSCAR_CNT`, replicates the cell as `1 x 1 x 160`, and writes:

- `model.xyz` for GPUMD,
- `basis.in` for pySED.

Copy or link `model.xyz` into `gpumd_run/` if needed before running GPUMD.

## Step 2: Run GPUMD

From the GPUMD folder:

```bash
cd ../gpumd_run
gpumd
```

The production part of `run.in` uses:

```text
ensemble       nve
dump_exyz      8     1
run            500000
```

This writes a `dump.xyz` trajectory every 8 MD steps. Therefore
`SED/input_SED.in` uses:

```text
total_num_steps = 500000
time_step = 1
output_data_stride = 8
dump_xyz_file = '../gpumd_run/dump.xyz'
```

## Step 3: Compute SED

From the SED folder:

```bash
cd ../SED
pysed input_SED.in
```

For the first calculation, use:

```text
plot_SED = 0
```

pySED writes:

- `CNT.SED`
- `CNT.Qpts`
- `CNT.THz`
- `CNT.Q_distances_and_labels`

## Step 4: Plot and Fit

After the SED data exist, set:

```text
plot_SED = 1
plot_slice = 1
lorentz = 1
```

Use `qpoint_slice_index`, `peak_height`, `peak_prominence`, and
`initial_guess_hwhm` to tune single-q-point fitting. When the single-q fit is
reasonable, use:

```text
lorentz_fit_all_qpoint = 1
```

The all-q fitting workflow writes `TOTAL-LORENTZ-Qpoints.Fre_lifetime`.

## Expected Result

The final SED figure is `SED/CNT-SED.svg`.

## Checks

- `num_atoms = 17920` must match `basis.in` and `dump.xyz`.
- `supercell_dim = 1 1 160` must match the generated supercell.
- The q-path `G-A` is along the tube axis: `q_path = 0 0 0  0 0 0.5`.
