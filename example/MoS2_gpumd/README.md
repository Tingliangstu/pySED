# Bulk MoS2 SED with GPUMD

This case emphasizes low-frequency out-of-plane modes of a
`12 x 12 x 16` supercell along G-A.

```bash
cd structure
python generate_lammps_data.py

cd ../gpumd_run
gpumd

cd ../SED
pySED input.in
```

The trajectory interval is 50 MD steps and the plotting range is 0--2 THz.
Compute first, then use `action = plot` or `action = fit`. The v3 fit uses
`peak_min_significance` and `fitting_function = auto`; the removed absolute
height/prominence thresholds are not used.

Historical v2.4 result files are preserved under
`../legacy/pySED-2.4/outputs/MoS2_gpumd/`.
