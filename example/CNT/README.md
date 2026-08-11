# Carbon nanotube SED with GPUMD

This one-dimensional example uses 17,920 atoms in a `1 x 1 x 160`
supercell and follows the G-A path.

```bash
cd structure
python generate_gpumd_xyz.py

cd ../gpumd_run
gpumd

cd ../SED
pySED input.in
```

The GPUMD run writes `dump.xyz` every 8 MD steps, matching
`output_data_stride = 8`. The supplied input starts with
`action = compute`; change it to `plot` or `fit` only after the numerical
SED files exist. The all-Q fit is enabled as an advanced demonstration.

The former v2.4 `CNT.SED`, frequency/Q files, and figure are archived beneath
`../legacy/pySED-2.4/outputs/CNT/SED/`.
