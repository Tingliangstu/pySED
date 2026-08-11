# Primitive-cell silicon SED with GPUMD

This three-dimensional example uses 16,000 atoms in a `20 x 20 x 20`
supercell and the G-X-U-K-G-L path.

```bash
cd structure
python generate_lammps_data.py

cd ../gpumd_run
gpumd

cd ../SED
pySED input.in
```

The authoritative sampling metadata is `total_num_steps = 100000` and
`output_data_stride = 20`, matching the supplied GPUMD production input.
The requested 5,000 frames divide evenly into five blocks. Q-path fractions
are written exactly.

`SED/compare_LD/` is an optional lattice-dynamics comparison. Historical
v2.4 SED and comparison outputs are under
`../legacy/pySED-2.4/outputs/Silicon_primitive_gpumd/`.
