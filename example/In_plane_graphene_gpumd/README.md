# In-plane graphene SED with GPUMD

This example computes the G-M-K-G dispersion of a `40 x 40 x 1` graphene
supercell.

```bash
cd structure
python generate_gpumd_xyz.py

cd ../gpumd_run
gpumd

cd ../SED
pySED input.in
```

`output_data_stride = 10` matches the GPUMD output interval. The K-point
coordinates are written as exact thirds in `input.in`. After computation,
set `action = plot` to render the map. The `SED/compare_LD/` scripts are
optional harmonic-reference tools and are not required by pySED.

Historical v2.4 results are under
`../legacy/pySED-2.4/outputs/In_plane_graphene_gpumd/`.
