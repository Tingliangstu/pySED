# Deterministic pySED validation example

This self-contained case needs Python and pySED, but no MD engine. It creates
eight atoms in a `2 x 2 x 2` simple-cubic supercell and 256 frames of a
damped plane wave at reduced Q = (1/2, 0, 0). The known center is 6.25 THz.

```bash
python generate_trajectory.py
pySED input.in
pySED input-plot.in
pySED input-fit.in
```

Run `pySED input-direct.in` to repeat the calculation without a NetCDF cache.
The cache and direct numerical results should agree to the tolerance in
`expected.json` (compare their arrays after accounting for the output prefix).

The four blocks each contain 64 frames, giving 1.5625 THz block resolution.
The expected peak tolerance is therefore one block-frequency bin. The
trajectory is generated rather than committed, so reruns are deterministic and
the repository stays small.
