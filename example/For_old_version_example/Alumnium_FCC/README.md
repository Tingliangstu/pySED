# Aluminum FCC SED with LAMMPS

The directory name retains its historical `Alumnium_FCC` spelling, but this
case has been migrated to pySED 3.0.

```bash
cd structure
python generate_lammps_data.py

cd ../lammps_run
bash run.sh
# or run LAMMPS directly with: lmp -in in.lammps

cd ../SED
pySED input.in
```

`in.lammps` writes one sorted `trajectory.dump` containing
`id type x y z vx vy vz`. The migrated input uses a `30 x 2 x 2`
supercell, 2,097,152 production steps, stride 32, and 64 blocks. Obsolete
`lat_params`, `num_qpoints`, separate position/velocity paths, and HDF5
controls were removed.

Historical v2.4 outputs are under
`../../legacy/pySED-2.4/outputs/For_old_version_example/Alumnium_FCC/`.
