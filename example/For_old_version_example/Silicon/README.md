# Conventional-cell silicon SED with LAMMPS

This historical case has been migrated to the pySED 3.0 trajectory and action
interfaces.

```bash
cd structure
python generate_lammps_data.py

cd ../lammps_run
bash run.sh
# or: lmp -in in.lammps

cd ../SED
pySED input.in
```

`in.lammps` writes one sorted `trajectory.dump` containing positions and
velocities. The v3 input corrects the primitive cell to exactly nine values,
sets `supercell_dim = 16 2 2`, and uses eight blocks so the 65,536 requested
frames divide evenly.

Historical v2.4 maps, lifetime data, and fit figures are under
`../../legacy/pySED-2.4/outputs/For_old_version_example/Silicon/`.
