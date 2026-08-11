# PTCDI-C8 SED with LAMMPS

This historical molecular-crystal case has a pySED 3.0 `SED/input.in`.
`lammps_run/in.lammps` is a plain-text combined-dump block to insert into the
production input contained in `lammps_run/lammps_run.rar`; it writes
`trajectory.dump` with positions and velocities in one sorted file.

The retained RAR asset is unmodified and has SHA-256:

```text
20616ce011c288765ddb3de2e0ed5fe00a4e529441a05baea157d582ed101a97  lammps_run.rar
```

After generating the trajectory:

```bash
cd SED
pySED input.in
```

The migrated cell is `30 x 4 x 4` with 88 basis atoms and 42,240 total
atoms. Historical v2.4 results are under
`../../legacy/pySED-2.4/outputs/For_old_version_example/PTCDI-C8/`.
