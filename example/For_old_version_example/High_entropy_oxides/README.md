# High-entropy oxide SED with LAMMPS

This migrated case accompanies Yao *et al.*, “Thermal conduction mechanism of
ferroelastic Zr–Y–Nb–Ta–Yb–O high-entropy oxides with glass-like thermal
conductivity,” *Journal of the American Ceramic Society* (2022),
[doi:10.1111/jace.18374](https://doi.org/10.1111/jace.18374).

Run `lammps_run/in.lammps` to create one sorted `trajectory.dump`, then:

```bash
cd SED
pySED input.in
```

The migrated input fixes the former relative paths, supplies the
`50 x 2 x 2` supercell, and replaces the malformed trailing input text and
v2.4-only keys. Peak-derived lifetimes for this strongly anharmonic system
should be interpreted cautiously.

Historical v2.4 outputs are under
`../../legacy/pySED-2.4/outputs/For_old_version_example/High_entropy_oxides/`.
