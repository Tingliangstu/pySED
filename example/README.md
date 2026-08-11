# pySED 3.0 examples

Every SED case uses the v3 control filename `input.in`. The recommended order
is:

1. generate the structure and `basis.in`;
2. run GPUMD or LAMMPS to create one trajectory containing positions,
   velocities, and cell information;
3. run `pySED input.in` with `action = compute`;
4. change the action to `plot` or `fit` and run the same command again.

The GPUMD cases are CNT, in-plane graphene, MoS2, and primitive-cell silicon.
The cases under `For_old_version_example/` originated with earlier pySED
releases but their inputs and LAMMPS dump commands have been migrated to v3.

`deterministic_sed/` is the smallest self-contained validation case. Its
standard-library generator writes a combined LAMMPS custom dump with a known
FFT-bin signal, so no external MD engine is required.

Generated v2.4 figures and numerical files are preserved under
`legacy/pySED-2.4/outputs/` with a manifest and SHA-256 checksums. They are
historical artifacts, not pySED 3.0 expected output.
