# pySED Examples

This directory contains example workflows for calculating phonon spectral
energy density (SED) with pySED. Reproducing at least one example is strongly
recommended before applying pySED to a new system.

## Standard Workflow

Each modern GPUMD example follows the same sequence:

1. Generate a GPUMD `model.xyz` supercell and a pySED `basis.in` file.
2. Run GPUMD to generate a production `dump.xyz` trajectory with coordinates
   and velocities.
3. Run pySED in compute mode with `plot_SED = 0`.
4. Run pySED again in plotting mode with `plot_SED = 1`.
5. Optionally tune Lorentzian peak fitting and compare with lattice dynamics.

## Recommended Examples

### 1D systems

- `CNT`
  - Carbon nanotube example.
  - Teaches q-path setup along a one-dimensional axis and all-q-point lifetime
    fitting.

### 2D systems

- `In_plane_graphene_gpumd`
  - In-plane graphene dispersion with GPUMD.
  - Includes a lattice-dynamics comparison workflow in `SED/compare_LD`.

- `MoS2_gpumd`
  - Low-frequency out-of-plane modes in layered MoS2.
  - Includes Lorentzian fitting outputs for the low-frequency range.

### 3D systems

- `Silicon_primitive_gpumd`
  - Bulk silicon with a multi-segment q-path.
  - Includes a lattice-dynamics comparison workflow in `SED/compare_LD`.

## Reference Dispersion

`Ref_Phonon_dispersion_from_phonopy` contains reference lattice-dynamics
workflows using phonopy. These are useful for checking whether the SED branches
agree with a harmonic reference dispersion.

## Older Examples

`For_old_version_example` contains older LAMMPS-based and legacy workflows. Use
the modern GPUMD examples above first unless you specifically need one of the
older cases.

## Tips

- Confirm that `num_atoms` matches both the trajectory and `basis.in`.
- Confirm that `output_data_stride` matches the MD dump stride.
- Use `plot_slice = 1` to inspect one q-point before fitting all q-points.
- Use `output_partial = 1` for atom-type and direction-resolved partial SED.
- pySED does not add LO-TO splitting by itself; the effect must already be
  present in the MD trajectory.

Author: pySED development team
