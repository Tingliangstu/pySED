# Changelog

## v3.0.0 - 2026-08-11

### Breaking interface changes

- Replaced the pySED 2.x parser with declarative parameter conversion and
  validation. The SED-only `action` and `device` workflow is now
  the public interface.
- Replaced split position/velocity inputs with one `trajectory_file` and
  automatic format detection. The default input filename is `input_SED.in`;
  other names must be passed explicitly on the command line.
- Removed 14 v2 parameters. Each removed key now fails with a targeted
  migration hint instead of being silently accepted.
- Changed the cache format/name to `.pysed.nc`, the public SED unit to eV/THz,
  the minimum Python version to 3.10, and the project license to
  GPL-3.0-or-later.

### Computation and trajectory I/O

- Added GPUMD extended-XYZ, LAMMPS combined custom-dump, and NetCDF readers,
  with direct streaming or reusable cache modes.
- Added block prefetch, persistent CPU multiprocessing with shared memory, and
  an optional CuPy single-GPU backend with timing and memory reporting.
- Added energy-preserving one-sided spectra, discrete `dt²` scaling, total and
  partial SED, and consistent eV/THz output.
- Preserved exact fractional q-path input and short-decimal normalization.
- Added local-group q-grid construction through `dump_group_supercell`.

### Fitting, documentation, and quality

- Added significance-based peak detection, Lorentz and velocity-DHO line
  shapes, AICc automatic model selection, small-amplitude-safe fitting, and
  combined lifetime outputs.
- Added a complete flat 45-parameter manual, workflow/output/unit chapters,
  migrated examples, a deterministic synthetic example, and a v2.4 migration
  guide.
- Added cross-platform tests, warnings-as-errors documentation builds, package
  builds, and installed-wheel CLI smoke tests.
- Fixed cache invalidation when `lammps_unit` changes, extended-XYZ property
  reordering, and false-success reporting for failed plot/fit actions.

## v2.4.0 - 2026-06-04

### Highlights

- Improved commensurate q-path handling and added exact fractional q-path
  coordinates such as `1/3`.
- Normalized short decimal approximations near common small-denominator
  high-symmetry points.
- Corrected restricted-triclinic LAMMPS supercell handling and primitive-cell
  consistency checks.
- Expanded the online manual, parameter pages, examples, theory notes, and
  practical workflow guidance.
