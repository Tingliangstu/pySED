# Changelog

## v2.4.0 - 2026-06-04

### ✨ Highlights

- pySED 2.4.0 improves commensurate q-path handling, fixes triclinic supercell box handling for LAMMPS-style workflows, and ships a much more complete online manual.

### ⚡ Core Fixes

- Added support for fractional `q_path` coordinates such as `1/3` in `input_SED.in`.
- Added tolerant handling of short decimal q-path approximations, so values close to small-denominator high-symmetry points, such as `0.33333`, are treated as the intended commensurate fractions. This fixes q-path precision issues around commensurate-point detection (#41).
- Fixed incorrect supercell box handling for LAMMPS-style restricted triclinic cells by deriving the full supercell matrix, converting it to the restricted representation, and applying the matching coordinate rotation for LAMMPS data output (#48).
- Improved primitive-cell consistency checks when comparing the expected supercell box with the trajectory box.

### 🧩 Workflow and Compatibility

- Fixed and clarified links in the MoS2 starting guide and example documentation.
- Expanded guidance for partial SED workflows, including atom-type and direction-resolved examples linked from issue #39.
- Expanded LO-TO splitting guidance with trajectory-dependence examples linked from issue #31.

### 🎨 Documentation and Examples

- Reworked the ReadTheDocs manual structure with clearer navigation for introduction, requirements, installation, starting guide, tutorials, examples, input parameters, theory, tips, publications, reference, and troubleshooting.
- Added a full `input_SED.in` parameter guide to the online manual, with dedicated pages for parameter syntax, defaults, examples, and related settings.
- Added theory and practical tips pages to the documentation, including SED formulas, commensurate q-point construction notes, Lorentzian lifetime notes, LO-TO splitting guidance, and run-script workflow tips.
- Expanded example documentation for CNT, graphene, MoS2, silicon, phonopy references, legacy LAMMPS workflows, and Jupyter tutorials.
- Updated the root README to point users to the expanded online manual and parameter guide.
- Updated the publication list with additional 2026 papers and corrected publication metadata.
