#!/usr/bin/env python
# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

"""Command-line entry point for pySED."""

from __future__ import annotations

import os
import sys

from pySED import __version__
from pySED.parser import read_input
from pySED.pipeline import run as run_pipeline


_LOGO_LINES = (
    r" _ __   _   _    _____  ______  _____",
    "| '_ \\ | | | |  / ____||  ____||  __ \\",
    r"| |_) || |_| | | (___  | |__   | |  | |",
    r"| .__/  \__, |  \___ \ |  __|  | |  | |",
    r"| |      __/ | ____) | | |____ | |__| |",
    r"|_|     |___/ |_____/  |______||_____/",
)
_BANNER_WIDTH = max(len(line) for line in _LOGO_LINES)
_BANNER_INDENT = "    "


LOGO = "\n".join(f"{_BANNER_INDENT}{line}" for line in _LOGO_LINES)


def _print_banner() -> None:
    """Print the original logo with its version below the right-hand side."""
    version = f"v-{__version__}"
    print(f"\n{LOGO}\n{_BANNER_INDENT}{version:>{_BANNER_WIDTH - 5}}\n")


def _configure_console_encoding() -> None:
    """Prevent locale-limited Windows terminals from aborting a run."""

    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(errors="replace")


def _print_run_summary(input_file, params) -> None:
    """Print the selected input, action, backend, and output prefix."""
    if params.backend == "cupy":
        backend = "CuPy (single GPU)"
    elif params.max_cores == 1:
        backend = "NumPy (1 CPU process)"
    else:
        backend = f"NumPy (up to {params.max_cores} CPU processes)"

    border = f"{_BANNER_INDENT}{'=' * _BANNER_WIDTH}"
    divider = f"{_BANNER_INDENT}{'-' * _BANNER_WIDTH}"
    print(border)
    print(f"{_BANNER_INDENT}Starting pySED task")
    print(divider)
    print(f"{_BANNER_INDENT}Input   : {input_file}")
    print(f"{_BANNER_INDENT}Method  : {params.method.upper()}")
    print(f"{_BANNER_INDENT}Action  : {params.action}")
    print(f"{_BANNER_INDENT}Backend : {backend}")
    print(f"{_BANNER_INDENT}Output  : {params.out_files_name}")
    print(f"{border}\n")


def _show_help() -> None:
    """Print command-line usage and exit successfully."""
    _print_banner()
    print(
        """  Usage
    pySED [input_file]
    pysed [input_file]
    pySED --help

  Description
    Calculate phonon spectral energy density (SED) and phonon
    lifetimes from molecular dynamics trajectories.

    Supports extended XYZ (.xyz), LAMMPS custom dump, and NetCDF formats.

    If no input_file is given, defaults to input_SED.in in the current
    directory. To use another file, pass its path on the command line.

  References
    [1] Liang et al., J. Appl. Phys. 138, 075101 (2025).
        https://doi.org/10.1063/5.0278798
    [2] Thomas et al., Phys. Rev. B 81, 081411 (2010).
        https://doi.org/10.1103/PhysRevB.81.081411

    Please cite these papers when using pySED in your research.

  Input parameters (input_SED.in)
    Write one parameter = value per line. Use # for comments.
    Boolean values accept True/False or 1/0.

    Workflow
      method                   : SED method. Default: sed.
      action                   : thinking, compute, plot, or fit.
                                 Default: thinking (next incomplete stage).
      backend                  : numpy (CPU) or cupy (single GPU).
                                 Default: numpy.

    MD simulation (required for computation)
      num_atoms                : Number of atoms in the trajectory.
      total_num_steps          : Number of MD steps.
      time_step                : MD time step in fs.
      output_data_stride       : Trajectory output interval in steps.
      num_blocks               : Number of blocks to average. Default: 5.

    Input and output
      trajectory_file          : Trajectory path. Default: dump.xyz.
      basis_lattice_file       : Basis lattice path. Default: basis.in.
      out_files_name           : Output prefix. Default: pySED.
      lammps_unit              : metal or real for LAMMPS text. Default: metal.
      trajectory_read_mode     : cache or direct. Default: cache.
      trajectory_prefetch      : Prefetch trajectory batches. Default: True.
      netcdf_batch_size        : Frames per text batch. Default: 64.
      netcdf_compression_level : Cache compression level (0-9). Default: 1.
      output_partial           : Save partial SED. Default: False.
      legacy_unit              : Use the temporary 2.x SED scale. Default:
                                 False.

      cache reuses a .pysed.nc file; direct streams the text trajectory.
      Native NetCDF is read directly in either mode.
      legacy_unit is for historical comparisons and will be removed in a
      future major release. Saved SED scale metadata controls plot/fit labels.

    Parallelization
      max_cores                : Maximum CPU processes. Default: 4.

    Crystal structure
      prim_unitcell            : Primitive cell, nine values in Angstrom.
      prim_axis                : Optional primitive-cell transformation.
      supercell_dim            : Replication along three axes. Default: 1 1 1.
      rescale_prim             : Rescale the primitive cell. Default: True.
      dump_group_supercell     : Optional local cell for a dumped atom group.

    Q-points and Brillouin zone
      num_qpaths               : Number of path segments. Default: 1.
      q_path_name              : Labels for path endpoints. Default: GA.
      q_path                   : Endpoints in fractional reciprocal coordinates.
                                 Example: 0 0 0  0.5 0 0.

    Plotting
      plot_cutoff_freq         : Upper plot frequency in THz. Default: None.
      plot_interval            : Frequency tick interval in THz. Default: 5.
      plot_color               : Matplotlib colormap. Default: RdBu_r.
      colorbar_min             : Lower log color limit. Default: None.
      colorbar_max             : Upper log color limit. Default: None.
      use_contourf             : Use filled contours. Default: False.
      plot_slice               : Plot a single q-point spectrum. Default: False.
      qpoint_slice_index       : Zero-based q-point index. Default: 0.
      plot_partial_SED         : Element and optional direction, e.g. Si x.
      if_show_figures          : Show figures interactively. Default: False.

    Peak detection and fitting
      fitting_function         : auto, lorentz, or dho. Default: auto.
      peak_min_significance    : Local-noise significance threshold. Default: 4.
      initial_guess_hwhm       : Initial peak HWHM in THz. Default: 0.001.
      peak_max_hwhm            : Maximum peak HWHM in THz. Default: 1e6.
      lorentz_fit_freq_min     : Lower fit frequency in THz. Default: None.
      lorentz_fit_freq_max     : Upper fit frequency in THz. Default: None.
      lorentz_fit_all_qpoint   : Fit all q-points. Default: False.
      modulate_factor          : Trim samples at each fit-range end. Default: 0.
      re_output_total_freq_lifetime
                                 Rebuild lifetime summary after a single-q
                                 refit. Default: False.

      Tips
        Use action = plot with plot_slice = True and qpoint_slice_index
        to inspect a spectrum before fitting.

        Increase peak_min_significance to reject more weak/noisy peaks;
        decrease it to detect weaker peaks. It is a dimensionless threshold.

        To refit one q-point, use action = fit, set qpoint_slice_index,
        keep lorentz_fit_all_qpoint = False, and adjust the fitting controls.
        Set re_output_total_freq_lifetime = True to update the summary.

  Examples
    pysed                       # Read input_SED.in in the current directory.
    pysed custom.in             # Read a named input file.
    pysed input_SED.in           # action = compute, plot, or fit in the file.
    pysed -h                    # Show this help.
"""
    )
    raise SystemExit(0)


def _detect_input_file() -> str:
    """Resolve the optional input-file argument or handle a help request."""
    if len(sys.argv) == 1:
        return "input_SED.in"

    if len(sys.argv) == 2:
        argument = sys.argv[1]
        if argument.lower() in {"h", "help", "-h", "--help", "-help"}:
            _show_help()
        return argument

    raise ValueError("pySED takes zero or one argument; run 'pySED -h' for help")


def main() -> int:
    """Run the command and return a process exit status."""

    _configure_console_encoding()
    try:
        input_file = _detect_input_file()
        if not os.path.isfile(input_file):
            raise FileNotFoundError(f"input file '{input_file}' was not found")
        _print_banner()
        params = read_input(input_file)
        _print_run_summary(input_file, params)
        run_pipeline(params)
    except (EOFError, ImportError, KeyError, OSError, RuntimeError, ValueError) as error:
        print(f"\nERROR: {error}\n", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted.\n", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
