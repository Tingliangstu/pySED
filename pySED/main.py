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
    r"              ____    _____   ____",
    " _ __   _   _/ ___|  | ____| |  _ \\",
    r"| '_ \ | | | \___ \  |  _|   | | | |",
    r"| |_) || |_| |___) | | |___  | |_| |",
    r"| .__/  \__, |____/  |_____| |____/",
    r"|_|     |___/",
)
_BANNER_WIDTH = max(len(line) for line in _LOGO_LINES)
_BANNER_INDENT = "  "


def _center_banner_line(text: str) -> str:
    """Center one line within the command-line banner width."""
    return f"{_BANNER_INDENT}{text:^{_BANNER_WIDTH}}".rstrip()


LOGO = "\n".join(f"{_BANNER_INDENT}{line}" for line in _LOGO_LINES)
TAGLINE = "\n".join(
    _center_banner_line(line)
    for line in (
        f"phonon spectral energy density   v{__version__}",
        "from molecular-dynamics trajectories",
        "Author: liangting.zj@gmail.com",
    )
)


def _print_banner() -> None:
    """Print the program logo, version, and author contact."""
    print(f"\n{LOGO}\n\n{TAGLINE}\n")


def _configure_console_encoding() -> None:
    """Use UTF-8 so status symbols survive Windows terminals and redirected logs."""

    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(encoding="utf-8", errors="replace")


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
    print(f"{_BANNER_INDENT}🚀 Starting pySED task")
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
        """USAGE:
    pySED [input.in]
    pySED -h

DESCRIPTION:
    pySED 3 computes phonon spectral energy density, plots reciprocal-space
    spectra, detects spectral peaks, and fits Lorentz or velocity-DHO line
    shapes. GPUMD extended XYZ, LAMMPS custom dumps, and NetCDF trajectories
    are supported.

CONTROL:
    method = sed
    action = thinking      # run the next incomplete stage
    action = compute       # calculate SED
    action = plot          # plot existing SED output
    action = fit           # fit existing SED output
    backend = numpy        # CPU, optionally using max_cores
    backend = cupy         # optional single-GPU backend
    legacy_unit = False    # True restores the temporary 2.x output scale

legacy_unit is for historical comparisons and will be removed in a future
major release. Saved SED scale metadata controls plot and fit labels.

TRAJECTORY:
    trajectory_file = trajectory.dump
    trajectory_read_mode = cache   # reusable .pysed.nc cache
    trajectory_read_mode = direct  # stream supported text formats

If no path is supplied, pySED reads input.in, then the legacy filename
input_SED.in. New projects should use input.in.

REFERENCE:
    Liang et al., J. Appl. Phys. 138, 075101 (2025),
    https://doi.org/10.1063/5.0278798
"""
    )
    raise SystemExit(0)


def _detect_input_file() -> str:
    """Resolve the optional input-file argument or handle a help request."""
    if len(sys.argv) == 1:
        if os.path.exists("input.in"):
            return "input.in"
        if os.path.exists("input_SED.in"):
            return "input_SED.in"
        raise FileNotFoundError(
            "No input.in or input_SED.in was found in the current directory."
        )

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
        print(f"\n  ✖  ERROR: {error}\n", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\n  ⚠  Interrupted.\n", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
