#!/usr/bin/env python
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
    print(f"\n{LOGO}\n\n{TAGLINE}\n")


def _print_run_summary(input_file, params) -> None:
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
