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
from importlib.resources import files
from pathlib import Path

from pySED import __version__
from pySED.parser import read_input
from pySED.pipeline import run as run_pipeline


_MANUAL_URL = "https://pysed.readthedocs.io/"

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
_LOGO_INDENT = _BANNER_INDENT


LOGO = "\n".join(f"{_LOGO_INDENT}{line}" for line in _LOGO_LINES)


def _print_banner() -> None:
    """Print the original logo with its version below the right-hand side."""
    version = f"v-{__version__}"
    print(f"\n{LOGO}\n{_LOGO_INDENT}{version:>{_BANNER_WIDTH - 5}}\n")


def _configure_console_encoding() -> None:
    """Use UTF-8 so status symbols survive Windows terminals and redirected logs."""

    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(encoding="utf-8", errors="replace")


def _print_run_summary(input_file, params) -> None:
    """Print the selected input, action, device, and output prefix."""
    if params.device == "gpu":
        device = "GPU (single CUDA device)"
    elif params.max_cores == 1:
        device = "CPU (1 process)"
    else:
        device = f"CPU (up to {params.max_cores} processes)"

    border = f"{_BANNER_INDENT}{'=' * _BANNER_WIDTH}"
    divider = f"{_BANNER_INDENT}{'-' * _BANNER_WIDTH}"
    print(border)
    print(f"{_BANNER_INDENT}🚀 Starting pySED task")
    print(divider)
    print(f"{_BANNER_INDENT}Input   : {input_file}")
    print(f"{_BANNER_INDENT}Action  : {params.action}")
    print(f"{_BANNER_INDENT}Device  : {device}")
    print(f"{_BANNER_INDENT}Output  : {params.out_files_name}")
    print(f"{border}\n")


def _show_help() -> None:
    """Print command-line usage and exit successfully."""
    _print_banner()
    print(
        f"""    🚀 USAGE:
        pySED [input_file]      or    pysed [input_file]
        pySED -h                or    pysed --help

    ℹ DESCRIPTION:
        pySED is a tool for calculating Phonon Spectral Energy Density (SED)
        and phonon lifetimes from Molecular Dynamics (MD) simulation trajectories.
        It supports GPUMD/LAMMPS trajectories output formats.

        If no input_file is given, defaults to `input_SED.in` in the current directory.

    ⚙ INPUT PARAMETERS (input_SED.in):
        Write one parameter = value per line. Use # for comments.
        Boolean values accept True/False or 1/0.

        🚀 To create a basic input file, run:
            pysed init

        For detailed parameter descriptions and usage, see the manual:
            {_MANUAL_URL}

    ▣ REFERENCE:
        [1] Liang et al., J. Appl. Phys. 138, 075101 (2025). (For pySED paper)
            https://doi.org/10.1063/5.0278798

        [2] Thomas et al., Phys. Rev. B 81, 081411 (2010). (For the SED method)
            https://doi.org/10.1103/PhysRevB.81.081411

        ** Please cite these papers when using pySED in your research. **
"""
    )
    raise SystemExit(0)


def _show_init_help() -> None:
    """Explain input-file generation without writing a file."""
    print(
        f"""    USAGE:
        pysed init

    DESCRIPTION:
        Create a basic input_SED.in in the current directory.
        Edit the example values to match your system, then run pysed.
        An existing input_SED.in is never overwritten.

        Online manual:
            {_MANUAL_URL}
"""
    )
    raise SystemExit(0)


def _generate_input_file() -> None:
    """Create the bundled starter input without overwriting an existing path."""
    template = files("pySED").joinpath("templates", "input_SED.in").read_text(
        encoding="utf-8"
    )
    path = Path("input_SED.in")
    try:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(template)
    except FileExistsError:
        raise FileExistsError(
            "input_SED.in already exists; no changes made. "
            "Edit that file or run 'pysed init' in a new directory."
        ) from None
    print(
        "Created input_SED.in.\n"
        "Edit the example values to match your system, then run pysed.\n"
        f"Online manual: {_MANUAL_URL}"
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
        if argument == "init":
            _generate_input_file()
        return argument

    if (
        len(sys.argv) == 3
        and sys.argv[1] == "init"
        and sys.argv[2].lower() in {"-h", "--help"}
    ):
        _show_init_help()

    raise ValueError(
        "Use 'pysed [input_file]' to run, 'pysed init' to create an input file, "
        "or 'pysed -h' for help."
    )


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
