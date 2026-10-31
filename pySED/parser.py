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

"""Read a pySED ``key = value`` input file."""

import shlex
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

from pySED.parameters import PARAMETERS, validate_common, validate_sed


LEGACY_PARAMETER_HINTS = {
    "compress": "use trajectory_read_mode = cache or direct",
    "dump_xyz_file": "use trajectory_file with one combined trajectory",
    "file_format": "the trajectory format is detected automatically",
    "lorentz": "use action = fit",
    "lorentz_fit_cutoff": (
        "use lorentz_fit_freq_min and lorentz_fit_freq_max"
    ),
    "num_splits": "use num_blocks",
    "output_hdf5": "the .pysed.nc cache path is managed automatically",
    "peak_height": "use peak_min_significance",
    "peak_prominence": "use peak_min_significance",
    "plot_SED": "use action = plot",
    "pos_file": "use trajectory_file with positions and velocities together",
    "use_parallel": "select backend and set max_cores",
    "vels_file": "use trajectory_file with positions and velocities together",
    "with_eigs": "this internal option is no longer an input parameter",
}


def _read_entries(input_file):
    """Read (line_number, key, values) records from a key = value input file.

    Ignore blank lines and comments, retain quoted tokens, and reject duplicate
    keys or malformed assignments. Line numbers are preserved so later conversion
    errors can identify the offending input line.
    """

    path = Path(input_file).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    entries = []
    first_occurrence = {}
    with path.open("r", encoding="utf-8-sig") as stream:
        for line_number, line in enumerate(stream, start=1):
            clean = line.strip()
            if not clean or clean.startswith("#"):
                continue
            if "=" not in clean:
                raise ValueError(
                    f"{path}:{line_number}: expected 'key = value'"
                )

            key, raw_values = clean.split("=", 1)
            key = key.strip()
            try:
                values = shlex.split(
                    raw_values,
                    comments=True,
                    posix=False,
                )
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: {error}"
                ) from error
            if not key or not values:
                raise ValueError(
                    f"{path}:{line_number}: parameter or value is missing"
                )
            if key in first_occurrence:
                raise ValueError(
                    f"{path}:{line_number}: parameter '{key}' is repeated; "
                    f"first defined on line {first_occurrence[key]}"
                )
            first_occurrence[key] = line_number
            entries.append((line_number, key, values))
    return entries


def read_input(input_file="input_SED.in"):
    """Return validated parameters with defaults and input-relative paths resolved.

    Convert each key using PARAMETERS, rejecting unknown keys and reporting
    migration hints for removed 2.x names. Resolve trajectory_file and
    basis_lattice_file relative to the input file's directory, then run common
    and SED validation. Return a SimpleNamespace used by the pipeline.
    """

    input_path = Path(input_file).expanduser().resolve()
    entries = _read_entries(input_path)
    defaults = {
        key: deepcopy(parameter.default)
        for key, parameter in PARAMETERS.items()
    }
    params = SimpleNamespace(input_file=str(input_path), **defaults)

    for line_number, key, values in entries:
        parameter = PARAMETERS.get(key)
        if parameter is None:
            if key in LEGACY_PARAMETER_HINTS:
                raise ValueError(
                    f"{input_path}:{line_number}: parameter '{key}' was "
                    f"removed in pySED 3.0; {LEGACY_PARAMETER_HINTS[key]}"
                )
            raise ValueError(
                f"{input_path}:{line_number}: parameter '{key}' is not "
                "valid for pySED"
            )

        try:
            value = parameter.read(values)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"{input_path}:{line_number}: invalid '{key}': {error}"
            ) from error
        setattr(params, key, value)

    for key in ("trajectory_file", "basis_lattice_file"):
        path = Path(getattr(params, key)).expanduser()
        if not path.is_absolute():
            path = input_path.parent / path
        setattr(params, key, str(path.resolve()))

    try:
        validate_common(params)
        validate_sed(params)
    except ValueError as error:
        raise ValueError(f"{input_path}: {error}") from error

    return params


__all__ = ["LEGACY_PARAMETER_HINTS", "read_input"]
