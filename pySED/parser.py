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

from pySED.parameters import common_params, methods, validate_common


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
    """Return ``(line_number, key, values)`` entries from an input file."""

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


def _select_method(entries):
    """Read the method before choosing its parameter table.

    Defaults to ``"sed"`` when no ``method = ...`` line is present,
    matching the ``method`` parameter default in ``common_params``.
    """

    method = common_params["method"].default
    for _, key, values in entries:
        if key == "method":
            method = values[0].strip("'\"").lower()

    if method not in methods:
        supported = ", ".join(methods)
        raise ValueError(
            f"unsupported method '{method}'; supported methods: {supported}"
        )
    return method


def read_input(input_file="input.in"):
    """Parse an input file into a simple parameter namespace."""

    input_path = Path(input_file).expanduser().resolve()
    entries = _read_entries(input_path)
    method = _select_method(entries)
    method_config = methods[method]
    repeated_names = set(common_params) & set(method_config["parameters"])
    if repeated_names:
        names = ", ".join(sorted(repeated_names))
        raise RuntimeError(
            "parameters are defined as both common and method-specific: "
            f"{names}"
        )
    parameter_table = {
        **common_params,
        **method_config["parameters"],
    }
    defaults = {
        key: deepcopy(parameter.default)
        for key, parameter in parameter_table.items()
    }
    params = SimpleNamespace(input_file=str(input_path), **defaults)

    for line_number, key, values in entries:
        parameter = parameter_table.get(key)
        if parameter is None:
            if key in LEGACY_PARAMETER_HINTS:
                raise ValueError(
                    f"{input_path}:{line_number}: parameter '{key}' was "
                    f"removed in pySED 3.0; {LEGACY_PARAMETER_HINTS[key]}"
                )
            raise ValueError(
                f"{input_path}:{line_number}: parameter '{key}' is not "
                f"valid for method '{method}'"
            )

        try:
            value = parameter.read(values)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"{input_path}:{line_number}: invalid '{key}': {error}"
            ) from error
        setattr(params, key, value)

    for key in ("trajectory_file", "basis_lattice_file"):
        if not hasattr(params, key):
            continue
        path = Path(getattr(params, key)).expanduser()
        if not path.is_absolute():
            path = input_path.parent / path
        setattr(params, key, str(path.resolve()))

    try:
        validate_common(params)
        method_config["validate"](params)
    except ValueError as error:
        raise ValueError(f"{input_path}: {error}") from error

    return params


__all__ = ["LEGACY_PARAMETER_HINTS", "read_input"]
