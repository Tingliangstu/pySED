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

"""Declarative input-parameter table for the pySED 3.x SED workflow."""

from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from pySED.structure.atoms import atomic_masses


@dataclass(frozen=True)
class Parameter:
    """One input parameter: its default value and conversion function."""

    default: object
    read: Callable


# ==================== Value readers ====================

def _one(values):
    """Require and return exactly one parameter value."""
    if len(values) != 1:
        raise ValueError("requires exactly 1 value")
    return values[0]


def read_text(values):
    """Read one string, removing surrounding quote characters."""
    return _one(values).strip("'\"")


def read_lower_text(values):
    """Read one string and convert it to lowercase."""
    return read_text(values).lower()


def read_int(values):
    """Read one integer parameter value."""
    return int(_one(values))


def read_float(values):
    """Read one floating-point parameter value."""
    return float(_one(values))


def read_bool(values):
    """Read a case-insensitive True/False value or its 1/0 equivalent."""
    value = _one(values).lower()
    if value not in {"0", "1", "true", "false"}:
        raise ValueError("must be 0 or 1, or True or False")
    return value in {"1", "true"}


def read_optional_float(values):
    """Read one floating-point value or the explicit None value."""
    value = _one(values)
    return None if value.lower() == "none" else float(value)


def read_int3(values):
    """Read exactly three integers into an array."""
    if len(values) != 3:
        raise ValueError("requires 3 values")
    return np.array(values, dtype=int)


def read_matrix3(values):
    """Read nine numbers into a 3x3 row-wise matrix."""
    if len(values) != 9:
        raise ValueError("requires 9 values")
    return np.array(values, dtype=float).reshape(3, 3)


def parse_float_or_fraction(value, max_simple_denominator=12, atol=1.0e-5):
    """Parse a q-path coordinate without losing an explicit fraction.

    Decimal spellings close to a familiar small-denominator fraction are
    normalized (``0.33333`` becomes ``1/3``); all other finite decimal input
    is retained exactly as written by :class:`fractions.Fraction`.
    """

    try:
        parsed = Fraction(value)
    except ZeroDivisionError as error:
        raise ValueError("fraction denominator must not be zero") from error
    if "/" in value:
        return parsed
    simple = parsed.limit_denominator(max_simple_denominator)
    if abs(parsed - simple) <= Fraction(str(atol)):
        return simple
    return parsed


def read_q_path(values):
    """Read q-path coordinates as exact or simplified fractions."""
    return np.array(
        [parse_float_or_fraction(value) for value in values],
        dtype=object,
    )


def read_partial_sed(values):
    """Read ``plot_partial_SED = element [x|y|z]``."""

    if len(values) not in {1, 2}:
        raise ValueError("requires an element and optional x, y, or z")

    element = values[0].strip("'\"").capitalize()
    if element not in atomic_masses:
        raise ValueError(
            f"element '{values[0]}' was not found in the periodic table"
        )

    direction = values[1].strip("'\"").lower() if len(values) > 1 else None
    if direction not in {None, "x", "y", "z"}:
        raise ValueError("direction must be x, y, or z")
    return True, element, direction


# ==================== Input parameters ====================

PARAMETERS = {
    # Control
    "action": Parameter("thinking", read_lower_text),
    "device": Parameter("cpu", read_lower_text),

    # Trajectory and output
    "trajectory_file": Parameter("dump.xyz", read_text),
    "out_files_name": Parameter("pySED", read_text),
    "legacy_unit": Parameter(False, read_bool),
    "lammps_unit": Parameter("metal", read_lower_text),
    "trajectory_read_mode": Parameter("cache", read_lower_text),
    "trajectory_prefetch": Parameter(True, read_bool),
    "netcdf_compression_level": Parameter(1, read_int),
    "netcdf_batch_size": Parameter(64, read_int),

    # Trajectory sampling
    "time_step": Parameter(0.0, read_float),
    "output_data_stride": Parameter(0, read_int),

    # Computation
    "num_blocks": Parameter(5, read_int),
    "max_cores": Parameter(4, read_int),

    # Structure
    "prim_unitcell": Parameter(None, read_matrix3),

    # MD simulation
    "num_atoms": Parameter(0, read_int),
    "total_num_steps": Parameter(0, read_int),

    # Structure
    "basis_lattice_file": Parameter("basis.in", read_text),
    "supercell_dim": Parameter(np.array([1, 1, 1]), read_int3),
    "dump_group_supercell": Parameter(None, read_matrix3),
    "prim_axis": Parameter(None, read_matrix3),
    "rescale_prim": Parameter(True, read_bool),

    # Q-points
    "num_qpaths": Parameter(1, read_int),
    "q_path_name": Parameter("GA", read_text),
    "q_path": Parameter(None, read_q_path),

    # Plot
    "plot_cutoff_freq": Parameter(None, read_optional_float),
    "plot_interval": Parameter(5.0, read_float),
    "plot_slice": Parameter(False, read_bool),
    "qpoint_slice_index": Parameter(0, read_int),
    "plot_color": Parameter("RdBu_r", read_text),
    "colorbar_min": Parameter(None, read_optional_float),
    "colorbar_max": Parameter(None, read_optional_float),
    "if_show_figures": Parameter(False, read_bool),
    "use_contourf": Parameter(False, read_bool),
    "output_partial": Parameter(False, read_bool),
    "plot_partial_SED": Parameter(
        (False, None, None),
        read_partial_sed,
    ),

    # Spectral peak fit
    "lorentz_fit_all_qpoint": Parameter(False, read_bool),
    "lorentz_fit_freq_min": Parameter(None, read_optional_float),
    "lorentz_fit_freq_max": Parameter(None, read_optional_float),
    "fitting_function": Parameter("auto", read_lower_text),
    "peak_min_significance": Parameter(4.0, read_float),
    "initial_guess_hwhm": Parameter(0.001, read_float),
    "peak_max_hwhm": Parameter(1e6, read_float),
    "modulate_factor": Parameter(0, read_int),
    "re_output_total_freq_lifetime": Parameter(False, read_bool),
}


supported_devices = {"cpu", "gpu"}
supported_actions = {"thinking", "compute", "plot", "fit"}


# ==================== Validation ====================

def validate_common(params):
    """Check workflow choices, trajectory sampling, and resource settings.

    Compute/thinking actions require positive time_step and output_data_stride;
    plot/fit may leave those unset. Validate device, trajectory mode, LAMMPS
    units, block/worker counts, and NetCDF tuning values. Invalid settings
    raise ValueError; SED geometry and fitting checks are in validate_sed.
    """
    if params.action not in supported_actions:
        supported = ", ".join(sorted(supported_actions))
        raise ValueError(
            f"action '{params.action}' is not available; supported actions: "
            f"{supported}"
        )
    if params.device not in supported_devices:
        supported = ", ".join(sorted(supported_devices))
        raise ValueError(
            f"device '{params.device}' is not available; "
            f"supported devices: {supported}"
        )
    if params.lammps_unit not in {"metal", "real"}:
        raise ValueError("lammps_unit must be 'metal' or 'real'")
    if params.trajectory_read_mode not in {"cache", "direct"}:
        raise ValueError(
            "trajectory_read_mode must be 'cache' or 'direct'"
        )

    needs_trajectory = params.action in {"thinking", "compute"}
    if needs_trajectory and params.time_step <= 0:
        raise ValueError("time_step must be positive for computation")
    if needs_trajectory and params.output_data_stride <= 0:
        raise ValueError(
            "output_data_stride must be positive for computation"
        )
    if params.time_step < 0:
        raise ValueError("time_step must not be negative")
    if params.output_data_stride < 0:
        raise ValueError("output_data_stride must not be negative")
    if params.num_blocks < 1:
        raise ValueError("num_blocks must be positive")
    if params.max_cores < 1:
        raise ValueError("max_cores must be positive")
    if not 0 <= params.netcdf_compression_level <= 9:
        raise ValueError("netcdf_compression_level must be between 0 and 9")
    if params.netcdf_batch_size < 1:
        raise ValueError("netcdf_batch_size must be positive")


def validate_qpoint_slice_index(params, qpoints=None):
    """Check that qpoint_slice_index is a valid zero-based selection.

    Without loaded qpoints, only reject negative values. When qpoints are
    provided, also enforce the available row count and report the valid range.
    """

    q_index = params.qpoint_slice_index
    if q_index < 0:
        raise ValueError(
            "qpoint_slice_index must be non-negative (zero-based)"
        )
    if qpoints is None:
        return

    num_qpoints = len(np.atleast_2d(qpoints))
    if q_index >= num_qpoints:
        raise ValueError(
            f"qpoint_slice_index = {q_index} is out of range: "
            f"{num_qpoints} q-points are available; valid zero-based "
            f"indices are 0 to {num_qpoints - 1}."
        )


def validate_sed(params):
    """Validate SED geometry, sampling, plotting, and peak-fit controls.

    Compute/thinking require a primitive cell, q path, basis file, and a frame
    count divisible into complete blocks. Plot/fit can use saved spectra without
    these computation inputs. Supplied controls are still validated for every
    action. Local-group cells must be nonsingular and use fixed primitive
    geometry. Invalid combinations raise ValueError.

    This function also unpacks params.plot_partial_SED from its parsed tuple
    into a boolean plus plot_partial_element and plot_partial_dir attributes.
    A selected element is checked against masses in basis.in.
    """

    needs_compute = params.action in {"thinking", "compute"}
    if needs_compute and params.num_atoms <= 0:
        raise ValueError("num_atoms must be positive for SED computation")
    if needs_compute and params.total_num_steps <= 0:
        raise ValueError(
            "total_num_steps must be positive for SED computation"
        )
    if params.num_atoms < 0:
        raise ValueError("num_atoms must not be negative")
    if params.total_num_steps < 0:
        raise ValueError("total_num_steps must not be negative")
    if params.num_qpaths < 1:
        raise ValueError("num_qpaths must be positive")
    validate_qpoint_slice_index(params)
    if np.any(params.supercell_dim <= 0):
        raise ValueError("supercell_dim values must be positive")
    for name in ("colorbar_min", "colorbar_max"):
        value = getattr(params, name)
        if value is not None and not np.isfinite(value):
            raise ValueError(f"{name} must be finite")
    if (
        params.colorbar_min is not None
        and params.colorbar_max is not None
        and params.colorbar_min >= params.colorbar_max
    ):
        raise ValueError(
            "colorbar_min must be smaller than colorbar_max"
        )
    if (
        not np.isfinite(params.peak_min_significance)
        or params.peak_min_significance <= 0
    ):
        raise ValueError("peak_min_significance must be positive and finite")
    if (
        params.lorentz_fit_freq_min is not None
        and (
            not np.isfinite(params.lorentz_fit_freq_min)
            or params.lorentz_fit_freq_min < 0
        )
    ):
        raise ValueError(
            "lorentz_fit_freq_min must be non-negative and finite"
        )
    if (
        params.lorentz_fit_freq_max is not None
        and (
            not np.isfinite(params.lorentz_fit_freq_max)
            or params.lorentz_fit_freq_max <= 0
        )
    ):
        raise ValueError("lorentz_fit_freq_max must be positive and finite")
    if (
        params.lorentz_fit_freq_min is not None
        and params.lorentz_fit_freq_max is not None
        and params.lorentz_fit_freq_min >= params.lorentz_fit_freq_max
    ):
        raise ValueError(
            "lorentz_fit_freq_min must be smaller than "
            "lorentz_fit_freq_max"
        )
    supported_fitting_functions = {"lorentz", "dho", "auto"}
    if params.fitting_function not in supported_fitting_functions:
        supported = ", ".join(sorted(supported_fitting_functions))
        raise ValueError(
            f"fitting_function must be one of: {supported}"
        )
    if (
        not np.isfinite(params.initial_guess_hwhm)
        or params.initial_guess_hwhm <= 0
    ):
        raise ValueError("initial_guess_hwhm must be positive and finite")
    if (
        not np.isfinite(params.peak_max_hwhm)
        or params.peak_max_hwhm <= 0
    ):
        raise ValueError("peak_max_hwhm must be positive and finite")
    if params.initial_guess_hwhm > params.peak_max_hwhm:
        raise ValueError(
            "initial_guess_hwhm must not exceed peak_max_hwhm"
        )
    if params.modulate_factor < 0:
        raise ValueError("modulate_factor must be non-negative")

    if needs_compute and params.prim_unitcell is None:
        raise ValueError("prim_unitcell is required for SED computation")
    if params.prim_unitcell is not None:
        if not np.all(np.isfinite(params.prim_unitcell)):
            raise ValueError("prim_unitcell values must be finite")
        if np.isclose(np.linalg.det(params.prim_unitcell), 0.0):
            raise ValueError("prim_unitcell must not be singular")
    if params.dump_group_supercell is not None:
        if not np.all(np.isfinite(params.dump_group_supercell)):
            raise ValueError(
                "dump_group_supercell must contain nine valid lattice-vector numbers"
            )
        if np.isclose(np.linalg.det(params.dump_group_supercell), 0.0):
            raise ValueError(
                "dump_group_supercell must define a 3D cell; its three "
                "lattice vectors cannot be coplanar"
            )
        if params.rescale_prim:
            raise ValueError(
                "dump_group_supercell is for a local group trajectory. "
                "Set rescale_prim = 0 and provide its matching prim_unitcell"
            )

    if needs_compute and params.q_path is None:
        raise ValueError("q_path is required for SED computation")
    if params.q_path is not None:
        expected = (params.num_qpaths + 1) * 3
        if params.q_path.size != expected:
            raise ValueError(
                f"q_path requires {expected} values for "
                f"num_qpaths = {params.num_qpaths}"
            )
        if not np.all(np.isfinite(np.asarray(params.q_path, dtype=float))):
            raise ValueError("q_path values must be finite")
        if len(params.q_path_name) != params.num_qpaths + 1:
            raise ValueError(
                "q_path_name requires "
                f"{params.num_qpaths + 1} labels for "
                f"num_qpaths = {params.num_qpaths}"
            )

    if needs_compute:
        if params.total_num_steps % params.output_data_stride:
            raise ValueError(
                "total_num_steps must be divisible by "
                "output_data_stride"
            )
        num_frames = (
            params.total_num_steps // params.output_data_stride
        )
        if num_frames < params.num_blocks:
            raise ValueError(
                "the number of trajectory frames must not be smaller "
                "than num_blocks"
            )
        if num_frames % params.num_blocks:
            raise ValueError(
                "the number of trajectory frames must be divisible by "
                "num_blocks"
            )
        if not Path(params.basis_lattice_file).is_file():
            raise ValueError(
                f"input structure '{params.basis_lattice_file}' "
                "was not found"
            )

    (
        params.plot_partial_SED,
        params.plot_partial_element,
        params.plot_partial_dir,
    ) = params.plot_partial_SED

    if params.plot_partial_SED:
        try:
            masses = np.loadtxt(
                params.basis_lattice_file,
                skiprows=2,
                usecols=3,
            )
        except OSError as error:
            raise ValueError(
                f"input structure '{params.basis_lattice_file}' "
                "was not found"
            ) from error

        unique_masses = np.unique(np.atleast_1d(masses))
        matches = np.flatnonzero(
            np.isclose(
                unique_masses,
                atomic_masses[params.plot_partial_element],
                rtol=0.0,
                atol=5.0e-2,
            )
        )
        if not len(matches):
            raise ValueError(
                f"element '{params.plot_partial_element}' was not found "
                f"in the input structure '{params.basis_lattice_file}'"
            )

__all__ = [
    "PARAMETERS",
    "parse_float_or_fraction",
    "validate_common",
    "validate_sed",
    "validate_qpoint_slice_index",
]
