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


"""
pySED pipeline — thinking-mode workflow with auto-detect progress.

Thinking mode selects the next needed stage from existing output files.
Explicit actions run the requested stage without an automatic skip.
"""
import os
import sys
import time
from copy import copy
from pathlib import Path

import numpy as np

from pySED.io.netcdf import is_netcdf
from pySED.io.prepare import prepare_trajectory
from pySED.io.schema import CELLS
from pySED.io.trajectory import prefetch_one, trajectory_block_source
from pySED.parameters import validate_qpoint_slice_index
from pySED import (
    FileIO,
    Lorentz,
    OutputPaths,
    Phonon,
    Plot_SED,
    construct_BZ,
)

_FIT_BLOCK_SEPARATOR = "  " + "-" * 58

# ── step detectors ──────────────────────────────────────────────

def _plot_data_ready(params):
    return all(
        path.is_file()
        for path in OutputPaths.sed_data_files(params.out_files_name)
    )


def _require_plot_data(params):
    """Raise a user-facing error when plot/fit prerequisites are missing."""

    required = OutputPaths.sed_data_files(params.out_files_name)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "SED output is incomplete; run action = compute first. Missing: "
            + ", ".join(missing)
        )


def _plot_ready(params):
    """Return whether the images requested by the current input exist."""

    if not OutputPaths.dispersion_figure(params).is_file():
        return False
    return (
        not getattr(params, "plot_slice", False)
        or OutputPaths.slice_figure(params).is_file()
    )


def _fit_ready(params):
    """Return whether the currently requested fit output exists."""

    if params.lorentz_fit_all_qpoint:
        output = OutputPaths.combined_lifetime_data()
    else:
        output = OutputPaths.qpoint_lifetime_data(
            params.qpoint_slice_index
        )
    return output.is_file()


def _print_detail(label, value):
    """Print one aligned detail row for the current pipeline step."""

    print(f"      {label:<16}: {value}")


def _sed_component_description(params):
    """Describe the total or partial SED selected for plotting."""

    if not getattr(params, "plot_partial_SED", False):
        return "total SED"
    direction = params.plot_partial_dir or "x+y+z (summed)"
    return f"{params.plot_partial_element}, direction {direction}"


# ── step runners ────────────────────────────────────────────────

def step_prepare_trajectory(params):
    """Resolve a text trajectory to direct reading or a NetCDF cache."""
    source = Path(params.trajectory_file)
    input_directory = Path(params.input_file).parent
    converted_path = input_directory / f"{source.name}.pysed.nc"

    t0 = time.perf_counter()
    prepared, source_format = prepare_trajectory(
        source=source,
        converted_path=converted_path,
        lammps_unit=params.lammps_unit,
        batch_size=params.netcdf_batch_size,
        compression_level=params.netcdf_compression_level,
        text_mode=params.trajectory_read_mode,
    )
    params.source_format = source_format
    params.trajectory_path = str(prepared)
    elapsed = time.perf_counter() - t0
    if prepared == source and is_netcdf(prepared):
        print(f"  OK Reading NetCDF trajectory: {prepared}\n")
    elif prepared == source:
        print(f"  OK Reading text trajectory directly: {prepared}")
        print("      NetCDF cache         : disabled\n")
    else:
        print(f"  OK Trajectory ready: {prepared} ({elapsed:.1f} s)\n")


def step_compute_sed(params):
    """Compute phonon SED."""
    step_prepare_trajectory(params)
    print("  >> Computing phonon SED ...")
    t0 = time.perf_counter()

    trajectory = trajectory_block_source(params)
    num_frames = params.total_num_steps // params.output_data_stride
    block_size = num_frames // params.num_blocks
    raw_blocks = trajectory.iter_blocks(block_size, params.num_blocks)
    prefetch_enabled = bool(
        params.trajectory_prefetch and params.num_blocks > 1
    )
    if prefetch_enabled:
        prefetch_status = (
            "ON (next block loads in background; may speed up SED)"
        )
    elif params.trajectory_prefetch:
        prefetch_status = "OFF (only one block to process)"
    else:
        prefetch_status = "OFF (blocks load when needed)"
    print(f"      Trajectory prefetch : {prefetch_status}")
    with prefetch_one(raw_blocks, enabled=prefetch_enabled) as blocks:
        try:
            first_block = next(blocks)
        except StopIteration as error:
            raise EOFError("trajectory contains no readable frames") from error

        dump_group_supercell = getattr(params, "dump_group_supercell", None)
        bz_cell = (
            np.asarray(dump_group_supercell, dtype=float)
            if dump_group_supercell is not None
            else first_block.data[CELLS][0]
        )
        BZ_info = construct_BZ.BZ_methods(params, bz_cell)
        sed = Phonon.spectral_energy_density(params, trajectory.info)
        setup_done = time.perf_counter()
        sed.compute_sed(
            params,
            BZ_info,
            first_block,
            blocks,
        )
        calculation_done = time.perf_counter()

    # Flush calculation progress before reporting output files.
    sys.stdout.flush()

    FileIO.write_output(sed, params, BZ_info)
    finished = time.perf_counter()
    elapsed = finished - t0
    print(f"  OK SED compute done ({elapsed:.1f} s)")
    print(f"      Initial read/setup : {setup_done - t0:.1f} s")
    print(
        "      SED calculation    : "
        f"{calculation_done - setup_done:.1f} s"
    )
    print(
        "      Output/finalization: "
        f"{finished - calculation_done:.1f} s"
    )


def step_plot_sed(params):
    """Plot SED dispersion and optional q-slice."""
    _require_plot_data(params)

    data = FileIO.load_data(params)
    if params.plot_slice:
        validate_qpoint_slice_index(params, data.qpoints)

    print("  >> Plotting SED dispersion")
    _print_detail("Component", _sed_component_description(params))
    Plot_SED.plot_bands(data, params)

    if params.plot_slice:
        q_index = params.qpoint_slice_index
        qpoint = data.qpoints[q_index]
        max_frequency, frequency_control = (
            Plot_SED.resolve_slice_frequency_limit(data.freq_fft, params)
        )
        input_name = Path(
            getattr(params, "input_file", "input.in")
        ).name
        frequency_parameter = (
            frequency_control or "plot_cutoff_freq"
        )
        frequency_guidance = (
            f"adjust with {frequency_parameter} in {input_name}"
        )
        print("\n  >> Plotting single-q SED")
        _print_detail("q-point index", f"{q_index} (zero-based)")
        _print_detail(
            "q-point",
            f"({qpoint[0]:.6g}, {qpoint[1]:.6g}, {qpoint[2]:.6g})",
        )
        _print_detail(
            "Frequency range",
            f"0 to {max_frequency:g} THz ({frequency_guidance})",
        )
        _print_detail("Component", _sed_component_description(params))
        _print_detail(
            "Y-axis",
            "SED (eV/THz), logarithmic scale",
        )
        Plot_SED.plot_slice(data, params)

    print("\n  OK SED plotting complete")


def _fit_one_qpoint(data, params):
    """Keep fitting, file output, and plotting as three explicit operations."""
    fit_result = Lorentz.lorentz(data, params)
    lifetime_file = FileIO.write_phonon_lifetime(fit_result, params)
    Plot_SED.plot_slice(data, params, fit_result=fit_result)
    print(f"      Fit figure           : {OutputPaths.slice_figure(params, fitting=True)}")
    print(f"\n  [OK] Lifetime data written: {lifetime_file}")
    return fit_result


def step_fit_sed(params):
    """Fit SED peaks and write linewidth-derived lifetimes."""
    _require_plot_data(params)
    lifetime_figure = None
    lifetime_figure_action = "written"

    if params.lorentz_fit_all_qpoint:
        print("  >> Fitting spectral peaks at all Q-points ...")
        data = FileIO.load_data(params)
        total_qpoints = len(data.q_distances)
        incomplete_qpoints = []
        print(f"      Total Q-points       : {total_qpoints}")
        for j in range(total_qpoints):
            print(f"\n{_FIT_BLOCK_SEPARATOR}")
            print(
                f"  >> Fitting Q-point #{j} "
                f"({j + 1}/{total_qpoints}) ..."
            )
            q_params = copy(params)
            q_params.qpoint_slice_index = j
            q_params.if_show_figures = False
            fit_result = _fit_one_qpoint(data, q_params)
            incomplete_peaks = [
                int(result["peak_number"])
                for result in fit_result.peak_fits
                if result["incomplete_peak_shapes"][0]
            ]
            if incomplete_peaks:
                incomplete_qpoints.append((j, incomplete_peaks))

        total_file, total_modes = FileIO.deal_total_fre_lifetime(
            params,
            total_qpoints,
        )
        if total_modes:
            frequency, lifetime = FileIO.read_lifetime_data(total_file)
            lifetime_figure = Plot_SED.plot_lifetime_summary(
                frequency,
                lifetime,
                params,
            )
        print("\n  >> All-Q fitting summary")
        print(f"      {'Q-points processed':<20}: {total_qpoints}")
        print(f"      {'Fitted modes':<20}: {total_modes}")
        print(f"      {'Lifetime data':<20}: {total_file}")
        print(
            f"      {'Q-point figures':<20}: "
            f"{OutputPaths.FITTING_QPOINT_DIRECTORY}"
        )
        print("\n  OK All-Q spectral fitting complete")
        if incomplete_qpoints:
            incomplete_peak_count = sum(
                len(peaks) for _, peaks in incomplete_qpoints
            )
            peak_word = "peak" if incomplete_peak_count == 1 else "peaks"
            qpoint_word = (
                "Q-point" if len(incomplete_qpoints) == 1 else "Q-points"
            )
            print("\n  >  Q-points to inspect")
            print(
                "      Incomplete line shapes: "
                f"{incomplete_peak_count} {peak_word} across "
                f"{len(incomplete_qpoints)} {qpoint_word}"
            )
            for qpoint_index, peak_numbers in incomplete_qpoints:
                fitted_peak_word = (
                    "peak" if len(peak_numbers) == 1 else "peaks"
                )
                peak_list = ", ".join(str(number) for number in peak_numbers)
                print(
                    f"      Q-point #{qpoint_index:<8}: "
                    f"{fitted_peak_word} {peak_list}"
                )
            print("      Check these figures in Fitting-Qpoint/.")
            print("      Adjust peak_min_significance if needed.")
        print(
            "\n  Tip: Check all fit figures. To refine one Q point, set its"
        )
        print(
            "       qpoint_slice_index, adjust peak_min_significance and"
        )
        print(
            "       fitting_function = auto | lorentz | dho, then set"
        )
        print("       lorentz_fit_all_qpoint = 0 and")
        print("       re_output_total_freq_lifetime = 1.")
    else:
        # single Q-point fit
        print(f"  >> Fitting Q-point #{params.qpoint_slice_index} ...")
        data = FileIO.load_data(params)
        validate_qpoint_slice_index(params, data.qpoints)
        _fit_one_qpoint(data, params)
        if getattr(params, "re_output_total_freq_lifetime", False):
            total_file, total_modes = FileIO.deal_total_fre_lifetime(
                params,
                len(data.q_distances),
            )
            if total_modes:
                frequency, lifetime = FileIO.read_lifetime_data(total_file)
                lifetime_figure = Plot_SED.plot_lifetime_summary(
                    frequency,
                    lifetime,
                    params,
                )
                lifetime_figure_action = "redrawn"
            print(
                "  OK Combined lifetime data rebuilt -> "
                f"{total_file} ({total_modes} modes)"
            )
        print("\n  OK Single-Q spectral fit done")

    if lifetime_figure is not None:
        print(
            f"\n  Lifetime figure {lifetime_figure_action}: {lifetime_figure}"
        )


def run(params):
    """Run one action in the SED-only workflow."""

    method = getattr(params, "method", "sed")
    action = getattr(params, "action", "thinking")

    if method != "sed":
        raise ValueError(
            f"unknown method '{method}'; the supported method is 'sed'"
        )

    steps = {
        "compute": step_compute_sed,
        "plot": step_plot_sed,
        "fit": step_fit_sed,
    }
    if action not in steps and action != "thinking":
        supported = ", ".join(("thinking", *steps))
        raise ValueError(
            f"unknown action '{action}'; supported actions: {supported}"
        )

    if action == "thinking":
        # ── thinking mode: detect what needs doing ──
        _run_thinking(params, steps)
    else:
        # ── explicit action ──
        steps[action](params)

    print("\n  OK pySED task finished successfully.\n")


def _run_thinking(params, steps):
    """Auto-detect progress and run the next needed step(s)."""

    if not _plot_data_ready(params):
        print("  thinking: incomplete SED data -> compute + plot")
        steps["compute"](params)
        steps["plot"](params)
        return

    if not _plot_ready(params):
        print("  thinking: SED data found, plot missing -> plot")
        steps["plot"](params)
        return

    if not _fit_ready(params):
        print("  thinking: SED data and plots ready -> fit")
        steps["fit"](params)
        return

    print("  OK All steps already complete!")
    print("      SED data, requested plots, and lifetime data all exist.")
    print("      To replot, use action = plot.")
    print("      To compute, use action = compute.")
    print("      To fit, use action = fit.")
