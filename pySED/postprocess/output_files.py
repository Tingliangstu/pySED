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

"""Single source for SED data, plot, fit, and lifetime output paths."""

from pathlib import Path


FITTING_QPOINT_DIRECTORY = Path("Fitting-Qpoint")
LIFETIME_DIRECTORY = Path("Lifetime")


def sed_data_files(prefix):
    """Return Paths for SED, reduced q points, THz frequencies, and path labels.

    The tuple order is (.SED, .Qpts, .THz, .Q_distances_and_labels). prefix
    may contain a directory. This function constructs names without creating
    files or directories.
    """
    return tuple(Path(str(prefix) + suffix) for suffix in (
        ".SED", ".Qpts", ".THz", ".Q_distances_and_labels",
    ))


def partial_sed_data(prefix, element, direction):
    """Return the data path for one element and Cartesian direction.

    Use direction x, y, or z. For prefix /work/sample, the file is placed in
    /work/sample_partial_SED/ and named sample.SED_<element>_<direction>.
    """
    return (
        Path(str(prefix) + "_partial_SED")
        / f"{Path(prefix).name}.SED_{element}_{direction}"
    )


def dispersion_figure(params):
    """Return the path for the selected total or element-resolved dispersion plot."""
    if getattr(params, "plot_partial_SED", False):
        direction = params.plot_partial_dir or "xyz"
        return (
            Path(params.out_files_name + "_partial_SED")
            / f"SED_{params.plot_partial_element}_{direction}.png"
        )
    return Path(params.out_files_name + "-SED.png")


def slice_figure(params, *, fitting=False):
    """Return the filename for a plain q-point slice or a fit overlay.

    Plain partial slices include the selected element and direction. Fitting
    uses the shared Fitting-Qpoint directory and the zero-based q index.
    The plotting caller creates directories and writes the figure.
    """
    q_index = params.qpoint_slice_index
    if fitting:
        return fitting_qpoint_figure(q_index)
    if getattr(params, "plot_partial_SED", False):
        direction = params.plot_partial_dir or "xyz"
        return (
            Path(params.out_files_name + "_partial_SED")
            / f"SED_{params.plot_partial_element}_{direction}-{q_index}-qpoint.png"
        )
    return Path(f"SED-{q_index}-qpoint.png")


def fitting_qpoint_figure(q_index):
    """Return the fit-overlay figure path for a zero-based q-point index."""
    return FITTING_QPOINT_DIRECTORY / f"Fitting-{q_index}-qpoint.png"


def qpoint_lifetime_data(q_index):
    """Return the frequency-lifetime data path for one q point."""
    return LIFETIME_DIRECTORY / f"Fitting-{q_index}-qpoint.Fre_lifetime"


def combined_lifetime_data():
    """Return the path for lifetime data combined across all q points."""
    return LIFETIME_DIRECTORY / "Fitting-All-Qpoints.Fre_lifetime"


def lifetime_summary_figure(params=None):
    """Return the lifetime summary path in the current or partial-SED directory."""
    filename = "Fitting-Frequency-Lifetime.png"
    if params is not None and getattr(params, "plot_partial_SED", False):
        return Path(params.out_files_name + "_partial_SED") / filename
    return Path(filename)
