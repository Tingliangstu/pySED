"""Single source for SED data, plot, fit, and lifetime output paths."""

from pathlib import Path


FITTING_QPOINT_DIRECTORY = Path("Fitting-Qpoint")
LIFETIME_DIRECTORY = Path("Lifetime")


def sed_data_files(prefix):
    """The four files required to reload one complete SED result."""
    return tuple(Path(str(prefix) + suffix) for suffix in (
        ".SED", ".Qpts", ".THz", ".Q_distances_and_labels",
    ))


def partial_sed_data(prefix, element, direction):
    """One element-labelled Cartesian component of the SED."""
    return (
        Path(str(prefix) + "_partial_SED")
        / f"{Path(prefix).name}.SED_{element}_{direction}"
    )


def dispersion_figure(params):
    """The selected total or element-resolved dispersion plot."""
    if getattr(params, "plot_partial_SED", False):
        direction = params.plot_partial_dir or "xyz"
        return (
            Path(params.out_files_name + "_partial_SED")
            / f"SED_{params.plot_partial_element}_{direction}.png"
        )
    return Path(params.out_files_name + "-SED.png")


def slice_figure(params, *, fitting=False):
    """A plain single-Q spectrum or an explicit fit overlay."""
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
    return FITTING_QPOINT_DIRECTORY / f"Fitting-{q_index}-qpoint.png"


def qpoint_lifetime_data(q_index):
    return LIFETIME_DIRECTORY / f"Fitting-{q_index}-qpoint.Fre_lifetime"


def combined_lifetime_data():
    return LIFETIME_DIRECTORY / "Fitting-All-Qpoints.Fre_lifetime"


def lifetime_summary_figure(params=None):
    """Return the summary figure beside the selected SED figure."""
    filename = "Fitting-Frequency-Lifetime.png"
    if params is not None and getattr(params, "plot_partial_SED", False):
        return Path(params.out_files_name + "_partial_SED") / filename
    return Path(filename)
