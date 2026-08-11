"""Output paths for SED peak fitting and lifetime summaries."""

from pathlib import Path


FITTING_QPOINT_DIRECTORY = Path("Fitting-Qpoint")
LIFETIME_DIRECTORY = Path("Lifetime")


def fitting_qpoint_figure(q_index):
    """Return the fitted SED figure path for one zero-based Q-point."""

    return FITTING_QPOINT_DIRECTORY / f"Fitting-{q_index}-qpoint.png"


def qpoint_lifetime_data(q_index):
    """Return the frequency-lifetime data path for one Q-point."""

    return LIFETIME_DIRECTORY / f"Fitting-{q_index}-qpoint.Fre_lifetime"


def combined_lifetime_data():
    """Return the combined frequency-lifetime data path."""

    return LIFETIME_DIRECTORY / "Fitting-All-Qpoints.Fre_lifetime"


def lifetime_summary_figure(params=None):
    """Return the summary figure beside the selected SED figure."""

    filename = "Fitting-Frequency-Lifetime.png"
    if params is not None and getattr(params, "plot_partial_SED", False):
        return Path(params.out_files_name + "_partial_SED") / filename
    return Path(filename)
