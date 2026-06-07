"""Plotting helpers for paper-ready scattering-map comparisons."""

from pathlib import Path

import numpy as np

from pySED.compare_experiment import MapComparison, ScatteringMap, linecut


def _import_pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required for pySED.scattering_plot") from exc
    return plt


def _map_extent(scattering_map):
    q = np.asarray(scattering_map.q_axis, dtype=float)
    energy = np.asarray(scattering_map.energy_axis, dtype=float)
    return [float(q[0]), float(q[-1]), float(energy[0]), float(energy[-1])]


def _save_if_requested(fig, save_path=None, dpi=300):
    if save_path is None:
        return None
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def plot_scattering_map(
    scattering_map,
    ax=None,
    title=None,
    cmap="magma",
    vmin=None,
    vmax=None,
    q_label="Q path",
    energy_label="Energy",
    colorbar=True,
    save_path=None,
    dpi=300,
):
    """Plot a single q-energy scattering map.

    Returns ``(fig, ax, image)`` so callers can further tune labels, ticks, or
    panel annotations before saving.
    """

    if not isinstance(scattering_map, ScatteringMap):
        raise TypeError("scattering_map must be a ScatteringMap")

    plt = _import_pyplot()
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    else:
        fig = ax.figure

    image = ax.imshow(
        scattering_map.intensity.T,
        origin="lower",
        aspect="auto",
        extent=_map_extent(scattering_map),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xlabel(q_label)
    ax.set_ylabel(energy_label)
    if title:
        ax.set_title(title)
    if colorbar:
        fig.colorbar(image, ax=ax, pad=0.02)

    _save_if_requested(fig, save_path=save_path, dpi=dpi)
    return fig, ax, image


def plot_linecuts(
    simulated_map,
    experimental_map=None,
    q_indices=None,
    ax=None,
    labels=None,
    energy_label="Energy",
    intensity_label="Intensity",
    save_path=None,
    dpi=300,
):
    """Plot energy line cuts at selected q indices."""

    if not isinstance(simulated_map, ScatteringMap):
        raise TypeError("simulated_map must be a ScatteringMap")
    if experimental_map is not None and not isinstance(experimental_map, ScatteringMap):
        raise TypeError("experimental_map must be a ScatteringMap")

    n_q = simulated_map.q_axis.size
    if q_indices is None:
        q_indices = [0, n_q // 2, n_q - 1] if n_q > 2 else list(range(n_q))
    q_indices = [int(index) for index in q_indices]
    if any(index < 0 or index >= n_q for index in q_indices):
        raise ValueError("q_indices must be valid simulated-map q indices")

    if experimental_map is not None:
        if experimental_map.q_axis.size != n_q:
            raise ValueError("experimental_map must be aligned to simulated_map before plotting linecuts")
        if experimental_map.energy_axis.shape != simulated_map.energy_axis.shape:
            raise ValueError("experimental_map energy axis must match simulated_map")

    plt = _import_pyplot()
    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 3.0), constrained_layout=True)
    else:
        fig = ax.figure

    label_prefix = labels or {"simulated": "sim", "experimental": "exp"}
    for index in q_indices:
        q_value = simulated_map.q_axis[index]
        sim_label = "%s q=%.4g" % (label_prefix.get("simulated", "sim"), q_value)
        ax.plot(
            simulated_map.energy_axis,
            linecut(simulated_map.intensity, index, axis=0),
            label=sim_label,
            linewidth=1.6,
        )
        if experimental_map is not None:
            exp_label = "%s q=%.4g" % (label_prefix.get("experimental", "exp"), q_value)
            ax.plot(
                experimental_map.energy_axis,
                linecut(experimental_map.intensity, index, axis=0),
                label=exp_label,
                linewidth=1.0,
                linestyle="--",
            )

    ax.set_xlabel(energy_label)
    ax.set_ylabel(intensity_label)
    ax.legend(frameon=False, fontsize="small")
    _save_if_requested(fig, save_path=save_path, dpi=dpi)
    return fig, ax


def plot_map_comparison(
    comparison,
    q_label="Q path",
    energy_label="Energy",
    cmap="magma",
    residual_cmap="coolwarm",
    q_indices=None,
    figsize=(10.0, 6.2),
    save_path=None,
    dpi=300,
):
    """Create a simulated/experimental/residual map and line-cut figure."""

    if not isinstance(comparison, MapComparison):
        raise TypeError("comparison must be a MapComparison")
    if comparison.simulated_map is None or comparison.experimental_map is None:
        raise ValueError("comparison must include simulated_map and experimental_map")

    plt = _import_pyplot()
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.82))
    ax_sim = fig.add_subplot(grid[0, 0])
    ax_exp = fig.add_subplot(grid[0, 1])
    ax_res = fig.add_subplot(grid[0, 2])
    ax_line = fig.add_subplot(grid[1, :])

    sim_map = ScatteringMap(
        comparison.simulated_map.q_axis,
        comparison.simulated_map.energy_axis,
        comparison.simulated,
        metadata=comparison.simulated_map.metadata,
    )
    exp_map = ScatteringMap(
        comparison.experimental_map.q_axis,
        comparison.experimental_map.energy_axis,
        comparison.experimental,
        metadata=comparison.experimental_map.metadata,
    )
    res_map = ScatteringMap(
        comparison.experimental_map.q_axis,
        comparison.experimental_map.energy_axis,
        comparison.residual,
        metadata={"source": "comparison residual"},
    )

    shared_max = np.nanmax(np.abs(np.concatenate([comparison.simulated.ravel(), comparison.experimental.ravel()])))
    shared_vmax = None if shared_max == 0 else float(shared_max)
    plot_scattering_map(
        sim_map,
        ax=ax_sim,
        title="Simulated",
        cmap=cmap,
        vmin=0.0 if shared_vmax is not None else None,
        vmax=shared_vmax,
        q_label=q_label,
        energy_label=energy_label,
    )
    plot_scattering_map(
        exp_map,
        ax=ax_exp,
        title="Experimental",
        cmap=cmap,
        vmin=0.0 if shared_vmax is not None else None,
        vmax=shared_vmax,
        q_label=q_label,
        energy_label=energy_label,
    )

    res_abs = np.nanmax(np.abs(comparison.residual))
    res_limit = None if res_abs == 0 else float(res_abs)
    plot_scattering_map(
        res_map,
        ax=ax_res,
        title="Residual",
        cmap=residual_cmap,
        vmin=-res_limit if res_limit is not None else None,
        vmax=res_limit,
        q_label=q_label,
        energy_label=energy_label,
    )
    plot_linecuts(
        sim_map,
        experimental_map=exp_map,
        q_indices=q_indices,
        ax=ax_line,
        energy_label=energy_label,
    )

    fig.suptitle(
        "RMSE %.4g, MAE %.4g, corr %.4g"
        % (comparison.rmse, comparison.mae, comparison.correlation),
        fontsize=11,
    )
    _save_if_requested(fig, save_path=save_path, dpi=dpi)
    return fig
