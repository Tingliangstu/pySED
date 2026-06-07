import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from pySED.compare_experiment import ScatteringMap, prepare_map_comparison
from pySED.scattering_plot import plot_linecuts, plot_map_comparison, plot_scattering_map


def _comparison():
    simulated = ScatteringMap(
        q_axis=np.array([0.0, 0.5, 1.0]),
        energy_axis=np.array([0.0, 1.0, 2.0, 3.0]),
        intensity=np.array(
            [
                [0.0, 1.0, 2.0, 1.0],
                [1.0, 3.0, 4.0, 2.0],
                [0.5, 1.5, 2.5, 1.0],
            ]
        ),
    )
    experimental = ScatteringMap(
        q_axis=simulated.q_axis,
        energy_axis=simulated.energy_axis,
        intensity=simulated.intensity * 0.8,
    )
    return prepare_map_comparison(simulated, experimental)


def test_plot_scattering_map_saves_png(tmp_path):
    comparison = _comparison()
    path = tmp_path / "map.png"

    fig, ax, image = plot_scattering_map(comparison.simulated_map, save_path=path)

    assert path.exists()
    assert path.stat().st_size > 0
    assert image.get_array().shape == (4, 3)
    assert ax.get_xlabel() == "Q path"
    fig.clear()


def test_plot_linecuts_saves_png(tmp_path):
    comparison = _comparison()
    path = tmp_path / "linecuts.png"

    fig, ax = plot_linecuts(
        comparison.simulated_map,
        experimental_map=comparison.experimental_map,
        q_indices=[0, 2],
        save_path=path,
    )

    assert path.exists()
    assert path.stat().st_size > 0
    assert len(ax.lines) == 4
    fig.clear()


def test_plot_map_comparison_saves_summary_png(tmp_path):
    comparison = _comparison()
    path = tmp_path / "comparison.png"

    fig = plot_map_comparison(comparison, save_path=path, q_indices=[1])

    assert path.exists()
    assert path.stat().st_size > 0
    assert len(fig.axes) >= 4
    fig.clear()


def test_plot_linecuts_requires_aligned_maps():
    comparison = _comparison()
    experimental = ScatteringMap(
        q_axis=np.array([0.0, 0.5]),
        energy_axis=comparison.experimental_map.energy_axis,
        intensity=np.ones((2, 4)),
    )

    with pytest.raises(ValueError, match="aligned"):
        plot_linecuts(comparison.simulated_map, experimental_map=experimental)
