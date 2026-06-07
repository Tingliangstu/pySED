import numpy as np
from scipy.fft import fft

from pySED.scattering_kernel import (
    circular_autocorrelation,
    density_amplitude,
    density_amplitude_from_phase,
    fft_power,
    frequency_grid,
    phase_factors,
    resolve_backend,
    standard_error,
    to_numpy,
)


def test_resolve_cpu_backend_and_to_numpy():
    backend = resolve_backend("cpu")
    assert backend.name == "cpu"
    arr = backend.xp.asarray([1.0, 2.0])
    np.testing.assert_allclose(to_numpy(arr, backend), [1.0, 2.0])


def test_cupy_backend_is_optional():
    try:
        backend = resolve_backend("cupy")
    except ImportError as exc:
        assert "CuPy backend requested" in str(exc)
    else:
        assert backend.name == "cupy"


def test_density_amplitude_matches_manual_phase_weighting():
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
            [[0.25, 0.0, 0.0], [0.75, 0.0, 0.0]],
        ]
    )
    q_vec = np.array([2 * np.pi, 0.0, 0.0])
    weights = np.array([1.0, 2.0])

    phase = phase_factors(positions, q_vec)
    expected = phase @ weights

    np.testing.assert_allclose(density_amplitude(positions, q_vec, weights), expected)
    np.testing.assert_allclose(density_amplitude_from_phase(phase, weights), expected)


def test_fft_power_and_circular_autocorrelation_are_consistent():
    signal = np.array([1.0 + 0.0j, 0.0 + 1.0j, -1.0 + 0.0j, 0.0 - 1.0j])

    direct = fft_power(signal, norm="forward")
    correlation = circular_autocorrelation(signal)
    via_correlation = fft(correlation, norm="forward").real

    np.testing.assert_allclose(direct, via_correlation)


def test_frequency_grid_and_standard_error():
    grid = frequency_grid(block_size=8, dt=1e-12, positive_only=True)
    np.testing.assert_allclose(grid.frequencies_thz, [0.0, 0.125, 0.25, 0.375])

    samples = np.array([[1.0, 2.0], [3.0, 6.0], [5.0, 10.0]])
    expected = np.std(samples, axis=0, ddof=1) / np.sqrt(3.0)
    np.testing.assert_allclose(standard_error(samples), expected)
