import numpy as np
from scipy.fft import fft

from pySED.dsf import (
    compute_coherent_dsf_via_correlation,
    compute_current_correlations,
    compute_dsf,
    compute_partial_dsf,
    xray_form_factor_weights,
)


def test_dsf_shapes_and_nonnegative_intensity():
    n_frames = 32
    dt = 1e-12
    positions = np.zeros((n_frames, 2, 3))
    positions[:, 1, 0] = 0.5
    qpoints = np.array([[2 * np.pi, 0.0, 0.0], [0.0, 2 * np.pi, 0.0]])

    result = compute_dsf(
        positions,
        qpoints,
        dt,
        atom_types=["Si", "Si"],
        experiment="neutron",
        num_blocks=2,
    )

    assert result.coherent.shape == (2, 8)
    assert result.incoherent.shape == (2, 8)
    assert result.total.shape == (2, 8)
    assert np.all(result.total >= 0)
    assert np.all(result.elastic_coherent >= 0)


def test_current_correlations_shapes():
    n_frames = 32
    dt = 1e-12
    positions = np.zeros((n_frames, 1, 3))
    velocities = np.zeros_like(positions)
    velocities[:, 0, 0] = 1.0
    qpoints = np.array([[2 * np.pi, 0.0, 0.0]])

    result = compute_current_correlations(positions, velocities, qpoints, dt)

    assert result.longitudinal.shape == (1, 16)
    assert result.transverse.shape == (1, 16)
    assert result.longitudinal_sem.shape == result.longitudinal.shape
    assert result.transverse_sem.shape == result.transverse.shape
    assert result.metadata["uncertainty"] == "standard error of the block-averaged spectra"
    assert result.longitudinal[0, 0] > 0


def test_dsf_matches_direct_density_fft_estimator():
    n_frames = 16
    n_atoms = 3
    dt = 2e-13
    time = np.arange(n_frames, dtype=float)
    positions = np.zeros((n_frames, n_atoms, 3), dtype=float)
    positions[:, 0, 0] = 0.05 * np.sin(2 * np.pi * time / n_frames)
    positions[:, 1, 1] = 0.10 * np.cos(4 * np.pi * time / n_frames)
    positions[:, 2, 2] = 0.15 * np.sin(6 * np.pi * time / n_frames)

    qpoints = np.array(
        [
            [2 * np.pi, 0.0, 0.0],
            [0.0, 2 * np.pi, 0.0],
        ]
    )
    weights = np.array(
        [
            [1.0, 2.0, 3.0],
            [0.5, 1.5, 2.5],
        ]
    )

    result = compute_dsf(
        positions,
        qpoints,
        dt,
        experiment="custom",
        coherent_weights=weights,
        num_blocks=1,
    )

    expected = []
    for q_index, q_vec in enumerate(qpoints):
        phase = np.exp(1j * np.tensordot(positions, q_vec, axes=([2], [0])))
        rho = phase @ weights[q_index]
        spectrum = np.abs(fft(rho, norm="forward")) ** 2
        expected.append(spectrum[: n_frames // 2].real / (n_frames * n_atoms))

    np.testing.assert_allclose(result.coherent, np.asarray(expected), rtol=1e-14, atol=1e-14)


def test_dsf_cpu_backend_matches_default_backend():
    positions = np.zeros((16, 2, 3), dtype=float)
    positions[:, 1, 0] = np.linspace(0.0, 0.5, 16)
    qpoints = np.array([[2 * np.pi, 0.0, 0.0]])
    weights = np.array([1.0, 2.0])

    default = compute_dsf(
        positions,
        qpoints,
        dt=1e-12,
        experiment="custom",
        coherent_weights=weights,
        num_blocks=2,
    )
    cpu = compute_dsf(
        positions,
        qpoints,
        dt=1e-12,
        experiment="custom",
        coherent_weights=weights,
        num_blocks=2,
        backend="cpu",
    )

    np.testing.assert_allclose(cpu.coherent, default.coherent)
    np.testing.assert_allclose(cpu.coherent_sem, default.coherent_sem)
    assert cpu.metadata["backend"] == "cpu"


def test_xray_dsf_uses_q_dependent_form_factors_by_default():
    n_frames = 16
    positions = np.zeros((n_frames, 1, 3), dtype=float)
    qpoints = np.array([[0.0, 0.0, 0.0], [8.0, 0.0, 0.0]])

    result = compute_dsf(
        positions,
        qpoints,
        dt=1e-12,
        atom_types=["C"],
        experiment="xray",
    )
    weights = xray_form_factor_weights(["C"], qpoints)
    expected_ratio = (weights[1, 0] / weights[0, 0]) ** 2

    assert result.coherent[1, 0] < result.coherent[0, 0]
    np.testing.assert_allclose(
        result.coherent[1, 0] / result.coherent[0, 0],
        expected_ratio,
        rtol=1e-14,
    )
    assert result.metadata["xray_weights"] == "Cromer-Mann form factors with atomic-number fallback"


def test_dsf_reports_block_standard_error():
    block_size = 8
    n_frames = 2 * block_size
    positions = np.zeros((n_frames, 1, 3), dtype=float)
    positions[block_size:, 0, 0] = np.arange(block_size, dtype=float) / block_size

    result = compute_dsf(
        positions,
        np.array([[2 * np.pi, 0.0, 0.0]]),
        dt=1e-12,
        experiment="custom",
        num_blocks=2,
    )

    assert result.coherent_sem.shape == result.coherent.shape
    assert result.incoherent_sem.shape == result.incoherent.shape
    assert result.total_sem.shape == result.total.shape
    assert result.elastic_coherent_sem.shape == result.elastic_coherent.shape
    assert result.coherent_sem[0, 0] > 0
    assert result.coherent_sem[0, 1] > 0
    np.testing.assert_allclose(result.total_sem, result.coherent_sem)


def test_direct_dsf_matches_correlation_function_reference():
    n_frames = 24
    n_atoms = 2
    dt = 5e-13
    time = np.arange(n_frames, dtype=float)
    positions = np.zeros((n_frames, n_atoms, 3), dtype=float)
    positions[:, 0, 0] = 0.08 * np.sin(2 * np.pi * time / 12)
    positions[:, 0, 1] = 0.03 * np.cos(2 * np.pi * time / 8)
    positions[:, 1, 0] = 0.25 + 0.05 * np.cos(2 * np.pi * time / 6)
    positions[:, 1, 2] = 0.04 * np.sin(2 * np.pi * time / 4)
    qpoints = np.array([[2 * np.pi, 0.0, 0.0], [0.0, 2 * np.pi, 2 * np.pi]])
    weights = np.array([[1.0, 2.0], [0.75, 1.25]])

    direct = compute_dsf(
        positions,
        qpoints,
        dt,
        experiment="custom",
        coherent_weights=weights,
        num_blocks=3,
    )
    reference = compute_coherent_dsf_via_correlation(
        positions,
        qpoints,
        dt,
        coherent_weights=weights,
        num_blocks=3,
    )

    np.testing.assert_allclose(direct.coherent, reference.coherent, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(direct.coherent_sem, reference.coherent_sem, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(direct.elastic_coherent, reference.elastic_coherent)


def test_partial_dsf_reconstructs_weighted_coherent_total():
    n_frames = 24
    dt = 1e-12
    time = np.arange(n_frames, dtype=float)
    positions = np.zeros((n_frames, 4, 3), dtype=float)
    positions[:, 0, 0] = 0.03 * np.sin(2 * np.pi * time / 12)
    positions[:, 1, 0] = 0.50 + 0.04 * np.cos(2 * np.pi * time / 8)
    positions[:, 2, 1] = 0.25 + 0.02 * np.sin(2 * np.pi * time / 6)
    positions[:, 3, 2] = 0.75 + 0.05 * np.cos(2 * np.pi * time / 4)
    atom_types = np.array(["Si", "Si", "O", "O"])
    qpoints = np.array([[2 * np.pi, 0.0, 0.0], [0.0, 2 * np.pi, 2 * np.pi]])
    species_weights = {"Si": 2.0, "O": 0.5}
    atom_weights = np.array([species_weights[str(atom)] for atom in atom_types])

    total = compute_dsf(
        positions,
        qpoints,
        dt,
        experiment="custom",
        coherent_weights=atom_weights,
        num_blocks=3,
    )
    partial = compute_partial_dsf(
        positions,
        qpoints,
        dt,
        atom_types=atom_types,
        species_weights=species_weights,
        num_blocks=3,
    )

    assert partial.species == ("Si", "O")
    assert partial.partial.shape == (2, 2, 2, 4)
    np.testing.assert_allclose(partial.weighted_total, total.coherent, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(partial.weighted_total_sem, total.coherent_sem, rtol=1e-13, atol=1e-15)


def test_partial_dsf_cross_terms_are_hermitian():
    positions = np.zeros((8, 2, 3), dtype=float)
    positions[:, 1, 0] = np.linspace(0.0, 0.5, 8)
    result = compute_partial_dsf(
        positions,
        np.array([[2 * np.pi, 0.0, 0.0]]),
        dt=1e-12,
        atom_types=["A", "B"],
    )

    np.testing.assert_allclose(result.partial[0, 0, 1], np.conjugate(result.partial[0, 1, 0]))
