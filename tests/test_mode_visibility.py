import numpy as np
import pytest

from pySED.mode_visibility import (
    build_one_phonon_decomposition_map_from_mode_spectra,
    build_one_phonon_map,
    build_one_phonon_map_from_mode_spectra,
    compute_one_phonon_visibility,
    compute_one_phonon_visibility_for_q_path,
)


def test_custom_one_phonon_visibility_reports_interference_and_direction_selection():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    masses = np.ones(2)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0]])
    eig = np.zeros((1, 2, 2, 3), dtype=complex)
    eig[0, 0, :, 0] = [1.0, 1.0] / np.sqrt(2.0)
    eig[0, 1, :, 0] = [1.0, -1.0] / np.sqrt(2.0)

    result = compute_one_phonon_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
        experiment="custom",
    )

    acoustic, optic = result.visibility[0, 0]
    assert optic > acoustic
    assert result.atom_interference[0, 0, 0] < 0.0
    assert result.atom_interference[0, 0, 1] > 0.0
    np.testing.assert_allclose(result.direction_visibility[0, 0, :, 1:], 0.0, atol=1e-14)
    np.testing.assert_allclose(np.abs(result.amplitude) ** 2, result.visibility)


def test_xray_one_phonon_visibility_uses_q_dependent_form_factors():
    primitive = np.eye(3)
    qpoints = np.array([[0.5, 0.0, 0.0]])
    g_vectors = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0

    result = compute_one_phonon_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        atom_types=["C"],
        experiment="xray",
    )

    assert result.scattering_weights.shape == (1, 2, 1)
    assert result.scattering_weights[0, 1, 0] < result.scattering_weights[0, 0, 0]
    assert result.metadata["weight_model"] == "Cromer-Mann X-ray form factors"


def test_neutron_one_phonon_visibility_requires_atom_types():
    eig = np.ones((1, 1, 1, 3), dtype=complex)

    with pytest.raises(ValueError, match="atom_types"):
        compute_one_phonon_visibility(
            np.array([[0.0, 0.0, 0.0]]),
            np.array([[1.0, 0.0, 0.0]]),
            np.eye(3),
            basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
            masses=np.ones(1),
            eigenvectors=eig,
            experiment="neutron",
        )


def test_frequency_power_applies_one_over_omega_visibility_factor():
    eig = np.zeros((1, 2, 1, 3), dtype=complex)
    eig[0, :, 0, 0] = 1.0

    result = compute_one_phonon_visibility(
        np.array([[1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0]]),
        np.eye(3),
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        experiment="custom",
        mode_frequencies=np.array([[1.0, 2.0]]),
        frequency_power=-1,
    )

    np.testing.assert_allclose(result.visibility[0, 0, 0] / result.visibility[0, 0, 1], 2.0)


def test_build_one_phonon_map_from_visibility():
    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0
    visibility = compute_one_phonon_visibility(
        np.array([[1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0]]),
        np.eye(3),
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        experiment="custom",
    )
    energy = np.linspace(0.0, 10.0, 101)

    scattering_map = build_one_phonon_map(
        visibility,
        mode_frequencies=np.array([[5.0]]),
        energy_axis=energy,
        broadening=0.2,
        broadening_kind="gaussian",
    )

    assert scattering_map.intensity.shape == (1, 1, 101)
    assert energy[np.argmax(scattering_map.intensity[0, 0])] == 5.0


def test_build_one_phonon_map_from_mode_spectra():
    eig = np.zeros((2, 2, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0
    eig[:, 1, 0, 1] = 1.0
    visibility = compute_one_phonon_visibility(
        np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        np.eye(3),
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        experiment="custom",
    )
    energy = np.array([1.0, 2.0, 3.0])
    mode_spectra = np.zeros((3, 2, 2), dtype=float)
    mode_spectra[:, 0, 0] = [1.0, 2.0, 3.0]
    mode_spectra[:, 0, 1] = [4.0, 5.0, 6.0]
    mode_spectra[:, 1, 0] = [7.0, 8.0, 9.0]
    mode_spectra[:, 1, 1] = [10.0, 11.0, 12.0]

    one_phonon_map = build_one_phonon_map_from_mode_spectra(visibility, mode_spectra, energy)
    expected = np.einsum("qgm,eqm->qge", visibility.visibility, mode_spectra)

    assert one_phonon_map.intensity.shape == (2, 2, 3)
    np.testing.assert_allclose(one_phonon_map.intensity, expected)


def test_pairwise_one_phonon_q_path_does_not_form_q_by_g_cartesian_product():
    primitive = np.eye(3)
    qpoints = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    eig = np.zeros((2, 1, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0

    product = compute_one_phonon_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        experiment="custom",
    )
    path = compute_one_phonon_visibility_for_q_path(
        qpoints,
        g_vectors,
        primitive,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        experiment="custom",
    )

    assert product.visibility.shape == (2, 2, 1)
    assert path.visibility.shape == (2, 1, 1)
    np.testing.assert_allclose(path.Q_reduced[:, 0, :], qpoints + g_vectors)
    assert path.metadata["pairwise_g_vectors"] is True


def test_pairwise_one_phonon_q_path_validates_shape():
    with pytest.raises(ValueError, match="pairwise"):
        compute_one_phonon_visibility_for_q_path(
            np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
            np.array([[1.0, 0.0, 0.0]]),
            np.eye(3),
            basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
            masses=np.ones(1),
            eigenvectors=np.ones((2, 1, 1, 3), dtype=complex),
            experiment="custom",
        )


def test_one_phonon_decomposition_map_preserves_atom_interference_identity():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    masses = np.ones(2)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0]])
    eig = np.zeros((1, 2, 2, 3), dtype=complex)
    eig[0, 0, :, 0] = [1.0, 1.0] / np.sqrt(2.0)
    eig[0, 1, :, 0] = [1.0, -1.0] / np.sqrt(2.0)

    visibility = compute_one_phonon_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
        experiment="custom",
    )
    energy = np.array([1.0, 2.0, 3.0])
    mode_spectra = np.zeros((3, 1, 2), dtype=float)
    mode_spectra[:, 0, 0] = [1.0, 2.0, 3.0]
    mode_spectra[:, 0, 1] = [4.0, 5.0, 6.0]

    diagnostic = build_one_phonon_decomposition_map_from_mode_spectra(
        visibility,
        mode_spectra,
        energy,
    )

    assert diagnostic.atom_intensity.shape == (1, 1, 2, 3)
    assert diagnostic.direction_intensity.shape == (1, 1, 3, 3)
    np.testing.assert_allclose(
        diagnostic.intensity,
        np.sum(diagnostic.atom_intensity, axis=2) + diagnostic.atom_interference_intensity,
    )


def test_build_one_phonon_map_from_mode_spectra_validates_shape():
    eig = np.ones((1, 1, 1, 3), dtype=complex)
    visibility = compute_one_phonon_visibility(
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[1.0, 0.0, 0.0]]),
        np.eye(3),
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        masses=np.ones(1),
        eigenvectors=eig,
        experiment="custom",
    )

    with pytest.raises(ValueError, match="mode_spectra"):
        build_one_phonon_map_from_mode_spectra(visibility, np.ones((1, 1)), np.ones(1))
