import numpy as np
import pytest

from pySED.eels import (
    build_eels_decomposition_map_from_mode_spectra,
    build_eels_map,
    build_eels_map_from_mode_spectra,
    compute_mode_visibility,
    compute_mode_visibility_decomposition,
    compute_mode_visibility_for_q_path,
)
from pySED.probe_weights import electron_form_factor_weights


def test_extended_zone_interference_selects_different_modes():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    masses = np.ones(2)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])

    eig = np.zeros((1, 2, 2, 3), dtype=complex)
    eig[0, 0, :, 0] = [1.0, 1.0] / np.sqrt(2.0)   # acoustic-like
    eig[0, 1, :, 0] = [1.0, -1.0] / np.sqrt(2.0)  # optic-like

    result = compute_mode_visibility(qpoints, g_vectors, primitive, basis_positions, masses, eig)

    acoustic_g1, optic_g1 = result.visibility[0, 0]
    acoustic_g2, optic_g2 = result.visibility[0, 1]
    assert optic_g1 > acoustic_g1
    assert acoustic_g2 > optic_g2


def test_visibility_decomposition_reports_atom_interference_and_direction_selection():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    masses = np.ones(2)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0]])

    eig = np.zeros((1, 2, 2, 3), dtype=complex)
    eig[0, 0, :, 0] = [1.0, 1.0] / np.sqrt(2.0)
    eig[0, 1, :, 0] = [1.0, -1.0] / np.sqrt(2.0)

    result = compute_mode_visibility_decomposition(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
    )

    assert result.atom_interference[0, 0, 0] < 0.0
    assert result.atom_interference[0, 0, 1] > 0.0
    np.testing.assert_allclose(result.direction_visibility[0, 0, :, 1:], 0.0, atol=1e-14)
    np.testing.assert_allclose(np.abs(result.amplitude) ** 2, result.visibility)


def test_eels_decomposition_map_preserves_atom_interference_identity():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    masses = np.ones(2)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0]])

    eig = np.zeros((1, 2, 2, 3), dtype=complex)
    eig[0, 0, :, 0] = [1.0, 1.0] / np.sqrt(2.0)
    eig[0, 1, :, 0] = [1.0, -1.0] / np.sqrt(2.0)

    decomposition = compute_mode_visibility_decomposition(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
    )
    energy = np.array([1.0, 2.0, 3.0])
    mode_spectra = np.zeros((3, 1, 2), dtype=float)
    mode_spectra[:, 0, 0] = [1.0, 2.0, 3.0]
    mode_spectra[:, 0, 1] = [4.0, 5.0, 6.0]

    diagnostic = build_eels_decomposition_map_from_mode_spectra(
        decomposition,
        mode_spectra,
        energy,
    )

    assert diagnostic.atom_intensity.shape == (1, 1, 2, 3)
    assert diagnostic.direction_intensity.shape == (1, 1, 3, 3)
    np.testing.assert_allclose(
        diagnostic.intensity,
        np.sum(diagnostic.atom_intensity, axis=2) + diagnostic.atom_interference_intensity,
    )


def test_build_eels_map_from_visibility():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0]])
    masses = np.ones(1)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0]])
    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0

    visibility = compute_mode_visibility(qpoints, g_vectors, primitive, basis_positions, masses, eig)
    energy = np.linspace(0, 10, 101)
    eels_map = build_eels_map(visibility, np.array([[5.0]]), energy, broadening=0.2)

    assert eels_map.intensity.shape == (1, 1, 101)
    assert energy[np.argmax(eels_map.intensity[0, 0])] == 5.0


def test_build_eels_map_from_branch_resolved_mode_spectra():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0]])
    masses = np.ones(1)
    qpoints = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    eig = np.zeros((2, 2, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0
    eig[:, 1, 0, 1] = 1.0

    visibility = compute_mode_visibility(qpoints, g_vectors, primitive, basis_positions, masses, eig)
    energy = np.array([1.0, 2.0, 3.0])
    mode_spectra = np.zeros((3, 2, 2), dtype=float)
    mode_spectra[:, 0, 0] = [1.0, 2.0, 3.0]
    mode_spectra[:, 0, 1] = [4.0, 5.0, 6.0]
    mode_spectra[:, 1, 0] = [7.0, 8.0, 9.0]
    mode_spectra[:, 1, 1] = [10.0, 11.0, 12.0]

    eels_map = build_eels_map_from_mode_spectra(visibility, mode_spectra, energy)
    expected = np.einsum("qgm,eqm->qge", visibility.visibility, mode_spectra)

    assert eels_map.intensity.shape == (2, 2, 3)
    np.testing.assert_allclose(eels_map.intensity, expected)


def test_pairwise_eels_q_path_does_not_form_q_by_g_cartesian_product():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0]])
    masses = np.ones(1)
    qpoints = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    eig = np.zeros((2, 1, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0

    product = compute_mode_visibility(qpoints, g_vectors, primitive, basis_positions, masses, eig)
    path = compute_mode_visibility_for_q_path(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
    )

    assert product.visibility.shape == (2, 2, 1)
    assert path.visibility.shape == (2, 1, 1)
    np.testing.assert_allclose(path.Q_reduced[:, 0, :], qpoints + g_vectors)
    assert path.metadata["pairwise_g_vectors"] is True


def test_pairwise_eels_q_path_validates_shape():
    with pytest.raises(ValueError, match="pairwise"):
        compute_mode_visibility_for_q_path(
            np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
            np.array([[1.0, 0.0, 0.0]]),
            np.eye(3),
            np.array([[0.0, 0.0, 0.0]]),
            np.ones(1),
            np.ones((2, 1, 1, 3), dtype=complex),
        )


def test_mott_bethe_electron_form_factor_weights_eels_visibility():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0]])
    masses = np.ones(1)
    qpoints = np.array([[0.5, 0.0, 0.0]])
    g_vectors = np.array([[0.0, 0.0, 0.0]])
    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0

    unit = compute_mode_visibility(qpoints, g_vectors, primitive, basis_positions, masses, eig)
    mott_bethe = compute_mode_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
        atom_types=["C"],
        electron_form_factor_model="mott-bethe",
    )
    expected_weight = electron_form_factor_weights(["C"], mott_bethe.Q_cartesian.reshape(-1, 3))[0, 0]

    np.testing.assert_allclose(mott_bethe.visibility, unit.visibility * expected_weight ** 2)
    assert mott_bethe.metadata["electron_form_factor_model"] == "mott-bethe"


def test_mott_bethe_eels_visibility_requires_atom_types():
    with pytest.raises(ValueError, match="atom_types"):
        compute_mode_visibility(
            np.array([[0.5, 0.0, 0.0]]),
            np.array([[0.0, 0.0, 0.0]]),
            np.eye(3),
            np.array([[0.0, 0.0, 0.0]]),
            np.ones(1),
            np.ones((1, 1, 1, 3), dtype=complex),
            electron_form_factor_model="mott-bethe",
        )


def test_explicit_electron_form_factors_override_model():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0]])
    masses = np.ones(1)
    qpoints = np.array([[0.5, 0.0, 0.0]])
    g_vectors = np.array([[0.0, 0.0, 0.0]])
    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0

    result = compute_mode_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
        electron_form_factors=np.array([2.0]),
        atom_types=["C"],
        electron_form_factor_model="mott-bethe",
    )
    unit = compute_mode_visibility(qpoints, g_vectors, primitive, basis_positions, masses, eig)

    np.testing.assert_allclose(result.visibility, unit.visibility * 4.0)
    assert result.metadata["electron_form_factor_model"] == "user supplied"


def test_build_eels_map_from_mode_spectra_validates_shape():
    visibility = compute_mode_visibility(
        np.array([[0.0, 0.0, 0.0]]),
        np.array([[1.0, 0.0, 0.0]]),
        np.eye(3),
        np.array([[0.0, 0.0, 0.0]]),
        np.ones(1),
        np.ones((1, 1, 1, 3), dtype=complex),
    )

    with pytest.raises(ValueError, match="mode_spectra"):
        build_eels_map_from_mode_spectra(visibility, np.ones((1, 1)), np.ones(1))
