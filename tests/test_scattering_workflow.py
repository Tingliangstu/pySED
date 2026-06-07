import numpy as np
import pytest

from pySED.compare_experiment import ScatteringMap
from pySED.eigen_sed import EigenvectorSet
from pySED.scattering_workflow import (
    compute_current_correlation_workflow,
    compute_dsf_workflow,
    compute_eels_workflow,
    compute_eels_workflow_from_angle_axis,
    compute_eels_workflow_from_q_path,
    compute_one_phonon_workflow,
    compute_one_phonon_workflow_from_q_path,
    compute_partial_dsf_workflow,
    current_correlation_to_scattering_map,
    dsf_to_scattering_map,
    eels_map_to_scattering_map,
    one_phonon_map_to_scattering_map,
    partial_dsf_to_scattering_map,
)
from pySED.electron_kinematics import electron_wavenumber_per_angstrom
from pySED.visibility_diagnostics import diagnose_mode_visibility


def _single_mode_inputs(qpoint=(0.0, 0.0, 0.0)):
    n_frames = 32
    dt = 1e-12
    target_bin = 4
    time = np.arange(n_frames) * dt
    freq_hz = target_bin / (n_frames * dt)

    velocities = np.zeros((n_frames, 1, 3), dtype=float)
    velocities[:, 0, 0] = np.cos(2 * np.pi * freq_hz * time)
    eig = np.zeros((1, 1, 1, 3), dtype=complex)
    eig[0, 0, 0, 0] = 1.0
    eigset = EigenvectorSet(
        qpoints_reduced=np.array([qpoint], dtype=float),
        frequencies=np.array([[freq_hz / 1e12]], dtype=float),
        eigenvectors=eig,
    )

    return {
        "velocities": velocities,
        "unitcell_vectors": np.zeros((1, 3), dtype=float),
        "basis_index": np.array([1]),
        "masses": np.array([1.0]),
        "primitive_cell": np.eye(3),
        "supercell_cell": np.eye(3),
        "eigenvectors": eigset,
        "basis_positions_reduced": np.array([[0.0, 0.0, 0.0]]),
        "g_vectors_reduced": np.array([[1.0, 0.0, 0.0]]),
        "dt": dt,
    }


def test_compute_eels_workflow_builds_scattering_map_and_comparison():
    inputs = _single_mode_inputs()
    preliminary = compute_eels_workflow(**inputs)
    experimental = ScatteringMap(
        preliminary.scattering_map.q_axis,
        preliminary.scattering_map.energy_axis,
        preliminary.scattering_map.intensity.copy(),
    )

    result = compute_eels_workflow(**inputs, experimental_map=experimental)

    assert result.eigen_sed.sed.shape[1:] == (1, 1)
    assert result.eels_map.intensity.shape[:2] == (1, 1)
    assert result.scattering_map.intensity.shape == (1, result.eigen_sed.frequencies_thz.size)
    assert result.comparison.rmse == pytest.approx(0.0)
    assert hasattr(result.visibility, "atom_visibility")
    diagnostics = diagnose_mode_visibility(result.visibility, q_advice=result.q_advice)
    assert diagnostics.dominant_reason.shape == (1, 1, 1)


def test_eels_workflow_rejects_non_commensurate_qpoints():
    inputs = _single_mode_inputs(qpoint=(0.25, 0.0, 0.0))

    with pytest.raises(ValueError, match="commensurate"):
        compute_eels_workflow(**inputs)


def test_eels_map_to_scattering_map_validates_g_index():
    inputs = _single_mode_inputs()
    result = compute_eels_workflow(**inputs)

    with pytest.raises(ValueError, match="g_index"):
        eels_map_to_scattering_map(result.eels_map, g_index=3)


def test_eels_workflow_can_output_mev_axis():
    inputs = _single_mode_inputs()
    result = compute_eels_workflow(**inputs, output_energy_unit="meV", sed_backend="cpu")

    np.testing.assert_allclose(
        result.scattering_map.energy_axis,
        result.eigen_sed.frequencies_thz * 4.135667696,
    )
    assert result.scattering_map.metadata["energy_unit"] == "meV"
    assert result.eigen_sed.metadata["backend"] == "cpu"


def test_eels_workflow_accepts_mott_bethe_electron_form_factor_model():
    inputs = _single_mode_inputs()
    result = compute_eels_workflow(
        **inputs,
        atom_types=["C"],
        electron_form_factor_model="mott-bethe",
    )

    assert result.visibility.metadata["electron_form_factor_model"] == "mott-bethe"
    assert result.eels_map.intensity.shape[:2] == (1, 1)
    assert np.all(np.isfinite(result.eels_map.intensity))


def test_eels_workflow_accepts_pairwise_extended_zone_q_path():
    n_frames = 32
    dt = 1e-12
    velocities = np.zeros((n_frames, 4, 3), dtype=float)
    velocities[:, :, 0] = 1.0
    eig = np.zeros((2, 1, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0
    eigset = EigenvectorSet(
        qpoints_reduced=np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
        frequencies=np.array([[0.0], [0.0]], dtype=float),
        eigenvectors=eig,
    )

    result = compute_eels_workflow(
        velocities=velocities,
        unitcell_vectors=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        ),
        basis_index=np.array([1, 1, 1, 1]),
        masses=np.array([1.0]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        eigenvectors=eigset,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        g_vectors_reduced=np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
        dt=dt,
        pairwise_g_vectors=True,
    )

    assert result.eels_map.intensity.shape[:2] == (2, 1)
    np.testing.assert_allclose(result.visibility.Q_reduced[:, 0, :], [[1.0, 0.0, 0.0], [2.25, 0.0, 0.0]])
    assert result.visibility.metadata["pairwise_g_vectors"] is True


def test_eels_workflow_from_experimental_q_path_builds_q_plan():
    n_frames = 32
    dt = 1e-12
    velocities = np.zeros((n_frames, 4, 3), dtype=float)
    velocities[:, :, 0] = 1.0
    eig = np.zeros((2, 1, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0
    eigset = EigenvectorSet(
        qpoints_reduced=np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
        frequencies=np.array([[0.0], [0.0]], dtype=float),
        eigenvectors=eig,
    )

    result = compute_eels_workflow_from_q_path(
        velocities=velocities,
        unitcell_vectors=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        ),
        basis_index=np.array([1, 1, 1, 1]),
        masses=np.array([1.0]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        eigenvectors=eigset,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        experimental_qpoints=np.array([[1.0 * 2.0 * np.pi, 0.0, 0.0], [2.25 * 2.0 * np.pi, 0.0, 0.0]]),
        dt=dt,
        qpoint_coordinates="cartesian",
    )

    np.testing.assert_allclose(result.q_plan.qpoints_reduced, [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    np.testing.assert_array_equal(result.q_plan.g_vectors_reduced, [[1, 0, 0], [2, 0, 0]])
    np.testing.assert_allclose(result.visibility.Q_reduced[:, 0, :], result.q_plan.Q_reduced)
    assert result.visibility.metadata["pairwise_g_vectors"] is True
    assert result.eigen_sed.metadata["q_plan"] == "experimental extended-zone Q path folded to commensurate q"


def test_eels_workflow_from_angle_axis_builds_q_plan():
    n_frames = 32
    dt = 1e-12
    beam_energy = 200000.0
    velocities = np.zeros((n_frames, 4, 3), dtype=float)
    velocities[:, :, 0] = 1.0
    eig = np.zeros((2, 1, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0
    eigset = EigenvectorSet(
        qpoints_reduced=np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
        frequencies=np.array([[0.0], [0.0]], dtype=float),
        eigenvectors=eig,
    )
    k0 = electron_wavenumber_per_angstrom(beam_energy)
    angle_axis_mrad = np.array([0.0, (0.25 * 2.0 * np.pi / k0) * 1000.0])

    result = compute_eels_workflow_from_angle_axis(
        velocities=velocities,
        unitcell_vectors=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        ),
        basis_index=np.array([1, 1, 1, 1]),
        masses=np.array([1.0]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        eigenvectors=eigset,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        angle_axis=angle_axis_mrad,
        beam_energy_ev=beam_energy,
        dt=dt,
        q_policy="strict",
    )

    np.testing.assert_allclose(result.q_plan.qpoints_reduced, [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    np.testing.assert_array_equal(result.q_plan.g_vectors_reduced, [[0, 0, 0], [0, 0, 0]])
    assert result.eigen_sed.metadata["q_plan_source"] == "electron scattering angle axis"
    assert result.visibility.metadata["beam_energy_ev"] == beam_energy


def test_eels_workflow_from_experimental_q_path_rejects_non_commensurate_strict():
    inputs = _single_mode_inputs()

    with pytest.raises(ValueError, match="non-commensurate"):
        compute_eels_workflow_from_q_path(
            velocities=inputs["velocities"],
            unitcell_vectors=inputs["unitcell_vectors"],
            basis_index=inputs["basis_index"],
            masses=inputs["masses"],
            primitive_cell=np.eye(3),
            supercell_cell=np.diag([4.0, 1.0, 1.0]),
            eigenvectors=inputs["eigenvectors"],
            basis_positions_reduced=inputs["basis_positions_reduced"],
            experimental_qpoints=np.array([[0.26, 0.0, 0.0]]),
            dt=inputs["dt"],
            qpoint_coordinates="reduced",
            q_policy="strict",
        )


def test_compute_one_phonon_workflow_builds_scattering_map_and_comparison():
    inputs = _single_mode_inputs()
    preliminary = compute_one_phonon_workflow(**inputs, experiment="custom", sed_backend="cpu")
    experimental = ScatteringMap(
        preliminary.scattering_map.q_axis,
        preliminary.scattering_map.energy_axis,
        preliminary.scattering_map.intensity.copy(),
    )

    result = compute_one_phonon_workflow(
        **inputs,
        experiment="custom",
        sed_backend="cpu",
        experimental_map=experimental,
    )

    assert result.eigen_sed.sed.shape[1:] == (1, 1)
    assert result.one_phonon_map.intensity.shape[:2] == (1, 1)
    assert result.scattering_map.intensity.shape == (1, result.eigen_sed.frequencies_thz.size)
    assert result.comparison.rmse == pytest.approx(0.0)
    assert result.eigen_sed.metadata["backend"] == "cpu"


def test_one_phonon_workflow_from_experimental_q_path_builds_q_plan():
    n_frames = 32
    dt = 1e-12
    velocities = np.zeros((n_frames, 4, 3), dtype=float)
    velocities[:, :, 0] = 1.0
    eig = np.zeros((2, 1, 1, 3), dtype=complex)
    eig[:, 0, 0, 0] = 1.0
    eigset = EigenvectorSet(
        qpoints_reduced=np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
        frequencies=np.array([[1.0], [1.0]], dtype=float),
        eigenvectors=eig,
    )

    result = compute_one_phonon_workflow_from_q_path(
        velocities=velocities,
        unitcell_vectors=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        ),
        basis_index=np.array([1, 1, 1, 1]),
        masses=np.array([1.0]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        eigenvectors=eigset,
        basis_positions_reduced=np.array([[0.0, 0.0, 0.0]]),
        experimental_qpoints=np.array([[1.0 * 2.0 * np.pi, 0.0, 0.0], [2.25 * 2.0 * np.pi, 0.0, 0.0]]),
        dt=dt,
        qpoint_coordinates="cartesian",
        experiment="custom",
    )

    np.testing.assert_allclose(result.q_plan.qpoints_reduced, [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    np.testing.assert_array_equal(result.q_plan.g_vectors_reduced, [[1, 0, 0], [2, 0, 0]])
    np.testing.assert_allclose(result.visibility.Q_reduced[:, 0, :], result.q_plan.Q_reduced)
    assert result.one_phonon_map.intensity.shape[:2] == (2, 1)
    assert result.visibility.metadata["pairwise_g_vectors"] is True


def test_one_phonon_workflow_from_experimental_q_path_rejects_non_commensurate_strict():
    inputs = _single_mode_inputs()

    with pytest.raises(ValueError, match="non-commensurate"):
        compute_one_phonon_workflow_from_q_path(
            velocities=inputs["velocities"],
            unitcell_vectors=inputs["unitcell_vectors"],
            basis_index=inputs["basis_index"],
            masses=inputs["masses"],
            primitive_cell=np.eye(3),
            supercell_cell=np.diag([4.0, 1.0, 1.0]),
            eigenvectors=inputs["eigenvectors"],
            basis_positions_reduced=inputs["basis_positions_reduced"],
            experimental_qpoints=np.array([[0.26, 0.0, 0.0]]),
            dt=inputs["dt"],
            qpoint_coordinates="reduced",
            q_policy="strict",
            experiment="custom",
        )


def test_one_phonon_map_to_scattering_map_validates_g_index():
    inputs = _single_mode_inputs()
    result = compute_one_phonon_workflow(**inputs, experiment="custom")

    with pytest.raises(ValueError, match="g_index"):
        one_phonon_map_to_scattering_map(result.one_phonon_map, g_index=2)


def test_dsf_workflow_rejects_non_commensurate_qpoints_in_strict_mode():
    positions = np.zeros((16, 1, 3), dtype=float)

    with pytest.raises(ValueError, match="commensurate"):
        compute_dsf_workflow(
            positions,
            qpoints=np.array([[0.20, 0.0, 0.0]]),
            primitive_cell=np.eye(3),
            supercell_cell=np.diag([4.0, 1.0, 1.0]),
            dt=1e-12,
            qpoint_coordinates="reduced",
            q_policy="strict",
            experiment="custom",
        )


def test_dsf_workflow_maps_nearest_commensurate_qpoints():
    positions = np.zeros((16, 1, 3), dtype=float)
    result = compute_dsf_workflow(
        positions,
        qpoints=np.array([[0.26, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="nearest",
        max_error_reduced=0.02,
        experiment="custom",
    )

    np.testing.assert_allclose(result.qpoints_reduced, [[0.25, 0.0, 0.0]])
    np.testing.assert_allclose(result.qpoints_cartesian, [[0.5 * np.pi, 0.0, 0.0]])
    assert result.dsf.metadata["q_policy"] == "nearest"


def test_dsf_workflow_preserves_extended_zone_cartesian_q_mapping():
    positions = np.zeros((16, 1, 3), dtype=float)
    cartesian_q = np.array([[1.26 * 2.0 * np.pi, 0.0, 0.0]])

    result = compute_dsf_workflow(
        positions,
        qpoints=cartesian_q,
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        qpoint_coordinates="cartesian",
        q_policy="nearest",
        experiment="custom",
    )

    np.testing.assert_allclose(result.qpoints_reduced, [[1.25, 0.0, 0.0]])
    np.testing.assert_allclose(result.qpoints_cartesian, [[2.5 * np.pi, 0.0, 0.0]])
    np.testing.assert_allclose(result.dsf.qpoints_cartesian, result.qpoints_cartesian)


def test_dsf_workflow_builds_scattering_map_and_comparison():
    positions = np.zeros((16, 1, 3), dtype=float)
    preliminary = compute_dsf_workflow(
        positions,
        qpoints=np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="strict",
        experiment="custom",
    )
    experimental = ScatteringMap(
        preliminary.scattering_map.q_axis,
        preliminary.scattering_map.energy_axis,
        preliminary.scattering_map.intensity.copy(),
    )

    result = compute_dsf_workflow(
        positions,
        qpoints=np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="strict",
        experiment="custom",
        experimental_map=experimental,
    )

    assert result.scattering_map.intensity.shape == result.dsf.total.shape
    assert result.scattering_map.metadata["component"] == "total"
    assert result.comparison.rmse == pytest.approx(0.0)


def test_dsf_to_scattering_map_can_select_component_and_energy_unit():
    positions = np.zeros((16, 1, 3), dtype=float)
    result = compute_dsf_workflow(
        positions,
        qpoints=np.array([[0.0, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.eye(3),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="strict",
        experiment="custom",
    )

    scattering_map = dsf_to_scattering_map(
        result.dsf,
        component="coherent",
        output_energy_unit="meV",
    )

    np.testing.assert_allclose(scattering_map.intensity, result.dsf.coherent)
    np.testing.assert_allclose(scattering_map.energy_axis, result.dsf.frequencies_thz * 4.135667696)
    assert scattering_map.metadata["component"] == "coherent"
    assert scattering_map.metadata["energy_unit"] == "meV"


def test_dsf_to_scattering_map_rejects_unknown_component():
    positions = np.zeros((16, 1, 3), dtype=float)
    result = compute_dsf_workflow(
        positions,
        qpoints=np.array([[0.0, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.eye(3),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="strict",
        experiment="custom",
    )

    with pytest.raises(ValueError, match="component"):
        dsf_to_scattering_map(result.dsf, component="bad")


def test_partial_dsf_workflow_builds_weighted_total_and_species_pair_maps():
    n_frames = 16
    dt = 1e-12
    positions = np.zeros((n_frames, 2, 3), dtype=float)
    positions[:, 1, 0] = np.linspace(0.0, 0.5, n_frames)

    result = compute_partial_dsf_workflow(
        positions,
        qpoints=np.array([[0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=dt,
        atom_types=["A", "B"],
        qpoint_coordinates="reduced",
        q_policy="strict",
        species_weights={"A": 1.0, "B": 2.0},
        num_blocks=2,
    )

    assert result.scattering_map.intensity.shape == result.partial_dsf.weighted_total.shape
    assert result.scattering_map.metadata["source"] == "pySED PartialDSFResult"
    assert result.partial_dsf.metadata["q_policy"] == "strict"

    pair_map = partial_dsf_to_scattering_map(
        result.partial_dsf,
        component="partial",
        species_pair=("A", "B"),
        complex_part="abs",
        output_energy_unit="meV",
    )

    assert pair_map.metadata["species_pair"] == ("A", "B")
    assert pair_map.metadata["complex_part"] == "abs"
    np.testing.assert_allclose(pair_map.energy_axis, result.partial_dsf.frequencies_thz * 4.135667696)
    assert np.all(pair_map.intensity >= 0.0)


def test_partial_dsf_to_scattering_map_rejects_missing_species_pair():
    positions = np.zeros((8, 2, 3), dtype=float)
    result = compute_partial_dsf_workflow(
        positions,
        qpoints=np.array([[0.0, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.eye(3),
        dt=1e-12,
        atom_types=["A", "B"],
        qpoint_coordinates="reduced",
        q_policy="strict",
    )

    with pytest.raises(ValueError, match="species_pair"):
        partial_dsf_to_scattering_map(result.partial_dsf, component="partial")


def test_current_correlation_workflow_builds_scattering_map_and_sem():
    n_frames = 16
    dt = 1e-12
    positions = np.zeros((n_frames, 1, 3), dtype=float)
    velocities = np.zeros_like(positions)
    velocities[:, 0, 0] = 1.0

    preliminary = compute_current_correlation_workflow(
        positions,
        velocities,
        qpoints=np.array([[0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=dt,
        qpoint_coordinates="reduced",
        q_policy="strict",
        num_blocks=2,
    )
    experimental = ScatteringMap(
        preliminary.scattering_map.q_axis,
        preliminary.scattering_map.energy_axis,
        preliminary.scattering_map.intensity.copy(),
    )
    result = compute_current_correlation_workflow(
        positions,
        velocities,
        qpoints=np.array([[0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=dt,
        qpoint_coordinates="reduced",
        q_policy="strict",
        num_blocks=2,
        experimental_map=experimental,
    )

    assert result.scattering_map.intensity.shape == result.current.longitudinal.shape
    assert result.current.longitudinal_sem.shape == result.current.longitudinal.shape
    assert result.current.metadata["q_policy"] == "strict"
    assert result.comparison.rmse == pytest.approx(0.0)

    transverse_map = current_correlation_to_scattering_map(
        result.current,
        component="transverse",
        output_energy_unit="meV",
    )

    np.testing.assert_allclose(transverse_map.intensity, result.current.transverse)
    np.testing.assert_allclose(transverse_map.energy_axis, result.current.frequencies_thz * 4.135667696)
    assert transverse_map.metadata["component"] == "transverse"
