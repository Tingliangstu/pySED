import numpy as np
import pytest

from pySED.compare_experiment import ScatteringMap
from pySED.q_advisor import (
    QAdvisor,
    estimate_selected_q_efficiency,
    qpoints_from_mesh,
    qpoints_from_path_axis,
    qpoints_from_path_segments,
    split_extended_zone_qpoints,
    suggest_diagonal_supercell,
)


def test_commensurate_validation_and_nearest_point():
    primitive = np.eye(3)
    supercell = np.diag([4.0, 4.0, 1.0])
    advisor = QAdvisor(primitive, supercell)

    assert advisor.is_commensurate([0.25, 0.0, 0.0])
    assert not advisor.is_commensurate([0.20, 0.0, 0.0])

    advice = advisor.nearest_commensurate([0.26, 0.0, 0.0])
    np.testing.assert_allclose(advice.nearest_reduced, [0.25, 0.0, 0.0])
    assert advice.error_cartesian > 0


def test_commensurate_path_uses_supercell_grid():
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))
    qpoints, fractions = advisor.build_commensurate_path([0, 0, 0], [0.5, 0, 0])

    np.testing.assert_allclose(qpoints, [[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0]])
    np.testing.assert_allclose(fractions, [0.0, 0.5, 1.0])


def test_suggest_diagonal_supercell_from_target_qpoints():
    qpoints = np.array([[1 / 3, 0, 0], [0, 0.25, 0], [0, 0, 0.5]])
    np.testing.assert_array_equal(suggest_diagonal_supercell(qpoints), [3, 4, 2])


def test_cartesian_experimental_qpoints_are_mapped_to_commensurate_grid():
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 4.0, 1.0]))
    cartesian_q = np.array([[0.26 * 2 * np.pi, 0.0, 0.0]])

    advice = advisor.advise_cartesian_qpoints(cartesian_q)[0]

    np.testing.assert_allclose(advice.nearest_reduced, [0.25, 0.0, 0.0])
    np.testing.assert_allclose(advice.nearest_cartesian, [0.5 * np.pi, 0.0, 0.0])
    assert advice.error_cartesian > 0


def test_experimental_path_report_summarizes_errors_and_exports_csv(tmp_path):
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 4.0, 1.0]))
    qpoints = np.array([[0.25, 0.0, 0.0], [0.20, 0.0, 0.0]])

    report = advisor.advise_experimental_path(qpoints, coordinates="reduced")

    assert report.num_points == 2
    assert report.num_commensurate == 1
    assert not report.all_commensurate
    np.testing.assert_allclose(report.nearest_reduced[1], [0.25, 0.0, 0.0])
    np.testing.assert_array_equal(report.suggested_diagonal_supercell, [20, 1, 1])
    assert report.summary()["max_error_cartesian"] > 0.0
    assert "suggested_diagonal_supercell" in report.to_table()

    csv_path = tmp_path / "q_report.csv"
    report.write_csv(csv_path)
    text = csv_path.read_text(encoding="utf-8")
    assert "requested_reduced_1" in text
    assert "False" in text


def test_nearest_commensurate_can_preserve_extended_zone_image():
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))

    folded = advisor.nearest_commensurate([1.26, 0.0, 0.0])
    preserved = advisor.nearest_commensurate([1.26, 0.0, 0.0], preserve_image=True)

    np.testing.assert_allclose(folded.nearest_reduced, [0.25, 0.0, 0.0])
    np.testing.assert_allclose(preserved.nearest_reduced, [1.25, 0.0, 0.0])
    np.testing.assert_allclose(preserved.nearest_cartesian, [2.5 * np.pi, 0.0, 0.0])


def test_experimental_path_report_preserves_cartesian_extended_zone_image():
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))
    cartesian_q = np.array([[1.26 * 2.0 * np.pi, 0.0, 0.0]])

    report = advisor.advise_experimental_path(cartesian_q, coordinates="cartesian")

    np.testing.assert_allclose(report.requested_reduced, [[1.26, 0.0, 0.0]])
    np.testing.assert_allclose(report.nearest_reduced, [[1.25, 0.0, 0.0]])
    np.testing.assert_allclose(report.nearest_cartesian, [[2.5 * np.pi, 0.0, 0.0]])
    np.testing.assert_array_equal(report.suggested_diagonal_supercell, [50, 1, 1])


def test_split_extended_zone_qpoints_returns_folded_q_and_integer_g():
    q_folded, g_vectors = split_extended_zone_qpoints(
        np.array([[1.25, 0.0, 0.0], [0.75, 0.0, 0.0]])
    )

    np.testing.assert_allclose(q_folded, [[0.25, 0.0, 0.0], [-0.25, 0.0, 0.0]])
    np.testing.assert_array_equal(g_vectors, [[1, 0, 0], [1, 0, 0]])
    np.testing.assert_allclose(q_folded + g_vectors, [[1.25, 0.0, 0.0], [0.75, 0.0, 0.0]])


def test_extended_zone_q_path_plan_maps_experimental_q_to_pairwise_q_and_g(tmp_path):
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))

    plan = advisor.plan_extended_zone_q_path(
        np.array([[1.26 * 2.0 * np.pi, 0.0, 0.0]]),
        coordinates="cartesian",
        q_policy="nearest",
        max_error_reduced=0.02,
    )

    np.testing.assert_allclose(plan.Q_reduced, [[1.25, 0.0, 0.0]])
    np.testing.assert_allclose(plan.qpoints_reduced, [[0.25, 0.0, 0.0]])
    np.testing.assert_array_equal(plan.g_vectors_reduced, [[1, 0, 0]])
    np.testing.assert_allclose(plan.Q_reduced, plan.qpoints_reduced + plan.g_vectors_reduced)

    phonopy_path = tmp_path / "phonopy_qpoints.dat"
    plan.write_phonopy_qpoints(phonopy_path)
    np.testing.assert_allclose(np.loadtxt(phonopy_path).reshape(1, 3), [[0.25, 0.0, 0.0]])

    csv_path = tmp_path / "extended_q_plan.csv"
    plan.write_csv(csv_path)
    text = csv_path.read_text(encoding="utf-8")
    assert "folded_q_reduced_1" in text
    assert "G_reduced_1" in text


def test_qpoints_from_path_axis_interpolates_vector_endpoints():
    qpoints = qpoints_from_path_axis(
        q_axis=np.array([10.0, 20.0, 30.0]),
        start_qpoint=np.array([0.0, 0.0, 0.0]),
        end_qpoint=np.array([0.5, 0.0, 0.0]),
    )

    np.testing.assert_allclose(qpoints, [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]])


def test_qpoints_from_mesh_generates_selected_q_grid_and_reshape():
    mesh = qpoints_from_mesh([-0.25, 0.25, 3], 0.0, [1.0, 1.5, 2])

    assert mesh.mesh_shape == (3, 1, 2)
    np.testing.assert_allclose(mesh.axes[0], [-0.25, 0.0, 0.25])
    np.testing.assert_allclose(
        mesh.qpoints_reduced,
        [
            [-0.25, 0.0, 1.0],
            [-0.25, 0.0, 1.5],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.5],
            [0.25, 0.0, 1.0],
            [0.25, 0.0, 1.5],
        ],
    )

    values = np.arange(6)
    np.testing.assert_array_equal(mesh.reshape_values(values), values.reshape(3, 1, 2))


def test_qpoints_from_path_segments_matches_selected_path_convention():
    path = qpoints_from_path_segments(
        starts=np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
        ends=np.array([[0.5, 0.0, 0.0], [0.5, 0.5, 0.0]]),
        steps=np.array([2, 2]),
    )

    np.testing.assert_allclose(
        path.qpoints_reduced,
        [
            [0.0, 0.0, 0.0],
            [0.25, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.25, 0.0],
        ],
    )
    np.testing.assert_allclose(path.path_axis, [0.0, 0.25, 0.5, 0.75])
    np.testing.assert_array_equal(path.segment_start_indices, [0, 2])
    np.testing.assert_array_equal(path.segment_end_indices, [1, 3])


def test_selected_mesh_points_can_be_checked_by_q_advisor():
    mesh = qpoints_from_mesh([0.20, 0.30, 2], 0.0, 0.0)
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))

    report = advisor.advise_experimental_path(mesh.qpoints_reduced, coordinates="reduced")

    assert not report.all_commensurate
    np.testing.assert_allclose(report.nearest_reduced[:, 0], [0.25, 0.25])
    np.testing.assert_array_equal(report.suggested_diagonal_supercell, [10, 1, 1])


def test_selected_q_efficiency_report_estimates_full_grid_savings():
    report = estimate_selected_q_efficiency(
        num_selected_qpoints=8,
        num_full_qpoints=64,
        num_frames=16,
        num_atoms=10,
    )

    assert report.selected_fraction == 0.125
    assert report.qpoint_reduction_factor == 8.0
    assert report.selected_phase_evaluations == 1280.0
    assert report.full_phase_evaluations == 10240.0
    assert report.selected_fft_work_units == 8 * 16 * 4
    assert "q-point reduction factor" in report.to_table()
    assert report.summary()["num_atoms"] == 10


def test_q_advisor_selected_q_efficiency_uses_supercell_grid_size():
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 4.0, 2.0]))
    mesh = qpoints_from_mesh([0.0, 0.25, 2], 0.0, 0.0)

    report = advisor.estimate_selected_q_efficiency(
        mesh.qpoints_reduced,
        num_frames=32,
        num_atoms=5,
    )

    assert report.num_full_qpoints == 32
    assert report.num_selected_qpoints == 2
    assert report.qpoint_reduction_factor == 16.0
    assert report.full_phase_evaluations / report.selected_phase_evaluations == 16.0


def test_map_q_axis_advisor_uses_calibrated_image_axis():
    scattering_map = ScatteringMap(
        q_axis=np.array([0.0, 1.0, 2.0]),
        energy_axis=np.array([0.0, 10.0]),
        intensity=np.ones((3, 2)),
    )
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))

    report = advisor.advise_map_q_axis(
        scattering_map,
        start_qpoint=np.array([0.0, 0.0, 0.0]),
        end_qpoint=np.array([0.5, 0.0, 0.0]),
        coordinates="reduced",
    )
    plan = advisor.plan_map_q_axis(
        scattering_map,
        start_qpoint=np.array([1.0, 0.0, 0.0]),
        end_qpoint=np.array([1.5, 0.0, 0.0]),
        coordinates="reduced",
        q_policy="strict",
    )

    assert report.all_commensurate
    np.testing.assert_allclose(report.requested_reduced[:, 0], [0.0, 0.25, 0.5])
    np.testing.assert_allclose(plan.Q_reduced[:, 0], [1.0, 1.25, 1.5])
    np.testing.assert_allclose(plan.qpoints_reduced[:, 0], [0.0, 0.25, 0.5])
    np.testing.assert_array_equal(plan.g_vectors_reduced[:, 0], [1, 1, 1])


def test_extended_zone_q_path_plan_strict_rejects_non_commensurate_q():
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 1.0, 1.0]))

    with pytest.raises(ValueError, match="non-commensurate"):
        advisor.plan_extended_zone_q_path(
            np.array([[0.26, 0.0, 0.0]]),
            coordinates="reduced",
            q_policy="strict",
        )


def test_write_phonopy_qpoints_rejects_non_commensurate_points(tmp_path):
    advisor = QAdvisor(np.eye(3), np.diag([4.0, 4.0, 1.0]))
    qpoints = np.array([[0.25, 0.0, 0.0], [0.20, 0.0, 0.0]])

    assert not advisor.validate_qpoints(qpoints)
    with pytest.raises(ValueError, match="commensurate"):
        advisor.write_phonopy_qpoints(qpoints, tmp_path / "qpoints.dat")
