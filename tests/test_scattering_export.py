import csv
import json

import numpy as np

from pySED.compare_experiment import ScatteringMap
from pySED.eigen_sed import EigenvectorSet
from pySED.scattering_export import workflow_export_summary, write_scattering_workflow_bundle
from pySED.scattering_workflow import compute_current_correlation_workflow
from pySED.scattering_workflow import compute_dsf_workflow, compute_eels_workflow_from_q_path
from pySED.scattering_workflow import compute_partial_dsf_workflow


def _eels_q_path_result_with_comparison():
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

    preliminary = compute_eels_workflow_from_q_path(
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
    experimental = ScatteringMap(
        preliminary.scattering_map.q_axis,
        preliminary.scattering_map.energy_axis,
        preliminary.scattering_map.intensity.copy(),
    )
    return compute_eels_workflow_from_q_path(
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
        experimental_map=experimental,
        sigma_q=0.0,
        sigma_energy=0.0,
    )


def test_scattering_workflow_bundle_exports_eels_artifacts(tmp_path):
    result = _eels_q_path_result_with_comparison()

    written = write_scattering_workflow_bundle(
        result,
        tmp_path,
        prefix="eels",
        atom_labels=["C"],
        linecut_q_indices=[0],
        linecut_energy_indices=[0],
    )

    expected_keys = {
        "summary",
        "q_plan",
        "phonopy_qpoints",
        "visibility_diagnostics",
        "scattering_map",
        "comparison_summary",
        "comparison_residual",
        "comparison_peaks",
        "comparison_linecuts",
    }
    assert expected_keys <= set(written)

    with open(written["summary"], "r", encoding="utf-8") as handle:
        summary = json.load(handle)
    with open(written["visibility_diagnostics"], "r", encoding="utf-8") as handle:
        diagnostic_rows = list(csv.DictReader(handle))

    assert summary["workflow_type"] == "EELSWorkflowResult"
    assert summary["has_q_plan"] is True
    assert summary["comparison"]["rmse"] < 1e-30
    assert diagnostic_rows[0]["dominant_atom_label"] == "C"
    np.testing.assert_allclose(np.loadtxt(written["phonopy_qpoints"]), [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])


def test_scattering_workflow_bundle_exports_dsf_q_advice(tmp_path):
    positions = np.zeros((16, 1, 3), dtype=float)
    result = compute_dsf_workflow(
        positions,
        qpoints=np.array([[0.26, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="nearest",
        experiment="custom",
    )

    summary = workflow_export_summary(result)
    written = write_scattering_workflow_bundle(result, tmp_path, prefix="dsf")

    assert summary["workflow_type"] == "DSFWorkflowResult"
    assert summary["q_advice"]["num_commensurate"] == 0
    assert "q_advice" in written
    assert "visibility_diagnostics" not in written
    assert "scattering_map" in written


def test_scattering_workflow_bundle_exports_partial_and_current_metadata(tmp_path):
    positions = np.zeros((16, 2, 3), dtype=float)
    positions[:, 1, 0] = np.linspace(0.0, 0.5, 16)
    velocities = np.zeros_like(positions)
    velocities[:, :, 0] = 1.0

    partial = compute_partial_dsf_workflow(
        positions,
        qpoints=np.array([[0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        atom_types=["A", "B"],
        qpoint_coordinates="reduced",
        q_policy="strict",
    )
    current = compute_current_correlation_workflow(
        positions,
        velocities,
        qpoints=np.array([[0.25, 0.0, 0.0]]),
        primitive_cell=np.eye(3),
        supercell_cell=np.diag([4.0, 1.0, 1.0]),
        dt=1e-12,
        qpoint_coordinates="reduced",
        q_policy="strict",
    )

    partial_summary = workflow_export_summary(partial)
    current_summary = workflow_export_summary(current)
    partial_written = write_scattering_workflow_bundle(partial, tmp_path / "partial", prefix="partial")
    current_written = write_scattering_workflow_bundle(current, tmp_path / "current", prefix="current")

    assert partial_summary["workflow_type"] == "PartialDSFWorkflowResult"
    assert partial_summary["partial_dsf_metadata"]["q_policy"] == "strict"
    assert current_summary["workflow_type"] == "CurrentCorrelationWorkflowResult"
    assert current_summary["current_metadata"]["uncertainty"] == "standard error of the block-averaged spectra"
    assert "q_advice" in partial_written
    assert "scattering_map" in partial_written
    assert "q_advice" in current_written
    assert "scattering_map" in current_written
