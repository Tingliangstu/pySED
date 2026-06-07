import numpy as np

import pySED


def test_scattering_objects_are_available_from_package_top_level():
    advisor = pySED.QAdvisor(np.eye(3), np.diag([2.0, 1.0, 1.0]))
    assert advisor.is_commensurate([0.5, 0.0, 0.0])

    assert pySED.build_eels_decomposition_map_from_mode_spectra is not None
    assert pySED.build_one_phonon_decomposition_map_from_mode_spectra is not None
    assert pySED.compute_coherent_dsf_via_correlation is not None
    assert pySED.compute_current_correlation_workflow is not None
    assert pySED.compute_dsf is not None
    assert pySED.compute_eels_workflow is not None
    assert pySED.compute_eels_workflow_from_angle_axis is not None
    assert pySED.compute_eels_workflow_from_q_path is not None
    assert pySED.dsf_to_scattering_map is not None
    assert pySED.compute_mode_visibility_for_q_path is not None
    assert pySED.compute_one_phonon_visibility_for_q_path is not None
    assert pySED.compute_one_phonon_workflow_from_q_path is not None
    assert pySED.compute_partial_dsf is not None
    assert pySED.compute_partial_dsf_workflow is not None
    assert pySED.angular_resolution_to_q_sigma is not None
    assert pySED.diagnose_mode_visibility is not None
    assert pySED.eels_qpoints_from_angle_axis is not None
    assert pySED.scattering_angles_to_qpoints is not None
    assert pySED.write_map_comparison_report is not None
    assert pySED.write_scattering_workflow_bundle is not None
    assert pySED.workflow_export_summary is not None
    assert pySED.run_scattering_from_input is not None
    assert pySED.comparison_summary is not None
    assert pySED.comparison_linecut_records is not None
    assert pySED.estimate_selected_q_efficiency is not None
    assert pySED.qpoints_from_mesh is not None
    assert pySED.qpoints_from_path_axis is not None
    assert pySED.qpoints_from_path_segments is not None
    assert pySED.split_extended_zone_qpoints is not None
    assert pySED.ScatteringMap is not None
    assert pySED.QPointMesh is not None
    assert pySED.QPathSegments is not None
    assert pySED.SelectedQEfficiencyReport is not None
    assert pySED.CurrentCorrelationResult is not None
    assert pySED.PartialDSFResult is not None
    assert pySED.current_correlation_to_scattering_map is not None
    assert pySED.partial_dsf_to_scattering_map is not None
