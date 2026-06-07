"""SED calculation code: pySED."""

from pySED.version import __version__  # noqa F401
from pySED.compare_experiment import MapComparison, ScatteringMap
from pySED.compare_experiment import comparison_linecut_records
from pySED.compare_experiment import comparison_peak_records, comparison_residual_map
from pySED.compare_experiment import comparison_summary, write_map_comparison_report
from pySED.dsf import CurrentCorrelationResult, DSFResult, PartialDSFResult
from pySED.dsf import compute_coherent_dsf_via_correlation
from pySED.dsf import compute_current_correlations, compute_dsf, compute_partial_dsf
from pySED.electron_kinematics import angular_resolution_to_q_sigma
from pySED.electron_kinematics import eels_qpoints_from_angle_axis
from pySED.electron_kinematics import electron_wavelength_angstrom
from pySED.electron_kinematics import scattering_angles_to_qpoints
from pySED.eels import EELSDecompositionMap, EELSMap, EELSVisibility
from pySED.eels import build_eels_decomposition_map_from_mode_spectra
from pySED.eels import compute_mode_visibility
from pySED.eels import compute_mode_visibility_for_q_path
from pySED.eigen_sed import EigenSEDResult, EigenvectorSet, compute_eigen_sed
from pySED.mode_visibility import OnePhononDecompositionMap, OnePhononMap
from pySED.mode_visibility import OnePhononVisibility
from pySED.mode_visibility import build_one_phonon_decomposition_map_from_mode_spectra
from pySED.mode_visibility import compute_one_phonon_visibility
from pySED.mode_visibility import compute_one_phonon_visibility_for_q_path
from pySED.q_advisor import ExtendedZoneQPathPlan, QAdvisor, QAdvisorReport
from pySED.q_advisor import QPointAdvice, QPointMesh, QPathSegments, SelectedQEfficiencyReport
from pySED.q_advisor import estimate_selected_q_efficiency
from pySED.q_advisor import qpoints_from_mesh, qpoints_from_path_axis
from pySED.q_advisor import qpoints_from_path_segments, split_extended_zone_qpoints
from pySED.scattering_export import workflow_export_summary
from pySED.scattering_export import write_scattering_workflow_bundle
from pySED.scattering_input import ScatteringInputError
from pySED.scattering_input import build_q_advisor_from_input
from pySED.scattering_input import build_scattering_qpoints_from_input
from pySED.scattering_input import run_q_advisor_from_input, run_scattering_from_input
from pySED.scattering_workflow import CurrentCorrelationWorkflowResult, PartialDSFWorkflowResult
from pySED.scattering_workflow import compute_current_correlation_workflow
from pySED.scattering_workflow import compute_dsf_workflow, compute_eels_workflow
from pySED.scattering_workflow import compute_eels_workflow_from_angle_axis
from pySED.scattering_workflow import compute_eels_workflow_from_q_path
from pySED.scattering_workflow import compute_one_phonon_workflow
from pySED.scattering_workflow import compute_one_phonon_workflow_from_q_path
from pySED.scattering_workflow import compute_partial_dsf_workflow
from pySED.scattering_workflow import current_correlation_to_scattering_map
from pySED.scattering_workflow import dsf_to_scattering_map
from pySED.scattering_workflow import partial_dsf_to_scattering_map
from pySED.visibility_diagnostics import ModeVisibilityDiagnostics
from pySED.visibility_diagnostics import diagnose_mode_visibility

__all__ = [
    "__version__",
    "CurrentCorrelationResult",
    "CurrentCorrelationWorkflowResult",
    "DSFResult",
    "EELSDecompositionMap",
    "EELSMap",
    "EELSVisibility",
    "EigenSEDResult",
    "EigenvectorSet",
    "ExtendedZoneQPathPlan",
    "MapComparison",
    "OnePhononDecompositionMap",
    "OnePhononMap",
    "OnePhononVisibility",
    "PartialDSFResult",
    "PartialDSFWorkflowResult",
    "QAdvisor",
    "QAdvisorReport",
    "QPointAdvice",
    "QPointMesh",
    "QPathSegments",
    "ScatteringInputError",
    "SelectedQEfficiencyReport",
    "estimate_selected_q_efficiency",
    "qpoints_from_mesh",
    "qpoints_from_path_axis",
    "qpoints_from_path_segments",
    "ScatteringMap",
    "ModeVisibilityDiagnostics",
    "build_eels_decomposition_map_from_mode_spectra",
    "build_q_advisor_from_input",
    "build_one_phonon_decomposition_map_from_mode_spectra",
    "build_scattering_qpoints_from_input",
    "angular_resolution_to_q_sigma",
    "compute_coherent_dsf_via_correlation",
    "compute_current_correlations",
    "compute_current_correlation_workflow",
    "comparison_linecut_records",
    "comparison_peak_records",
    "comparison_residual_map",
    "comparison_summary",
    "compute_dsf",
    "compute_dsf_workflow",
    "dsf_to_scattering_map",
    "eels_qpoints_from_angle_axis",
    "electron_wavelength_angstrom",
    "compute_eels_workflow",
    "compute_eels_workflow_from_angle_axis",
    "compute_eels_workflow_from_q_path",
    "compute_eigen_sed",
    "compute_mode_visibility",
    "compute_mode_visibility_for_q_path",
    "compute_one_phonon_visibility",
    "compute_one_phonon_visibility_for_q_path",
    "compute_one_phonon_workflow",
    "compute_one_phonon_workflow_from_q_path",
    "compute_partial_dsf",
    "compute_partial_dsf_workflow",
    "current_correlation_to_scattering_map",
    "diagnose_mode_visibility",
    "partial_dsf_to_scattering_map",
    "run_q_advisor_from_input",
    "run_scattering_from_input",
    "scattering_angles_to_qpoints",
    "split_extended_zone_qpoints",
    "write_map_comparison_report",
    "workflow_export_summary",
    "write_scattering_workflow_bundle",
]
