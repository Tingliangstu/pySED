import json

import numpy as np

from pySED.validation import (
    build_dsf_validation_report,
    compare_arrays,
    compare_external_dsf_reference,
    compare_reference_spectrum,
    compute_pysed_validation_pair,
    export_dsf_validation_case,
    load_dsf_validation_case,
    package_status,
    run_internal_dsf_reference_validation,
    synthetic_dsf_validation_case,
    write_validation_report,
)


def test_internal_dsf_validation_report_passes():
    report = run_internal_dsf_reference_validation()

    assert report["passed"]
    assert report["checks"]["coherent"]["max_abs"] < 1e-12
    assert report["num_qpoints"] == 2


def test_build_dsf_validation_report_contains_environment_and_external_status():
    report = build_dsf_validation_report()

    assert report["passed"]
    assert "python" in report["environment"]
    assert "packages" in report["external"]
    assert any(package["module"] == "dynasor" for package in report["external"]["packages"])


def test_write_validation_report_json(tmp_path):
    path = tmp_path / "validation.json"
    report = build_dsf_validation_report()
    write_validation_report(report, path)

    with open(path, "r", encoding="utf-8") as handle:
        loaded = json.load(handle)
    assert loaded["passed"] is True
    assert loaded["internal"]["name"] == "pysed_direct_vs_explicit_correlation"


def test_export_and_load_dsf_validation_case(tmp_path):
    path = tmp_path / "case.npz"
    case = synthetic_dsf_validation_case(num_frames=16)

    export_dsf_validation_case(path, case=case)
    loaded = load_dsf_validation_case(path)
    data = np.load(path)

    np.testing.assert_allclose(loaded.positions, case.positions)
    np.testing.assert_allclose(loaded.qpoints_cartesian, case.qpoints_cartesian)
    assert loaded.num_blocks == case.num_blocks
    assert "pysed_direct_coherent" in data.files
    assert "pysed_reference_coherent" in data.files


def test_compare_arrays_reports_shape_mismatch():
    result = compare_arrays(np.zeros((2,)), np.zeros((3,)))

    assert not result["passed"]
    assert result["reason"] == "shape mismatch"


def test_compare_reference_spectrum_exact_match():
    spectrum = np.array([[0.0, 1.0, 0.25], [0.5, 0.0, 0.75]])

    comparison = compare_reference_spectrum(spectrum, spectrum, normalization="none")

    assert comparison.passed
    assert comparison.rmse == 0.0
    assert comparison.normalized_rmse == 0.0
    assert comparison.correlation == 1.0


def test_compare_reference_spectrum_least_squares_scale():
    reference = np.array([[0.0, 1.0, 2.0], [3.0, 0.0, 1.0]])
    candidate = 2.0 * reference

    comparison = compare_reference_spectrum(candidate, reference, normalization="least_squares")

    assert comparison.passed
    assert comparison.scale_factor == 0.5
    assert comparison.rmse == 0.0


def test_compare_reference_spectrum_interpolates_frequency_axis():
    candidate_frequencies = np.array([0.0, 1.0, 2.0])
    reference_frequencies = np.array([0.0, 2.0])
    candidate = np.array([[0.0, 1.0, 2.0]])
    reference = np.array([[0.0, 2.0]])

    comparison = compare_reference_spectrum(
        candidate,
        reference,
        candidate_frequencies=candidate_frequencies,
        reference_frequencies=reference_frequencies,
        normalization="none",
    )

    assert comparison.passed
    assert comparison.rmse == 0.0
    assert comparison.peak_frequency_mae == 0.0


def test_compare_reference_spectrum_shape_mismatch_without_frequencies():
    comparison = compare_reference_spectrum(np.zeros((2, 3)), np.zeros((2, 4)))

    assert not comparison.passed
    assert comparison.reason == "shape mismatch and no frequency grids supplied for interpolation"


def test_compare_reference_spectrum_non_finite_values_fail():
    comparison = compare_reference_spectrum(np.array([1.0, np.nan]), np.array([1.0, 0.0]))

    assert not comparison.passed
    assert comparison.reason == "candidate or reference spectrum contains non-finite values"


def test_compare_external_dsf_reference_npz_self_match(tmp_path):
    path = tmp_path / "external_reference.npz"
    case = synthetic_dsf_validation_case()
    direct, _ = compute_pysed_validation_pair(case)
    np.savez(path, coherent=direct.coherent, frequencies_thz=direct.frequencies_thz)

    comparison = compare_external_dsf_reference(
        path,
        reference_key="coherent",
        reference_frequency_key="frequencies_thz",
        case=case,
        normalization="none",
        rmse_tolerance=1e-12,
    )

    assert comparison["passed"]
    assert comparison["normalized_rmse"] == 0.0
    assert comparison["reference_frequency_key"] == "frequencies_thz"


def test_compare_external_dsf_reference_missing_frequency_key_raises(tmp_path):
    path = tmp_path / "external_reference.npz"
    np.savez(path, coherent=np.zeros((2, 32)))

    try:
        compare_external_dsf_reference(path, reference_frequency_key="missing")
    except KeyError as exc:
        assert "reference_frequency_key missing not found" in str(exc)
    else:
        raise AssertionError("missing reference frequency key should raise KeyError")


def test_package_status_for_missing_module():
    status = package_status("definitely_missing_pysed_validation_module")

    assert status["installed"] is False
