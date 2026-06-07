"""Validation helpers for scattering implementations."""

from dataclasses import dataclass
import importlib
import importlib.util
import json
import platform
import sys

import numpy as np

from pySED.dsf import compute_coherent_dsf_via_correlation, compute_dsf
from pySED.version import __version__


@dataclass
class DSFValidationCase:
    """Small deterministic DSF validation problem."""

    positions: np.ndarray
    qpoints_cartesian: np.ndarray
    coherent_weights: np.ndarray
    dt: float
    num_blocks: int = 4
    description: str = "synthetic deterministic trajectory"


@dataclass
class SpectrumComparison:
    """Metrics for comparing a pySED spectrum with an external reference."""

    passed: bool
    candidate_shape: tuple
    reference_shape: tuple
    normalization: str
    scale_factor: float
    rmse: float
    mae: float
    max_abs: float
    normalized_rmse: float
    correlation: float
    peak_frequency_mae: float = None
    reason: str = None

    def to_dict(self):
        return {
            "passed": bool(self.passed),
            "candidate_shape": tuple(self.candidate_shape),
            "reference_shape": tuple(self.reference_shape),
            "normalization": self.normalization,
            "scale_factor": float(self.scale_factor),
            "rmse": float(self.rmse),
            "mae": float(self.mae),
            "max_abs": float(self.max_abs),
            "normalized_rmse": float(self.normalized_rmse),
            "correlation": float(self.correlation),
            "peak_frequency_mae": None
            if self.peak_frequency_mae is None
            else float(self.peak_frequency_mae),
            "reason": self.reason,
        }


def synthetic_dsf_validation_case(num_frames=64):
    """Return a deterministic trajectory and weighted Q points."""

    time = np.arange(num_frames, dtype=float)
    positions = np.zeros((num_frames, 3, 3), dtype=float)
    positions[:, 0, 0] = 0.05 * np.sin(2 * np.pi * time / 16)
    positions[:, 1, 0] = 0.50 + 0.08 * np.cos(2 * np.pi * time / 8)
    positions[:, 1, 1] = 0.04 * np.sin(2 * np.pi * time / 16)
    positions[:, 2, 2] = 0.25 + 0.03 * np.cos(2 * np.pi * time / 4)
    qpoints = np.array([[2 * np.pi, 0.0, 0.0], [0.0, 2 * np.pi, 2 * np.pi]])
    weights = np.array([[1.0, 2.0, 1.5], [0.5, 1.5, 2.5]])
    return DSFValidationCase(
        positions=positions,
        qpoints_cartesian=qpoints,
        coherent_weights=weights,
        dt=1e-13,
        num_blocks=4,
    )


def compute_pysed_validation_pair(case=None):
    """Compute pySED direct and explicit-correlation DSF for one case."""

    if case is None:
        case = synthetic_dsf_validation_case()
    direct = compute_dsf(
        case.positions,
        case.qpoints_cartesian,
        case.dt,
        experiment="custom",
        coherent_weights=case.coherent_weights,
        num_blocks=case.num_blocks,
    )
    reference = compute_coherent_dsf_via_correlation(
        case.positions,
        case.qpoints_cartesian,
        case.dt,
        coherent_weights=case.coherent_weights,
        num_blocks=case.num_blocks,
    )
    return direct, reference


def compare_arrays(candidate, reference, atol=1e-12, rtol=1e-12):
    """Return scalar comparison metrics between two arrays."""

    candidate = np.asarray(candidate)
    reference = np.asarray(reference)
    if candidate.shape != reference.shape:
        return {
            "passed": False,
            "reason": "shape mismatch",
            "candidate_shape": tuple(candidate.shape),
            "reference_shape": tuple(reference.shape),
        }
    delta = candidate - reference
    max_abs = float(np.nanmax(np.abs(delta))) if delta.size else 0.0
    reference_norm = float(np.nanmax(np.abs(reference))) if reference.size else 0.0
    scale = atol + rtol * reference_norm
    return {
        "passed": bool(max_abs <= scale),
        "max_abs": max_abs,
        "reference_max_abs": reference_norm,
        "atol": float(atol),
        "rtol": float(rtol),
        "threshold": float(scale),
    }


def _interpolate_last_axis(values, source_frequencies, target_frequencies):
    values = np.asarray(values, dtype=float)
    source = np.asarray(source_frequencies, dtype=float)
    target = np.asarray(target_frequencies, dtype=float)
    if source.ndim != 1 or target.ndim != 1:
        raise ValueError("frequency arrays must be one-dimensional")
    if values.shape[-1] != source.size:
        raise ValueError("source_frequencies length must match spectrum last axis")

    flat = values.reshape(-1, values.shape[-1])
    output = np.zeros((flat.shape[0], target.size), dtype=float)
    order = np.argsort(source)
    source_sorted = source[order]
    for row_index, row in enumerate(flat):
        output[row_index] = np.interp(target, source_sorted, row[order])
    return output.reshape(values.shape[:-1] + (target.size,))


def _least_squares_scale(candidate, reference):
    numerator = float(np.sum(candidate * reference))
    denominator = float(np.sum(candidate * candidate))
    if denominator == 0.0:
        return 1.0
    return numerator / denominator


def _peak_frequency_mae(candidate, reference, candidate_frequencies, reference_frequencies):
    if candidate_frequencies is None or reference_frequencies is None:
        return None
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    candidate = candidate.reshape(-1, candidate.shape[-1])
    reference = reference.reshape(-1, reference.shape[-1])
    candidate_freq = np.asarray(candidate_frequencies, dtype=float)
    reference_freq = np.asarray(reference_frequencies, dtype=float)
    candidate_peaks = candidate_freq[np.argmax(candidate, axis=1)]
    reference_peaks = reference_freq[np.argmax(reference, axis=1)]
    return float(np.mean(np.abs(candidate_peaks - reference_peaks)))


def compare_reference_spectrum(
    candidate,
    reference,
    candidate_frequencies=None,
    reference_frequencies=None,
    normalization="least_squares",
    rmse_tolerance=None,
):
    """Compare a spectrum against an external reference spectrum.

    The last axis is treated as frequency.  If both frequency grids are supplied
    and differ, the reference is interpolated onto the candidate grid before
    scalar metrics are computed.
    """

    candidate_arr = np.asarray(candidate, dtype=float)
    reference_arr = np.asarray(reference, dtype=float)
    original_reference = reference_arr
    original_reference_frequencies = reference_frequencies

    if candidate_frequencies is not None and reference_frequencies is not None:
        candidate_freq = np.asarray(candidate_frequencies, dtype=float)
        reference_freq = np.asarray(reference_frequencies, dtype=float)
        if candidate_arr.shape[:-1] != reference_arr.shape[:-1]:
            return SpectrumComparison(
                passed=False,
                candidate_shape=candidate_arr.shape,
                reference_shape=reference_arr.shape,
                normalization=str(normalization),
                scale_factor=1.0,
                rmse=np.nan,
                mae=np.nan,
                max_abs=np.nan,
                normalized_rmse=np.nan,
                correlation=np.nan,
                reason="non-frequency dimensions do not match",
            )
        if candidate_freq.shape != reference_freq.shape or not np.allclose(candidate_freq, reference_freq):
            reference_arr = _interpolate_last_axis(reference_arr, reference_freq, candidate_freq)
            reference_frequencies = candidate_freq
    elif candidate_arr.shape != reference_arr.shape:
        return SpectrumComparison(
            passed=False,
            candidate_shape=candidate_arr.shape,
            reference_shape=reference_arr.shape,
            normalization=str(normalization),
            scale_factor=1.0,
            rmse=np.nan,
            mae=np.nan,
            max_abs=np.nan,
            normalized_rmse=np.nan,
            correlation=np.nan,
            reason="shape mismatch and no frequency grids supplied for interpolation",
        )

    if candidate_arr.shape != reference_arr.shape:
        return SpectrumComparison(
            passed=False,
            candidate_shape=candidate_arr.shape,
            reference_shape=reference_arr.shape,
            normalization=str(normalization),
            scale_factor=1.0,
            rmse=np.nan,
            mae=np.nan,
            max_abs=np.nan,
            normalized_rmse=np.nan,
            correlation=np.nan,
            reason="shape mismatch after interpolation",
        )

    if not np.all(np.isfinite(candidate_arr)) or not np.all(np.isfinite(reference_arr)):
        return SpectrumComparison(
            passed=False,
            candidate_shape=candidate_arr.shape,
            reference_shape=reference_arr.shape,
            normalization=str(normalization),
            scale_factor=1.0,
            rmse=np.nan,
            mae=np.nan,
            max_abs=np.nan,
            normalized_rmse=np.nan,
            correlation=np.nan,
            reason="candidate or reference spectrum contains non-finite values",
        )

    normalization_key = str(normalization).lower()
    if normalization_key in ("none", "raw"):
        scale_factor = 1.0
        candidate_scaled = candidate_arr
        reference_scaled = reference_arr
    elif normalization_key in ("least_squares", "least-squares", "scale"):
        scale_factor = _least_squares_scale(candidate_arr, reference_arr)
        candidate_scaled = scale_factor * candidate_arr
        reference_scaled = reference_arr
    elif normalization_key == "max":
        candidate_max = float(np.nanmax(np.abs(candidate_arr))) if candidate_arr.size else 0.0
        reference_max = float(np.nanmax(np.abs(reference_arr))) if reference_arr.size else 0.0
        scale_factor = reference_max / candidate_max if candidate_max != 0.0 else 1.0
        candidate_scaled = scale_factor * candidate_arr
        reference_scaled = reference_arr
    else:
        raise ValueError("normalization must be 'none', 'max', or 'least_squares'")

    delta = candidate_scaled - reference_scaled
    rmse = float(np.sqrt(np.nanmean(delta ** 2))) if delta.size else 0.0
    mae = float(np.nanmean(np.abs(delta))) if delta.size else 0.0
    max_abs = float(np.nanmax(np.abs(delta))) if delta.size else 0.0
    reference_norm = float(np.sqrt(np.nanmean(reference_scaled ** 2))) if reference_scaled.size else 0.0
    normalized_rmse = rmse / reference_norm if reference_norm != 0.0 else rmse

    cand_flat = candidate_scaled.ravel()
    ref_flat = reference_scaled.ravel()
    if cand_flat.size < 2 or np.std(cand_flat) == 0.0 or np.std(ref_flat) == 0.0:
        correlation = 1.0 if np.allclose(cand_flat, ref_flat) else 0.0
    else:
        correlation = float(np.corrcoef(cand_flat, ref_flat)[0, 1])

    peak_mae = _peak_frequency_mae(
        candidate_scaled,
        original_reference,
        candidate_frequencies,
        original_reference_frequencies,
    )
    passed = True if rmse_tolerance is None else normalized_rmse <= float(rmse_tolerance)
    return SpectrumComparison(
        passed=bool(passed),
        candidate_shape=candidate_arr.shape,
        reference_shape=reference_arr.shape,
        normalization=normalization_key,
        scale_factor=scale_factor,
        rmse=rmse,
        mae=mae,
        max_abs=max_abs,
        normalized_rmse=normalized_rmse,
        correlation=correlation,
        peak_frequency_mae=peak_mae,
    )


def run_internal_dsf_reference_validation(case=None, tolerance=1e-12):
    """Validate pySED's direct estimator against its explicit correlation path."""

    if case is None:
        case = synthetic_dsf_validation_case()
    direct, reference = compute_pysed_validation_pair(case)
    coherent = compare_arrays(direct.coherent, reference.coherent, atol=tolerance, rtol=tolerance)
    coherent_sem = compare_arrays(
        direct.coherent_sem,
        reference.coherent_sem,
        atol=tolerance,
        rtol=tolerance,
    )
    elastic = compare_arrays(
        direct.elastic_coherent,
        reference.elastic_coherent,
        atol=tolerance,
        rtol=tolerance,
    )
    return {
        "name": "pysed_direct_vs_explicit_correlation",
        "passed": bool(coherent["passed"] and coherent_sem["passed"] and elastic["passed"]),
        "case": case.description,
        "num_frames": int(case.positions.shape[0]),
        "num_atoms": int(case.positions.shape[1]),
        "num_qpoints": int(case.qpoints_cartesian.shape[0]),
        "num_blocks": int(case.num_blocks),
        "dt_seconds": float(case.dt),
        "checks": {
            "coherent": coherent,
            "coherent_sem": coherent_sem,
            "elastic_coherent": elastic,
        },
    }


def package_status(module_name):
    """Return installation metadata for an optional validation dependency."""

    spec = importlib.util.find_spec(module_name)
    if spec is None:
        return {"module": module_name, "installed": False}
    try:
        module = importlib.import_module(module_name)
        version = getattr(module, "__version__", "unknown")
        path = getattr(module, "__file__", "unknown")
    except Exception as exc:
        return {
            "module": module_name,
            "installed": True,
            "importable": False,
            "error": "%s: %s" % (type(exc).__name__, exc),
        }
    return {
        "module": module_name,
        "installed": True,
        "importable": True,
        "version": str(version),
        "path": str(path),
    }


def optional_scattering_package_report():
    """Report optional package availability for external DSF validation."""

    packages = [
        package_status("dynasor"),
        package_status("psf"),
        package_status("pynamic_structure_factor"),
        package_status("pynamic"),
    ]
    notes = {
        "dynasor": (
            "Use dynasor.compute_dynamic_structure_factors on the same trajectory "
            "and Q points, then compare weighted Sqw on a matched frequency grid."
        ),
        "pynamic": (
            "Use pynamic-structure-factor's direct selected-Q workflow with the "
            "same Q points, weights, block/window convention, and normalization."
        ),
    }
    return {"packages": packages, "notes": notes}


def validation_environment():
    """Return reproducibility metadata for a validation report."""

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pysed": __version__,
    }


def build_dsf_validation_report(case=None, tolerance=1e-12):
    """Build a structured DSF validation report."""

    if case is None:
        case = synthetic_dsf_validation_case()
    internal = run_internal_dsf_reference_validation(case, tolerance=tolerance)
    external = optional_scattering_package_report()
    return {
        "environment": validation_environment(),
        "internal": internal,
        "external": external,
        "passed": bool(internal["passed"]),
    }


def write_validation_report(report, path):
    """Write a JSON validation report."""

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)


def export_dsf_validation_case(path, case=None, include_pysed_reference=True):
    """Export a deterministic validation case as an ``.npz`` file."""

    if case is None:
        case = synthetic_dsf_validation_case()
    payload = {
        "positions": case.positions,
        "qpoints_cartesian": case.qpoints_cartesian,
        "coherent_weights": case.coherent_weights,
        "dt": np.asarray(case.dt),
        "num_blocks": np.asarray(case.num_blocks),
        "description": np.asarray(case.description),
    }
    if include_pysed_reference:
        direct, reference = compute_pysed_validation_pair(case)
        payload.update(
            {
                "pysed_direct_coherent": direct.coherent,
                "pysed_reference_coherent": reference.coherent,
                "pysed_frequencies_thz": direct.frequencies_thz,
            }
        )
    np.savez(path, **payload)


def load_dsf_validation_case(path):
    """Load a validation case exported by :func:`export_dsf_validation_case`."""

    data = np.load(path, allow_pickle=False)
    return DSFValidationCase(
        positions=data["positions"],
        qpoints_cartesian=data["qpoints_cartesian"],
        coherent_weights=data["coherent_weights"],
        dt=float(data["dt"]),
        num_blocks=int(data["num_blocks"]),
        description=str(data["description"]),
    )


def compare_external_dsf_reference(
    reference_path,
    reference_key="coherent",
    reference_frequency_key=None,
    case=None,
    component="coherent",
    normalization="least_squares",
    rmse_tolerance=None,
):
    """Compare pySED DSF with an external dynasor/pynamic-style ``.npz`` output."""

    if case is None:
        case = synthetic_dsf_validation_case()
    direct, _ = compute_pysed_validation_pair(case)
    if not hasattr(direct, component):
        raise ValueError("candidate component %s is not available" % component)
    candidate = getattr(direct, component)
    candidate_frequencies = direct.frequencies_thz

    data = np.load(reference_path, allow_pickle=False)
    if reference_key not in data.files:
        raise KeyError("reference_key %s not found in %s" % (reference_key, reference_path))
    reference = data[reference_key]

    if reference_frequency_key is None:
        for key in ("frequencies_thz", "frequency_thz", "frequencies", "frequency", "omega"):
            if key in data.files:
                reference_frequency_key = key
                break
    elif reference_frequency_key not in data.files:
        raise KeyError(
            "reference_frequency_key %s not found in %s" % (reference_frequency_key, reference_path)
        )
    reference_frequencies = data[reference_frequency_key] if reference_frequency_key is not None else None

    comparison = compare_reference_spectrum(
        candidate,
        reference,
        candidate_frequencies=candidate_frequencies,
        reference_frequencies=reference_frequencies,
        normalization=normalization,
        rmse_tolerance=rmse_tolerance,
    )
    payload = comparison.to_dict()
    payload.update(
        {
            "name": "pysed_vs_external_dsf_reference",
            "reference_path": str(reference_path),
            "reference_key": reference_key,
            "reference_frequency_key": reference_frequency_key,
            "candidate_component": component,
            "candidate_frequency_unit": "THz",
        }
    )
    return payload
