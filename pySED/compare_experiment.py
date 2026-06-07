"""Utilities for comparing simulated scattering maps with experiments."""

from dataclasses import dataclass
from pathlib import Path
import csv
import json

import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter, gaussian_filter1d

from pySED.quantum_correction import apply_quantum_correction


@dataclass
class ScatteringMap:
    q_axis: np.ndarray
    energy_axis: np.ndarray
    intensity: np.ndarray
    metadata: dict = None

    def __post_init__(self):
        self.q_axis = np.asarray(self.q_axis, dtype=float)
        self.energy_axis = np.asarray(self.energy_axis, dtype=float)
        self.intensity = np.asarray(self.intensity, dtype=float)
        if self.intensity.ndim != 2:
            raise ValueError("ScatteringMap intensity must have shape (num_q, num_energy)")
        if self.intensity.shape != (self.q_axis.size, self.energy_axis.size):
            raise ValueError("intensity shape must match q_axis and energy_axis lengths")
        if self.metadata is None:
            self.metadata = {}


@dataclass
class MapComparison:
    simulated: np.ndarray
    experimental: np.ndarray
    residual: np.ndarray
    rmse: float
    mae: float
    correlation: float
    simulated_map: object = None
    experimental_map: object = None


def load_intensity_map(path, delimiter=None):
    """Load an experimental intensity map from .npy or text/CSV."""

    path = str(path)
    if path.lower().endswith(".npy"):
        return np.load(path)
    return np.loadtxt(path, delimiter=delimiter)


def _as_float_image(image):
    arr = np.asarray(image)
    if np.issubdtype(arr.dtype, np.integer):
        arr = arr.astype(float) / float(np.iinfo(arr.dtype).max)
    else:
        arr = arr.astype(float)
    return arr


def image_to_grayscale(image, channel_weights=(0.299, 0.587, 0.114)):
    """Convert a grayscale/RGB/RGBA image array to a 2D intensity array."""

    arr = _as_float_image(image)
    if arr.ndim == 2:
        return arr
    if arr.ndim != 3 or arr.shape[2] < 3:
        raise ValueError("image must be a 2D grayscale array or an RGB/RGBA array")
    weights = np.asarray(channel_weights, dtype=float)
    if weights.shape != (3,):
        raise ValueError("channel_weights must contain three RGB weights")
    weights = weights / np.sum(weights)
    return np.tensordot(arr[:, :, :3], weights, axes=([2], [0]))


def crop_image(image, crop=None):
    """Crop an image with ``(row_start, row_stop, col_start, col_stop)``."""

    arr = np.asarray(image)
    if crop is None:
        return arr
    if len(crop) != 4:
        raise ValueError("crop must be (row_start, row_stop, col_start, col_stop)")
    row_start, row_stop, col_start, col_stop = [int(value) for value in crop]
    if row_start < 0 or col_start < 0 or row_stop <= row_start or col_stop <= col_start:
        raise ValueError("crop bounds must be positive and ordered")
    if row_stop > arr.shape[0] or col_stop > arr.shape[1]:
        raise ValueError("crop bounds exceed image dimensions")
    return arr[row_start:row_stop, col_start:col_stop]


def scattering_map_from_image(
    image,
    q_range,
    energy_range,
    crop=None,
    invert=False,
    origin="upper",
    normalization=None,
    channel_weights=(0.299, 0.587, 0.114),
    metadata=None,
):
    """Calibrate an image array into a :class:`ScatteringMap`.

    ``q_range`` is mapped left-to-right across image columns.  ``energy_range``
    is mapped bottom-to-top for ``origin="upper"`` images, matching ordinary
    image files where row zero is the top of the image.
    """

    if len(q_range) != 2 or len(energy_range) != 2:
        raise ValueError("q_range and energy_range must each contain two values")
    if origin not in ("upper", "lower"):
        raise ValueError("origin must be 'upper' or 'lower'")

    cropped = crop_image(image, crop=crop)
    intensity = image_to_grayscale(cropped, channel_weights=channel_weights)
    if invert:
        intensity = np.nanmax(intensity) - intensity
    if origin == "upper":
        intensity = intensity[::-1, :]

    if normalization is not None:
        intensity = normalize_intensity(intensity, normalization)

    n_energy, n_q = intensity.shape
    q_axis = np.linspace(float(q_range[0]), float(q_range[1]), n_q)
    energy_axis = np.linspace(float(energy_range[0]), float(energy_range[1]), n_energy)
    map_metadata = dict(metadata or {})
    map_metadata.update(
        {
            "source": "calibrated image",
            "crop": crop,
            "origin": origin,
            "inverted": bool(invert),
            "normalization": normalization,
        }
    )
    return ScatteringMap(q_axis, energy_axis, intensity.T, map_metadata)


def load_scattering_map_image(
    path,
    q_range,
    energy_range,
    crop=None,
    invert=False,
    origin="upper",
    normalization=None,
    channel_weights=(0.299, 0.587, 0.114),
    metadata=None,
):
    """Load an image file and calibrate it into a :class:`ScatteringMap`."""

    try:
        from matplotlib import image as mpimg
    except ImportError as exc:
        raise ImportError("matplotlib is required to load image intensity maps") from exc

    image = mpimg.imread(path)
    map_metadata = dict(metadata or {})
    map_metadata["image_path"] = str(path)
    return scattering_map_from_image(
        image,
        q_range=q_range,
        energy_range=energy_range,
        crop=crop,
        invert=invert,
        origin=origin,
        normalization=normalization,
        channel_weights=channel_weights,
        metadata=map_metadata,
    )


def normalize_intensity(intensity, method="max", eps=1e-30):
    arr = np.asarray(intensity, dtype=float)
    if method == "none" or method is None:
        return arr.copy()
    if method == "max":
        scale = np.nanmax(np.abs(arr))
        return arr / (scale + eps)
    if method == "area":
        scale = np.nansum(np.abs(arr))
        return arr / (scale + eps)
    if method == "zscore":
        return (arr - np.nanmean(arr)) / (np.nanstd(arr) + eps)
    raise ValueError("method must be 'none', 'max', 'area', or 'zscore'")


def interpolate_energy_axis(intensity, source_energy, target_energy, axis=-1):
    source = np.asarray(source_energy, dtype=float)
    target = np.asarray(target_energy, dtype=float)
    func = interp1d(source, intensity, axis=axis, bounds_error=False, fill_value=0.0)
    return func(target)


def align_intensity_map(intensity, source_q, source_energy, target_q, target_energy):
    """Interpolate a 2D q-energy map onto target experimental axes."""

    arr = np.asarray(intensity, dtype=float)
    source_q = np.asarray(source_q, dtype=float)
    source_energy = np.asarray(source_energy, dtype=float)
    target_q = np.asarray(target_q, dtype=float)
    target_energy = np.asarray(target_energy, dtype=float)

    if arr.shape != (source_q.size, source_energy.size):
        raise ValueError("intensity shape must match source_q and source_energy lengths")

    energy_aligned = interpolate_energy_axis(arr, source_energy, target_energy, axis=1)
    if source_q.size == 1:
        q_aligned = np.zeros((target_q.size, target_energy.size), dtype=float)
        matching = np.isclose(target_q, source_q[0])
        q_aligned[matching] = energy_aligned[0]
        return q_aligned

    q_func = interp1d(source_q, energy_aligned, axis=0, bounds_error=False, fill_value=0.0)
    return q_func(target_q)


def _axis_sigma_points(axis, sigma):
    if sigma is None or sigma == 0:
        return 0.0
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError("axis must contain at least two points")
    diffs = np.diff(axis)
    if np.any(diffs == 0):
        raise ValueError("axis values must be distinct")
    spacing = float(np.mean(np.abs(diffs)))
    return abs(float(sigma)) / spacing


def broaden_energy(intensity, sigma_points, axis=-1):
    return gaussian_filter1d(intensity, sigma=sigma_points, axis=axis, mode="nearest")


def convolve_q_resolution(intensity, sigma_q_points, sigma_energy_points=0.0):
    sigma = [sigma_q_points] * np.asarray(intensity).ndim
    sigma[-1] = sigma_energy_points
    return gaussian_filter(intensity, sigma=sigma, mode="nearest")


def convolve_resolution_units(scattering_map, sigma_q=None, sigma_energy=None):
    """Apply Gaussian q/energy resolution using physical axis units."""

    if not isinstance(scattering_map, ScatteringMap):
        raise TypeError("scattering_map must be a ScatteringMap")
    sigma_q_points = _axis_sigma_points(scattering_map.q_axis, sigma_q)
    sigma_energy_points = _axis_sigma_points(scattering_map.energy_axis, sigma_energy)
    intensity = gaussian_filter(
        scattering_map.intensity,
        sigma=(sigma_q_points, sigma_energy_points),
        mode="nearest",
    )
    metadata = dict(scattering_map.metadata)
    metadata.update(
        {
            "sigma_q": sigma_q,
            "sigma_energy": sigma_energy,
            "sigma_q_points": sigma_q_points,
            "sigma_energy_points": sigma_energy_points,
        }
    )
    return ScatteringMap(scattering_map.q_axis, scattering_map.energy_axis, intensity, metadata)


def apply_quantum_correction_to_map(scattering_map, temperature, energy_unit="thz"):
    """Apply the harmonic classical-to-quantum correction to a map."""

    if not isinstance(scattering_map, ScatteringMap):
        raise TypeError("scattering_map must be a ScatteringMap")
    intensity = apply_quantum_correction(
        scattering_map.intensity,
        scattering_map.energy_axis,
        temperature,
        unit=energy_unit,
        axis_index=1,
    )
    metadata = dict(scattering_map.metadata)
    metadata.update(
        {
            "quantum_correction": "harmonic classical-to-quantum detailed-balance factor",
            "quantum_temperature_kelvin": float(temperature),
            "energy_unit": energy_unit,
        }
    )
    return ScatteringMap(scattering_map.q_axis, scattering_map.energy_axis, intensity, metadata)


def align_to_experiment(simulated_map, experimental_map):
    """Return the simulated map interpolated onto the experimental axes."""

    if not isinstance(simulated_map, ScatteringMap) or not isinstance(experimental_map, ScatteringMap):
        raise TypeError("simulated_map and experimental_map must be ScatteringMap instances")
    aligned = align_intensity_map(
        simulated_map.intensity,
        simulated_map.q_axis,
        simulated_map.energy_axis,
        experimental_map.q_axis,
        experimental_map.energy_axis,
    )
    metadata = dict(simulated_map.metadata)
    metadata.update({"aligned_to": "experimental axes"})
    return ScatteringMap(experimental_map.q_axis, experimental_map.energy_axis, aligned, metadata)


def compare_maps(simulated, experimental, normalization="max"):
    sim = normalize_intensity(simulated, normalization)
    exp = normalize_intensity(experimental, normalization)
    if sim.shape != exp.shape:
        raise ValueError("simulated and experimental maps must have the same shape")
    residual = sim - exp
    rmse = float(np.sqrt(np.nanmean(residual ** 2)))
    mae = float(np.nanmean(np.abs(residual)))
    sim_flat = sim.ravel()
    exp_flat = exp.ravel()
    if np.nanstd(sim_flat) == 0 or np.nanstd(exp_flat) == 0:
        corr = np.nan
    else:
        corr = float(np.corrcoef(sim_flat, exp_flat)[0, 1])
    return MapComparison(sim, exp, residual, rmse, mae, corr)


def prepare_map_comparison(
    simulated_map,
    experimental_map,
    sigma_q=None,
    sigma_energy=None,
    normalization="max",
    quantum_temperature=None,
    energy_unit="thz",
):
    """Align, optionally quantum-correct, broaden, normalize, and compare maps."""

    aligned = align_to_experiment(simulated_map, experimental_map)
    if quantum_temperature is not None:
        aligned = apply_quantum_correction_to_map(
            aligned,
            temperature=quantum_temperature,
            energy_unit=energy_unit,
        )
    broadened = convolve_resolution_units(aligned, sigma_q=sigma_q, sigma_energy=sigma_energy)
    comparison = compare_maps(
        broadened.intensity,
        experimental_map.intensity,
        normalization=normalization,
    )
    comparison.simulated_map = broadened
    comparison.experimental_map = experimental_map
    return comparison


def linecut(intensity, index, axis=0):
    return np.take(np.asarray(intensity), indices=index, axis=axis)


def track_peak_positions(intensity, energy_axis, axis=-1):
    arr = np.asarray(intensity)
    energy = np.asarray(energy_axis)
    peak_indices = np.argmax(arr, axis=axis)
    return energy[peak_indices]


def _finite_float_or_none(value):
    value = float(value)
    return value if np.isfinite(value) else None


def _json_safe(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    return value


def _comparison_axes(comparison):
    if comparison.experimental_map is not None:
        q_axis = np.asarray(comparison.experimental_map.q_axis, dtype=float)
        energy_axis = np.asarray(comparison.experimental_map.energy_axis, dtype=float)
    else:
        q_axis = np.arange(np.asarray(comparison.experimental).shape[0], dtype=float)
        energy_axis = np.arange(np.asarray(comparison.experimental).shape[1], dtype=float)
    return q_axis, energy_axis


def comparison_summary(comparison):
    """Return JSON-ready scalar metrics for a map comparison."""

    simulated = np.asarray(comparison.simulated, dtype=float)
    experimental = np.asarray(comparison.experimental, dtype=float)
    residual = np.asarray(comparison.residual, dtype=float)
    q_axis, energy_axis = _comparison_axes(comparison)
    summary = {
        "rmse": _finite_float_or_none(comparison.rmse),
        "mae": _finite_float_or_none(comparison.mae),
        "correlation": _finite_float_or_none(comparison.correlation),
        "simulated_shape": list(simulated.shape),
        "experimental_shape": list(experimental.shape),
        "residual_max_abs": _finite_float_or_none(np.nanmax(np.abs(residual))) if residual.size else 0.0,
        "residual_mean": _finite_float_or_none(np.nanmean(residual)) if residual.size else 0.0,
        "q_min": _finite_float_or_none(np.nanmin(q_axis)) if q_axis.size else None,
        "q_max": _finite_float_or_none(np.nanmax(q_axis)) if q_axis.size else None,
        "energy_min": _finite_float_or_none(np.nanmin(energy_axis)) if energy_axis.size else None,
        "energy_max": _finite_float_or_none(np.nanmax(energy_axis)) if energy_axis.size else None,
        "has_simulated_map": comparison.simulated_map is not None,
        "has_experimental_map": comparison.experimental_map is not None,
    }
    if comparison.simulated_map is not None:
        summary["simulated_metadata"] = _json_safe(dict(comparison.simulated_map.metadata))
    if comparison.experimental_map is not None:
        summary["experimental_metadata"] = _json_safe(dict(comparison.experimental_map.metadata))
    return summary


def comparison_peak_records(comparison):
    """Return peak-tracking rows along the energy axis for every q point."""

    q_axis, energy_axis = _comparison_axes(comparison)
    simulated = np.asarray(comparison.simulated, dtype=float)
    experimental = np.asarray(comparison.experimental, dtype=float)
    if simulated.shape != experimental.shape:
        raise ValueError("simulated and experimental arrays must have the same shape")
    if simulated.shape != (q_axis.size, energy_axis.size):
        raise ValueError("comparison arrays must match q and energy axes")

    sim_peak_index = np.argmax(simulated, axis=1)
    exp_peak_index = np.argmax(experimental, axis=1)
    records = []
    for index, q_value in enumerate(q_axis):
        sim_idx = int(sim_peak_index[index])
        exp_idx = int(exp_peak_index[index])
        sim_energy = float(energy_axis[sim_idx])
        exp_energy = float(energy_axis[exp_idx])
        records.append(
            {
                "q_index": int(index),
                "q_value": float(q_value),
                "simulated_peak_energy": sim_energy,
                "experimental_peak_energy": exp_energy,
                "peak_energy_delta": sim_energy - exp_energy,
                "simulated_peak_intensity": float(simulated[index, sim_idx]),
                "experimental_peak_intensity": float(experimental[index, exp_idx]),
            }
        )
    return records


def _coerce_indices(indices, axis_size, name):
    if indices is None:
        return []
    if np.isscalar(indices):
        indices = [indices]
    output = []
    for value in indices:
        index = int(value)
        if index < 0:
            index += axis_size
        if index < 0 or index >= axis_size:
            raise IndexError("%s index %d is out of bounds" % (name, int(value)))
        output.append(index)
    return output


def comparison_linecut_records(comparison, q_indices=None, energy_indices=None):
    """Return fixed-q and fixed-energy line-cut records.

    The records use the normalized arrays stored in
    :class:`pySED.compare_experiment.MapComparison`, so they match the residual
    and scalar metrics used in the paper-ready comparison.
    """

    q_axis, energy_axis = _comparison_axes(comparison)
    simulated = np.asarray(comparison.simulated, dtype=float)
    experimental = np.asarray(comparison.experimental, dtype=float)
    residual = np.asarray(comparison.residual, dtype=float)
    if simulated.shape != experimental.shape or simulated.shape != residual.shape:
        raise ValueError("comparison arrays must have the same shape")
    if simulated.shape != (q_axis.size, energy_axis.size):
        raise ValueError("comparison arrays must match q and energy axes")

    q_indices = _coerce_indices(q_indices, q_axis.size, "q")
    energy_indices = _coerce_indices(energy_indices, energy_axis.size, "energy")
    records = []
    for q_index in q_indices:
        for energy_index, energy_value in enumerate(energy_axis):
            records.append(
                {
                    "cut_axis": "energy_at_q",
                    "cut_index": int(q_index),
                    "q_index": int(q_index),
                    "energy_index": int(energy_index),
                    "q_value": float(q_axis[q_index]),
                    "energy_value": float(energy_value),
                    "simulated": float(simulated[q_index, energy_index]),
                    "experimental": float(experimental[q_index, energy_index]),
                    "residual": float(residual[q_index, energy_index]),
                }
            )
    for energy_index in energy_indices:
        for q_index, q_value in enumerate(q_axis):
            records.append(
                {
                    "cut_axis": "q_at_energy",
                    "cut_index": int(energy_index),
                    "q_index": int(q_index),
                    "energy_index": int(energy_index),
                    "q_value": float(q_value),
                    "energy_value": float(energy_axis[energy_index]),
                    "simulated": float(simulated[q_index, energy_index]),
                    "experimental": float(experimental[q_index, energy_index]),
                    "residual": float(residual[q_index, energy_index]),
                }
            )
    return records


def map_value_records(scattering_map, value_name="intensity"):
    """Return flattened q-energy map records for CSV export."""

    if not isinstance(scattering_map, ScatteringMap):
        raise TypeError("scattering_map must be a ScatteringMap")
    records = []
    for i_q, q_value in enumerate(scattering_map.q_axis):
        for i_energy, energy_value in enumerate(scattering_map.energy_axis):
            records.append(
                {
                    "q_index": int(i_q),
                    "energy_index": int(i_energy),
                    "q_value": float(q_value),
                    "energy_value": float(energy_value),
                    value_name: float(scattering_map.intensity[i_q, i_energy]),
                }
            )
    return records


def comparison_residual_map(comparison):
    """Return the residual as a calibrated :class:`ScatteringMap`."""

    q_axis, energy_axis = _comparison_axes(comparison)
    metadata = {
        "source": "pySED map comparison residual",
        "definition": "normalized simulated intensity minus normalized experimental intensity",
        "rmse": _finite_float_or_none(comparison.rmse),
        "mae": _finite_float_or_none(comparison.mae),
        "correlation": _finite_float_or_none(comparison.correlation),
    }
    return ScatteringMap(q_axis, energy_axis, comparison.residual, metadata)


def _write_records_csv(path, records):
    records = list(records)
    fieldnames = list(records[0].keys()) if records else []
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(records)


def write_map_comparison_report(
    comparison,
    output_dir,
    prefix="comparison",
    write_residual=True,
    write_peaks=True,
    linecut_q_indices=None,
    linecut_energy_indices=None,
):
    """Write paper-ready comparison metrics and tables.

    The report consists of a JSON summary plus optional flattened CSV tables for
    the residual map, q-resolved peak tracking, and fixed-q/fixed-energy line
    cuts.
    """

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    prefix = str(prefix)
    written = {}

    summary_path = output / ("%s_summary.json" % prefix)
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(comparison_summary(comparison), handle, indent=2, sort_keys=True)
    written["summary"] = str(summary_path)

    if write_residual:
        residual_path = output / ("%s_residual.csv" % prefix)
        _write_records_csv(
            residual_path,
            map_value_records(comparison_residual_map(comparison), value_name="residual"),
        )
        written["residual"] = str(residual_path)

    if write_peaks:
        peaks_path = output / ("%s_peaks.csv" % prefix)
        _write_records_csv(peaks_path, comparison_peak_records(comparison))
        written["peaks"] = str(peaks_path)

    linecut_records = comparison_linecut_records(
        comparison,
        q_indices=linecut_q_indices,
        energy_indices=linecut_energy_indices,
    )
    if linecut_records:
        linecuts_path = output / ("%s_linecuts.csv" % prefix)
        _write_records_csv(linecuts_path, linecut_records)
        written["linecuts"] = str(linecuts_path)

    return written
