import numpy as np

from pySED.compare_experiment import (
    ScatteringMap,
    align_intensity_map,
    apply_quantum_correction_to_map,
    compare_maps,
    comparison_linecut_records,
    comparison_peak_records,
    comparison_residual_map,
    comparison_summary,
    convolve_resolution_units,
    crop_image,
    image_to_grayscale,
    load_scattering_map_image,
    normalize_intensity,
    prepare_map_comparison,
    scattering_map_from_image,
    track_peak_positions,
    write_map_comparison_report,
)


def test_compare_maps_identical_after_normalization():
    arr = np.array([[0.0, 1.0], [2.0, 3.0]])
    comparison = compare_maps(arr, arr * 2.0, normalization="max")

    assert comparison.rmse == 0.0
    assert comparison.mae == 0.0
    assert comparison.correlation > 0.999


def test_peak_tracking_returns_energy_values():
    energy = np.array([0.0, 1.0, 2.0])
    intensity = np.array([[0.0, 5.0, 1.0], [2.0, 1.0, 0.0]])
    np.testing.assert_allclose(track_peak_positions(intensity, energy), [1.0, 0.0])


def test_normalize_area():
    arr = np.array([1.0, 1.0, 2.0])
    np.testing.assert_allclose(normalize_intensity(arr, "area"), arr / 4.0)


def test_image_to_grayscale_and_crop():
    image = np.zeros((3, 4, 3), dtype=float)
    image[:, :, 0] = 1.0
    gray = image_to_grayscale(image)
    cropped = crop_image(gray, crop=(1, 3, 1, 4))

    np.testing.assert_allclose(gray, 0.299)
    assert cropped.shape == (2, 3)


def test_scattering_map_from_image_calibrates_axes_origin_and_inversion():
    image = np.array(
        [
            [0.0, 0.25, 0.50],
            [0.50, 0.75, 1.00],
        ]
    )

    scattering_map = scattering_map_from_image(
        image,
        q_range=(0.0, 2.0),
        energy_range=(10.0, 20.0),
        invert=True,
        origin="upper",
    )

    np.testing.assert_allclose(scattering_map.q_axis, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(scattering_map.energy_axis, [10.0, 20.0])
    expected_after_invert_and_flip = np.array(
        [
            [0.50, 0.25, 0.00],
            [1.00, 0.75, 0.50],
        ]
    )
    np.testing.assert_allclose(scattering_map.intensity, expected_after_invert_and_flip.T)
    assert scattering_map.metadata["inverted"] is True


def test_load_scattering_map_image_from_png(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    image = np.array([[0.0, 0.5], [0.75, 1.0]])
    path = tmp_path / "experimental.png"
    plt.imsave(path, image, cmap="gray", vmin=0.0, vmax=1.0)

    scattering_map = load_scattering_map_image(
        path,
        q_range=(-0.1, 0.1),
        energy_range=(0.0, 50.0),
        normalization="max",
    )

    assert scattering_map.intensity.shape == (2, 2)
    np.testing.assert_allclose(scattering_map.q_axis, [-0.1, 0.1])
    np.testing.assert_allclose(scattering_map.energy_axis, [0.0, 50.0])
    assert scattering_map.metadata["image_path"].endswith("experimental.png")
    assert np.nanmax(scattering_map.intensity) == 1.0


def test_align_intensity_map_interpolates_q_and_energy_axes():
    source_q = np.array([0.0, 1.0])
    source_energy = np.array([0.0, 2.0])
    intensity = np.array([[0.0, 2.0], [10.0, 12.0]])

    aligned = align_intensity_map(
        intensity,
        source_q,
        source_energy,
        target_q=np.array([0.5]),
        target_energy=np.array([1.0]),
    )

    np.testing.assert_allclose(aligned, [[6.0]])


def test_resolution_convolution_uses_axis_units():
    scattering_map = ScatteringMap(
        q_axis=np.array([0.0, 0.5, 1.0]),
        energy_axis=np.array([0.0, 2.0, 4.0]),
        intensity=np.eye(3),
    )

    broadened = convolve_resolution_units(scattering_map, sigma_q=0.5, sigma_energy=2.0)

    assert broadened.intensity.shape == scattering_map.intensity.shape
    assert broadened.metadata["sigma_q_points"] == 1.0
    assert broadened.metadata["sigma_energy_points"] == 1.0
    assert broadened.intensity[0, 1] > 0


def test_prepare_map_comparison_aligns_broadens_and_compares():
    simulated = ScatteringMap(
        q_axis=np.array([0.0, 1.0]),
        energy_axis=np.array([0.0, 2.0]),
        intensity=np.array([[0.0, 2.0], [10.0, 12.0]]),
    )
    experimental = ScatteringMap(
        q_axis=np.array([0.5]),
        energy_axis=np.array([1.0]),
        intensity=np.array([[6.0]]),
    )

    comparison = prepare_map_comparison(
        simulated,
        experimental,
        sigma_q=0.0,
        sigma_energy=0.0,
        normalization="none",
    )

    np.testing.assert_allclose(comparison.simulated, [[6.0]])
    np.testing.assert_allclose(comparison.residual, [[0.0]])
    assert comparison.rmse == 0.0
    assert comparison.simulated_map.q_axis[0] == 0.5


def test_apply_quantum_correction_to_map_updates_intensity_and_metadata():
    scattering_map = ScatteringMap(
        q_axis=np.array([0.0]),
        energy_axis=np.array([0.0, 10.0]),
        intensity=np.ones((1, 2)),
    )

    corrected = apply_quantum_correction_to_map(scattering_map, temperature=300.0, energy_unit="meV")

    assert corrected.intensity[0, 0] == 1.0
    assert corrected.intensity[0, 1] > 1.0
    assert corrected.metadata["quantum_temperature_kelvin"] == 300.0


def test_prepare_map_comparison_can_apply_quantum_correction():
    simulated = ScatteringMap(
        q_axis=np.array([0.0]),
        energy_axis=np.array([0.0, 10.0]),
        intensity=np.ones((1, 2)),
    )
    experimental = ScatteringMap(
        q_axis=np.array([0.0]),
        energy_axis=np.array([0.0, 10.0]),
        intensity=np.ones((1, 2)),
    )

    comparison = prepare_map_comparison(
        simulated,
        experimental,
        normalization="none",
        quantum_temperature=300.0,
        energy_unit="meV",
    )

    assert comparison.simulated[0, 1] > comparison.experimental[0, 1]
    assert comparison.simulated_map.metadata["quantum_correction"].startswith("harmonic")


def test_comparison_report_exports_summary_residual_and_peaks(tmp_path):
    simulated = ScatteringMap(
        q_axis=np.array([0.0, 1.0]),
        energy_axis=np.array([0.0, 5.0, 10.0]),
        intensity=np.array([[0.0, 2.0, 1.0], [0.0, 1.0, 3.0]]),
        metadata={"numpy_value": np.float64(2.0), "array_value": np.array([1, 2])},
    )
    experimental = ScatteringMap(
        q_axis=np.array([0.0, 1.0]),
        energy_axis=np.array([0.0, 5.0, 10.0]),
        intensity=np.array([[0.0, 1.0, 2.0], [0.0, 1.0, 3.0]]),
    )
    comparison = prepare_map_comparison(
        simulated,
        experimental,
        sigma_q=0.0,
        sigma_energy=0.0,
        normalization="none",
    )

    summary = comparison_summary(comparison)
    residual_map = comparison_residual_map(comparison)
    peaks = comparison_peak_records(comparison)
    linecuts = comparison_linecut_records(comparison, q_indices=[0], energy_indices=[2])
    written = write_map_comparison_report(
        comparison,
        tmp_path,
        prefix="eels",
        linecut_q_indices=[0],
        linecut_energy_indices=[2],
    )

    assert summary["simulated_shape"] == [2, 3]
    assert summary["residual_max_abs"] == 1.0
    assert residual_map.metadata["definition"].startswith("normalized simulated")
    assert peaks[0]["simulated_peak_energy"] == 5.0
    assert peaks[0]["experimental_peak_energy"] == 10.0
    assert len(linecuts) == 5
    assert linecuts[0]["cut_axis"] == "energy_at_q"
    assert linecuts[-1]["cut_axis"] == "q_at_energy"
    assert set(written) == {"summary", "residual", "peaks", "linecuts"}

    import csv
    import json

    with open(written["summary"], "r", encoding="utf-8") as handle:
        loaded_summary = json.load(handle)
    with open(written["peaks"], "r", encoding="utf-8") as handle:
        peak_rows = list(csv.DictReader(handle))
    with open(written["residual"], "r", encoding="utf-8") as handle:
        residual_rows = list(csv.DictReader(handle))
    with open(written["linecuts"], "r", encoding="utf-8") as handle:
        linecut_rows = list(csv.DictReader(handle))

    assert loaded_summary["rmse"] == summary["rmse"]
    assert loaded_summary["simulated_metadata"]["numpy_value"] == 2.0
    assert loaded_summary["simulated_metadata"]["array_value"] == [1, 2]
    assert peak_rows[0]["peak_energy_delta"] == "-5.0"
    assert len(residual_rows) == 6
    assert len(linecut_rows) == 5
    assert linecut_rows[0]["energy_value"] == "0.0"
