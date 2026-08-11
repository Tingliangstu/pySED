"""Tests for independent zero-background SED peak fitting."""

import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

from pySED.parser import read_input
from pySED.sed.FileIO import (
    hwhm_to_lifetime_ps,
    write_phonon_lifetime,
)
from pySED.sed.Plot_SED import velocity_dho as plotted_velocity_dho
from pySED.sed.Lorentz import (
    _find_sed_peaks,
    _fit_single_peak,
    _local_peak_bases,
    _lorentzian,
    _peak_fit_interval,
    _select_lorentz_frequency_range,
    _velocity_dho,
)


def _synthetic_sed():
    frequency = np.linspace(0.0, 20.0, 2001)
    baseline = 2.0e-4 * np.exp(-frequency / 15.0)
    peak_1 = 2.0e-2 / (1.0 + ((frequency - 4.0) / 0.12) ** 2)
    peak_2 = 1.5e-3 / (1.0 + ((frequency - 13.0) / 0.20) ** 2)
    ripple = 1.0 + 0.01 * np.sin(2.0 * np.pi * frequency / 0.37)
    return frequency, (baseline + peak_1 + peak_2) * ripple


class IndependentPeakFitTests(unittest.TestCase):
    def test_hwhm_conversion_matches_reported_time_convention(self):
        hwhm = np.array([0.01, 0.02, 0.0])

        lifetime = hwhm_to_lifetime_ps(hwhm)

        np.testing.assert_allclose(
            lifetime[:2], 1.0 / (2.0 * np.pi * hwhm[:2])
        )
        self.assertTrue(np.isnan(lifetime[2]))

    def test_significance_detection_is_invariant_to_linear_scaling(self):
        _, sed = _synthetic_sed()

        peaks, significance = _find_sed_peaks(
            sed, peak_min_significance=4.0
        )
        scaled_peaks, scaled_significance = _find_sed_peaks(
            1000.0 * sed, peak_min_significance=4.0
        )

        np.testing.assert_array_equal(peaks, scaled_peaks)
        np.testing.assert_allclose(
            significance,
            scaled_significance,
            rtol=1.0e-10,
            atol=1.0e-10,
        )

    def test_peak_detection_defaults_to_significance_four(self):
        _, sed = _synthetic_sed()

        default_peaks, default_significance = _find_sed_peaks(sed)
        explicit_peaks, explicit_significance = _find_sed_peaks(
            sed,
            peak_min_significance=4.0,
        )

        np.testing.assert_array_equal(default_peaks, explicit_peaks)
        np.testing.assert_allclose(
            default_significance,
            explicit_significance,
        )

    def test_lorentz_frequency_range_is_inclusive(self):
        frequency = np.arange(0.0, 8.0)
        sed = 10.0 + frequency

        selected_frequency, selected_sed = _select_lorentz_frequency_range(
            frequency, sed, freq_min=2.0, freq_max=5.0
        )

        np.testing.assert_array_equal(
            selected_frequency, np.array([2.0, 3.0, 4.0, 5.0])
        )
        np.testing.assert_array_equal(
            selected_sed, np.array([12.0, 13.0, 14.0, 15.0])
        )

    def test_parser_supports_only_auto_lorentz_dho_line_shape(self):
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            for function in ("lorentz", "dho", "auto"):
                with self.subTest(function=function):
                    input_file.write_text(
                        f"action = fit\nfitting_function = {function}\n",
                        encoding="utf-8",
                    )
                    params = read_input(input_file)
                    self.assertEqual(params.fitting_function, function)

            for legacy_key in (
                "fit_baseline_model",
                "fit_peak_strategy",
                "peak_height",
                "peak_prominence",
            ):
                input_file.write_text(
                    f"action = fit\n{legacy_key} = auto\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(ValueError, legacy_key):
                    read_input(input_file)

            input_file.write_text("action = fit\n", encoding="utf-8")
            params = read_input(input_file)
            self.assertEqual(params.fitting_function, "auto")
            self.assertEqual(params.peak_min_significance, 4.0)

            for invalid in ("0", "-1", "nan"):
                input_file.write_text(
                    "action = fit\n"
                    f"peak_min_significance = {invalid}\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(
                    ValueError,
                    "peak_min_significance must be positive and finite",
                ):
                    read_input(input_file)

    def test_parser_reads_and_validates_lorentz_frequency_range(self):
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text(
                "action = fit\n"
                "lorentz_fit_freq_min = 2.0\n"
                "lorentz_fit_freq_max = 5.0\n",
                encoding="utf-8",
            )
            params = read_input(input_file)
            self.assertEqual(params.lorentz_fit_freq_min, 2.0)
            self.assertEqual(params.lorentz_fit_freq_max, 5.0)

    def test_prominence_base_window_uses_all_points_between_bases(self):
        self.assertEqual(_peak_fit_interval(10, 30, 20, 100, 0), (10, 29))
        self.assertEqual(_peak_fit_interval(10, 30, 20, 100, 2), (12, 27))

    def test_significance_window_uses_complete_smoothed_peak_basin(self):
        valley_signal = np.array(
            [8.0, 6.0, 7.0, 1.0, 5.0, 9.0, 6.0, 2.0, 7.0]
        )
        peaks = np.array([1, 5, 8])

        bases = _local_peak_bases(
            valley_signal,
            peaks,
            peak_number=1,
        )

        self.assertEqual(bases, (3, 7))

    def test_outer_peak_basin_uses_smoothed_edge_valley(self):
        valley_signal = np.array([1.0, 8.0, 6.0, 2.0, 9.0, 5.0, 0.0])
        peaks = np.array([1, 4])

        first_bases = _local_peak_bases(
            valley_signal,
            peaks,
            peak_number=0,
        )
        last_bases = _local_peak_bases(
            valley_signal,
            peaks,
            peak_number=1,
        )

        self.assertEqual(first_bases, (0, 3))
        self.assertEqual(last_bases, (3, 6))

    def test_velocity_dho_uses_peak_amplitude_and_equivalent_hwhm(self):
        center = 2.0
        amplitude = 3.0
        hwhm = 0.25
        damping_width = 2.0 * hwhm
        root = np.sqrt(4.0 * center**2 + damping_width**2)
        half_max_frequencies = np.array(
            [
                0.5 * (root - damping_width),
                0.5 * (root + damping_width),
            ]
        )

        self.assertAlmostEqual(
            _velocity_dho(center, center, amplitude, hwhm), amplitude
        )
        np.testing.assert_allclose(
            _velocity_dho(half_max_frequencies, center, amplitude, hwhm),
            0.5 * amplitude,
        )
        np.testing.assert_allclose(
            plotted_velocity_dho(
                half_max_frequencies,
                center,
                amplitude,
                hwhm,
            ),
            _velocity_dho(half_max_frequencies, center, amplitude, hwhm),
        )
        self.assertAlmostEqual(np.diff(half_max_frequencies)[0], 2.0 * hwhm)

    def test_lorentz_fit_recovers_zero_background_peak(self):
        frequency = np.linspace(1.0, 3.0, 401)
        expected = np.array([2.0, 1.0, 0.08])
        sed = _lorentzian(frequency, *expected)
        peak = int(np.argmax(sed))

        result = _fit_single_peak(
            frequency,
            sed,
            peak,
            0,
            len(frequency) - 1,
            initial_hwhm=0.01,
            peak_max_hwhm=1.0,
            fitting_function="lorentz",
        )

        self.assertEqual(result["fitting_function"], "lorentz")
        self.assertFalse(result["incomplete_peak_shapes"][0])
        np.testing.assert_allclose(
            result["peak_parameters"][0], expected, rtol=1.0e-6, atol=1.0e-8
        )

    def test_fit_is_invariant_to_small_sed_amplitude_scale(self):
        frequency = np.linspace(1.0, 3.0, 401)
        expected = np.array([2.0, 1.0e-6, 0.08])
        sed = _lorentzian(frequency, *expected)
        peak = int(np.argmax(sed))

        result = _fit_single_peak(
            frequency,
            sed,
            peak,
            0,
            len(frequency) - 1,
            initial_hwhm=0.20,
            peak_max_hwhm=1.0,
            fitting_function="lorentz",
        )

        np.testing.assert_allclose(
            result["peak_parameters"][0], expected, rtol=1.0e-6, atol=1.0e-10
        )

    def test_fit_flags_a_peak_without_both_half_maximum_points(self):
        frequency = np.linspace(1.0, 3.0, 401)
        expected = np.array([2.0, 1.0, 0.20])
        sed = _lorentzian(frequency, *expected)
        peak = int(np.argmax(sed))
        start = int(np.argmin(np.abs(frequency - 1.90)))

        result = _fit_single_peak(
            frequency,
            sed,
            peak,
            start,
            len(frequency) - 1,
            initial_hwhm=0.05,
            peak_max_hwhm=1.0,
            fitting_function="lorentz",
        )

        self.assertTrue(result["incomplete_peak_shapes"][0])

    def test_auto_selects_dho_for_broad_zero_background_velocity_spectrum(self):
        frequency = np.linspace(0.05, 4.0, 500)
        expected = np.array([1.0, 2.0, 0.4])
        sed = _velocity_dho(frequency, *expected)
        peak = int(np.argmax(sed))

        result = _fit_single_peak(
            frequency,
            sed,
            peak,
            0,
            len(frequency) - 1,
            initial_hwhm=expected[2],
            peak_max_hwhm=2.0,
            fitting_function="auto",
        )

        self.assertEqual(result["fitting_function"], "dho")
        np.testing.assert_allclose(
            result["peak_parameters"][0], expected, rtol=1.0e-6, atol=1.0e-8
        )

    def test_auto_prefers_lorentz_in_weak_damping_limit(self):
        frequency = np.linspace(1.0, 3.0, 500)
        expected = np.array([2.0, 1.0, 0.02])
        sed = _velocity_dho(frequency, *expected)
        peak = int(np.argmax(sed))

        result = _fit_single_peak(
            frequency,
            sed,
            peak,
            0,
            len(frequency) - 1,
            initial_hwhm=expected[2],
            peak_max_hwhm=1.0,
            fitting_function="auto",
        )

        self.assertEqual(result["fitting_function"], "lorentz")

    def test_overdamped_dho_is_written_without_diagnostic_side_files(self):
        fit = SimpleNamespace(
            popt=np.array([[1.0, 2.0, 0.2], [0.1, 1.0, 0.2]]),
            fit_models=np.array(["dho", "dho"], dtype=object),
        )
        params = SimpleNamespace(qpoint_slice_index=0)

        with TemporaryDirectory() as directory:
            original = os.getcwd()
            try:
                os.chdir(directory)
                lifetime_file = write_phonon_lifetime(fit, params)
                lifetime_lines = Path(lifetime_file).read_text(
                    encoding="utf-8"
                ).splitlines()
                side_files = [
                    *Path(directory).glob("*.params"),
                    *Path(directory).glob("*.error"),
                    *Path(directory).glob("*.models"),
                ]
            finally:
                os.chdir(original)

        self.assertEqual(
            Path(lifetime_file),
            Path("Lifetime") / "Fitting-0-qpoint.Fre_lifetime",
        )
        self.assertEqual(side_files, [])
        self.assertEqual(len(lifetime_lines), 4)
        self.assertTrue(lifetime_lines[2].startswith("1.000000 "))
        self.assertTrue(lifetime_lines[3].startswith("0.100000 "))


if __name__ == "__main__":
    unittest.main()
