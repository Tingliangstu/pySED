"""CLI parser-to-pipeline trajectory integration tests."""

import importlib.util
import io
import os
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np

from pySED.parser import read_input
from pySED.parameters import validate_qpoint_slice_index
from pySED.postprocess.output_files import sed_data_files
from pySED.pipeline import (
    _run_thinking,
    run,
    step_compute_sed,
    step_fit_sed,
    step_plot_sed,
    step_prepare_trajectory,
)
from pySED.postprocess import FileIO, Lorentz
from pySED.sed.construct_BZ import BZ_methods
from pySED.postprocess.Plot_SED import (
    _prepare_sed_for_log_scale,
    _resolve_color_limits,
    plot_bands,
    plot_lifetime_summary,
    plot_slice,
    resolve_slice_frequency_limit,
    resolve_slice_frequency_start,
)

NETCDF4_AVAILABLE = importlib.util.find_spec("netCDF4") is not None

GPUMD_TEXT = """\
1
Lattice="5 0 0 0 5 0 0 0 5" Properties=species:S:1:pos:R:3:vel:R:3
Si 0 0 0 0.001 0.002 0.003
"""


class ParserTests(unittest.TestCase):
    def test_plot_fails_when_required_sed_outputs_are_missing(self):
        with TemporaryDirectory() as directory:
            params = SimpleNamespace(
                out_files_name=str(Path(directory) / "missing"),
            )

            with self.assertRaisesRegex(FileNotFoundError, "SED"):
                step_plot_sed(params)

            for path in sed_data_files(params.out_files_name)[:3]:
                path.touch()
            with self.assertRaisesRegex(FileNotFoundError, "Q_distances_and_labels"):
                step_plot_sed(params)

    def test_fit_failure_propagates_without_a_false_success_message(self):
        with TemporaryDirectory() as directory:
            params = SimpleNamespace(
                action="fit",
                out_files_name=str(Path(directory) / "missing"),
            )
            output = io.StringIO()

            with (
                redirect_stdout(output),
                self.assertRaisesRegex(FileNotFoundError, "SED"),
            ):
                run(params)

        self.assertNotIn("finished successfully", output.getvalue())

    def test_thinking_fits_after_sed_data_and_requested_plots_exist(self):
        """Thinking should fit only after data and plot outputs are ready."""

        with TemporaryDirectory() as directory:
            previous_directory = os.getcwd()
            os.chdir(directory)
            try:
                output = Path(directory) / "SrTiO3"
                for suffix in (".SED", ".Qpts", ".THz", ".Q_distances_and_labels"):
                    output.with_suffix(suffix).touch()
                output.with_name("SrTiO3-SED.png").touch()

                params = SimpleNamespace(
                    out_files_name=str(output),
                    plot_partial_SED=False,
                    plot_slice=False,
                    qpoint_slice_index=3,
                    lorentz_fit_all_qpoint=False,
                )
                steps = {
                    "compute": Mock(),
                    "plot": Mock(),
                    "fit": Mock(),
                }

                _run_thinking(params, steps)

                steps["compute"].assert_not_called()
                steps["plot"].assert_not_called()
                steps["fit"].assert_called_once_with(params)
            finally:
                os.chdir(previous_directory)

    def test_explicit_fit_overwrites_existing_all_qpoint_results(self):
        """Explicit fitting processes every q point without mutating user input."""
        data = SimpleNamespace(q_distances=np.array([0.0, 0.5]))
        params = SimpleNamespace(
            action="fit", lorentz_fit_all_qpoint=True,
            qpoint_slice_index=0, if_show_figures=True,
        )
        original_params = vars(params).copy()
        total_file = "Lifetime/Fitting-All-Qpoints.Fre_lifetime"
        with (
            patch("pySED.pipeline._require_plot_data"),
            patch("pySED.pipeline._fit_ready", return_value=True) as ready,
            patch("pySED.pipeline.FileIO.load_data", return_value=data) as load,
            patch("pySED.pipeline._fit_one_qpoint",
                  return_value=SimpleNamespace(peak_fits=[])) as fit,
            patch("pySED.pipeline.FileIO.deal_total_fre_lifetime",
                  return_value=(total_file, 2)) as combine,
            patch("pySED.pipeline.FileIO.read_lifetime_data",
                  return_value=(np.array([1.0]), np.array([2.0]))),
            patch("pySED.pipeline.Plot_SED.plot_lifetime_summary") as plot,
            redirect_stdout(io.StringIO()),
        ):
            run(params)

        ready.assert_not_called()
        load.assert_called_once_with(params)
        self.assertEqual(fit.call_count, 2)
        for q_index, call in enumerate(fit.call_args_list):
            self.assertIs(call.args[0], data)
            self.assertIsNot(call.args[1], params)
            self.assertEqual(call.args[1].qpoint_slice_index, q_index)
            self.assertFalse(call.args[1].if_show_figures)
        combine.assert_called_once_with(params, 2)
        plot.assert_called_once()
        self.assertEqual(vars(params), original_params)

    def test_single_q_refit_can_rebuild_combined_lifetime_data(self):
        data = SimpleNamespace(
            qpoints=np.zeros((3, 3)),
            q_distances=np.arange(3.0),
        )
        params = SimpleNamespace(
            lorentz_fit_all_qpoint=False,
            qpoint_slice_index=1,
            re_output_total_freq_lifetime=True,
        )
        output = io.StringIO()

        with (
            patch("pySED.pipeline._require_plot_data"),
            patch(
                "pySED.pipeline.FileIO.load_data",
                return_value=data,
            ),
            patch("pySED.pipeline._fit_one_qpoint"),
            patch(
                "pySED.pipeline.FileIO.deal_total_fre_lifetime",
                return_value=(
                    "Lifetime/Fitting-All-Qpoints.Fre_lifetime",
                    12,
                ),
            ) as combine,
            patch(
                "pySED.pipeline.FileIO.read_lifetime_data",
                return_value=(np.array([1.0]), np.array([2.0])),
            ),
            patch(
                "pySED.pipeline.Plot_SED.plot_lifetime_summary",
                return_value="Fitting-Frequency-Lifetime.png",
            ) as plot_lifetime,
            redirect_stdout(output),
        ):
            step_fit_sed(params)

        combine.assert_called_once_with(params, 3)
        self.assertIn("Combined lifetime data rebuilt", output.getvalue())
        self.assertIn(
            "Lifetime/Fitting-All-Qpoints.Fre_lifetime (12 modes)",
            output.getvalue(),
        )
        self.assertIn("Single-Q spectral fit done", output.getvalue())
        plot_lifetime.assert_called_once()
        self.assertIn(
            "Lifetime figure redrawn: Fitting-Frequency-Lifetime.png",
            output.getvalue(),
        )
        self.assertFalse(hasattr(params, "lifetime_figure_action"))

    def test_lifetime_outputs_are_grouped_and_combined(self):
        with TemporaryDirectory() as directory:
            previous_directory = os.getcwd()
            os.chdir(directory)
            try:
                lifetime_directory = Path("Lifetime")
                lifetime_directory.mkdir()
                header = "Generated by pySED\nFrequency (THz) lifetime (ps)\n"
                (lifetime_directory / "Fitting-0-qpoint.Fre_lifetime").write_text(
                    header + "1.000000 2.00000000\n",
                    encoding="utf-8",
                )
                (lifetime_directory / "Fitting-1-qpoint.Fre_lifetime").write_text(
                    header + "3.000000 4.00000000\n",
                    encoding="utf-8",
                )
                combined_file = (
                    lifetime_directory / "Fitting-All-Qpoints.Fre_lifetime"
                )
                combined_file.write_text(
                    header + "99.000000 99.00000000\n",
                    encoding="utf-8",
                )

                output_file, total_modes = FileIO.deal_total_fre_lifetime(
                    SimpleNamespace(),
                    2,
                )
                frequency, lifetime = FileIO.read_lifetime_data(output_file)
                combined_text = combined_file.read_text(encoding="utf-8")
            finally:
                os.chdir(previous_directory)

        self.assertEqual(
            Path(output_file),
            Path("Lifetime") / "Fitting-All-Qpoints.Fre_lifetime",
        )
        self.assertEqual(total_modes, 2)
        np.testing.assert_allclose(frequency, [1.0, 3.0])
        np.testing.assert_allclose(lifetime, [2.0, 4.0])
        self.assertNotIn("99.000000", combined_text)

    def test_lifetime_summary_uses_frequency_and_log_lifetime_axes(self):
        with TemporaryDirectory() as directory:
            previous_directory = os.getcwd()
            os.chdir(directory)
            try:
                with (
                    patch("matplotlib.figure.Figure.savefig") as savefig,
                    patch("pySED.postprocess.Plot_SED.plt.close"),
                ):
                    output_path = plot_lifetime_summary(
                        np.array([1.0, 2.0, 3.0]),
                        np.array([10.0, 1.0, 0.1]),
                    )
                    figure = plt.gcf()
            finally:
                os.chdir(previous_directory)

        self.assertEqual(
            Path(output_path),
            Path("Fitting-Frequency-Lifetime.png"),
        )
        self.assertEqual(figure.axes[0].get_xlabel(), "Frequency (THz)")
        self.assertEqual(figure.axes[0].get_ylabel(), "Lifetime")
        self.assertEqual(figure.axes[0].get_yscale(), "log")
        self.assertFalse(
            any(
                line.get_visible()
                for line in (
                    figure.axes[0].get_xgridlines()
                    + figure.axes[0].get_ygridlines()
                )
            )
        )
        savefig.assert_called_once()
        plt.close(figure)

    def test_dispersion_colorbar_uses_dimensionless_reference_ratio(self):
        """The heatmap label should state its physical reference value."""

        data = SimpleNamespace(
            sed_avg=np.exp(
                np.array(
                    [
                        [-4.0, -3.5],
                        [-3.0, -2.5],
                        [-2.0, -1.5],
                    ]
                )
            ),
            qpoints=np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
            freq_fft=np.array([0.0, 1.0, 2.0]),
            q_distances=np.array([0.0, 1.0]),
            q_labels={0.0: "G", 1.0: "X"},
        )
        params = SimpleNamespace(
            colorbar_min=-4.0,
            colorbar_max=0.0,
            input_file="/work/SrTiO3.in",
            plot_color="RdBu_r",
            plot_interval=1.0,
            num_qpaths=1,
            use_contourf=False,
            plot_cutoff_freq=2.0,
            plot_partial_SED=False,
            out_files_name="SrTiO3",
            if_show_figures=False,
        )

        for unit in ("eV/THz", "legacy 2.x scale"):
            data.sed_unit = unit
            with (
                patch("pySED.postprocess.Plot_SED.plt.savefig"),
                redirect_stdout(io.StringIO()) as output,
            ):
                plot_bands(data, params)
                figure = plt.gcf()
            expected_label = (
                r"$\ln\!\left[\mathrm{SED}(\mathbf{q},\omega)"
                r"/(\mathrm{eV/THz})\right]$"
                if unit == "eV/THz" else "ln[SED / (legacy 2.x scale)]"
            )
            self.assertEqual(figure.axes[1].get_ylabel(), expected_label)
            self.assertIn(f"ln[SED / ({unit})]", output.getvalue())
            plt.close(figure)

    def test_slice_frequency_limit_uses_the_matching_plot_control(self):
        """Plain plots and Lorentz overlays should use their own cutoffs."""

        thz = np.array([0.0, 5.0, 10.0, 15.0])
        params = SimpleNamespace(
            plot_cutoff_freq=12.0,
            lorentz_fit_freq_min=2.0,
            lorentz_fit_freq_max=8.0,
        )
        self.assertEqual(
            resolve_slice_frequency_limit(thz, params),
            (12.0, "plot_cutoff_freq"),
        )

        self.assertEqual(
            resolve_slice_frequency_limit(thz, params, fitting=True),
            (8.0, "lorentz_fit_freq_max"),
        )
        self.assertEqual(
            resolve_slice_frequency_start(thz, params, fitting=True),
            (2.0, "lorentz_fit_freq_min"),
        )

    def test_qpoint_slice_index_uses_zero_based_data_bounds(self):
        """Slice indices must stay within the loaded q-point range."""

        data = SimpleNamespace(
            qpoints=np.zeros((10, 3)),
        )

        params = SimpleNamespace(qpoint_slice_index=10)
        with self.assertRaisesRegex(
            ValueError,
            r"qpoint_slice_index = 10 is out of range: "
            r"10 q-points are available; valid zero-based "
            r"indices are 0 to 9\.",
        ):
            validate_qpoint_slice_index(params, data.qpoints)

        for q_index in (0, 9):
            params = SimpleNamespace(qpoint_slice_index=q_index)
            validate_qpoint_slice_index(params, data.qpoints)

    def test_sed_color_scale_uses_natural_log_limits(self):
        """Dispersion colors should reproduce the original ln rendering."""

        params = SimpleNamespace(
            colorbar_min=None,
            colorbar_max=None,
        )
        sed = np.exp(np.array([-22.8, -16.0, -10.2]))
        vmin, vmax = _resolve_color_limits(
            sed,
            params,
        )
        self.assertEqual(vmin, -22.0)
        self.assertEqual(vmax, -10.0)

        rendered = _prepare_sed_for_log_scale(
            np.exp(np.array([[-30.0, -16.0, -5.0]])),
            vmin,
            vmax,
        )
        np.testing.assert_allclose(
            rendered,
            np.array([[vmin, -16.0, vmax]]),
        )

    def test_colorbar_limits_accept_natural_log_values(self):
        """Input validation should accept negative ln-scale limits."""

        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text(
                "action = plot\n"
                "colorbar_min = -12\n"
                "colorbar_max = 0\n",
                encoding="utf-8",
            )

            params = read_input(input_file)
            self.assertEqual(params.colorbar_min, -12.0)
            self.assertEqual(params.colorbar_max, 0.0)

    def test_partial_sed_files_use_element_names_and_ev_per_thz(self):
        """Partial outputs should be readable by their public element names."""

        with TemporaryDirectory() as directory:
            previous_directory = os.getcwd()
            os.chdir(directory)
            try:
                params = SimpleNamespace(
                    out_files_name="SrTiO3",
                    output_partial=1,
                    plot_partial_SED=True,
                    plot_partial_element="O",
                    plot_partial_dir="y",
                )
                phonons = SimpleNamespace(
                    sed_avg=np.arange(4 * 2 * 3 * 3, dtype=float).reshape(
                        4,
                        2,
                        3,
                        3,
                    ),
                    type_symbols=("O", "Ti", "Sr"),
                    freq_fft=np.arange(4, dtype=float),
                )
                lattice = SimpleNamespace(
                    reduced_qpoints=np.zeros((2, 3)),
                    q_distances=(0.0, 1.0),
                    q_labels=((0.0, "G"), (1.0, "X")),
                )
                partial_dir = Path("SrTiO3_partial_SED")
                partial_dir.mkdir()
                legacy_oxygen_y = (
                    partial_dir / "SrTiO3.SED_type1_y"
                )
                legacy_oxygen_y.write_text("old", encoding="utf-8")

                FileIO.write_output(phonons, params, lattice)

                oxygen_y = partial_dir / "SrTiO3.SED_O_y"
                self.assertTrue(oxygen_y.is_file())
                self.assertEqual(legacy_oxygen_y.read_text(encoding="utf-8"), "old")

                loaded = FileIO.load_data(params)
                np.testing.assert_allclose(
                    loaded.sed_avg,
                    phonons.sed_avg[:, :, 0, 1],
                )
                for sed_file in (Path("SrTiO3.SED"), oxygen_y):
                    header = sed_file.read_text(
                        encoding="utf-8"
                    ).splitlines()[0]
                    self.assertIn("unit: eV/THz", header)
                    self.assertIn(
                        "each column corresponds to one q-point "
                        "(same order as the .Qpts file)",
                        header,
                    )
            finally:
                os.chdir(previous_directory)

    def test_single_qpoint_output_reloads_for_plotting_and_fitting(self):
        """Total, element, and direction spectra retain their q-point axis."""
        frequency = np.linspace(1.0, 3.0, 401)
        spectrum = 1.0 / (1.0 + ((frequency - 2.0) / 0.08) ** 2)
        lattice = SimpleNamespace(
            reduced_qpoints=np.array([[0.25, 0.0, 0.0]]),
            q_distances=(0.0,), q_labels=((0.0, "Q"),),
        )
        with TemporaryDirectory() as directory:
            previous_directory = os.getcwd()
            try:
                os.chdir(directory)
                params = SimpleNamespace(
                    out_files_name="single", output_partial=False,
                    plot_partial_SED=False, plot_partial_element="Si",
                    plot_partial_dir=None, qpoint_slice_index=0,
                    lorentz_fit_freq_min=1.0, lorentz_fit_freq_max=3.0,
                    fitting_function="lorentz", peak_min_significance=4.0,
                    modulate_factor=0, initial_guess_hwhm=0.05,
                    peak_max_hwhm=1.0, if_show_figures=False,
                )
                phonons = SimpleNamespace(
                    sed_avg=spectrum[:, None], freq_fft=frequency,
                    type_symbols=("Si",),
                )
                for output_partial, legacy in (
                    (False, False), (False, True), (True, False), (True, True),
                ):
                    params.output_partial = output_partial
                    phonons.legacy_unit = legacy
                    params.legacy_unit = not legacy
                    phonons.sed_avg = spectrum[:, None]
                    if output_partial:
                        phonons.sed_avg = (
                            spectrum[:, None, None, None] * np.array([1, 2, 3])
                        )
                    FileIO.write_output(phonons, params, lattice)
                    selections = [(False, None, 6 if output_partial else 1)]
                    if output_partial:
                        selections += [(True, None, 6), (True, "x", 1),
                                       (True, "y", 2), (True, "z", 3)]
                    for partial, direction, amplitude in selections:
                        with self.subTest(partial=partial, direction=direction,
                                          output_partial=output_partial):
                            params.plot_partial_SED = partial
                            params.plot_partial_dir = direction
                            loaded = FileIO.load_data(params)
                            self.assertEqual(loaded.legacy_unit, legacy)
                            unit = "legacy 2.x scale" if legacy else "eV/THz"
                            self.assertEqual(loaded.sed_unit, unit)
                            np.testing.assert_allclose(
                                loaded.sed_avg, amplitude * spectrum[:, None],
                            )
                            np.testing.assert_array_equal(
                                loaded.qpoints, lattice.reduced_qpoints,
                            )
                            np.testing.assert_allclose(loaded.freq_fft, frequency)
                            validate_qpoint_slice_index(params, loaded.qpoints)
                            original_params = vars(params).copy()
                            original_sed = loaded.sed_avg.copy()
                            original_files = set(Path(directory).rglob("*"))
                            with redirect_stdout(io.StringIO()):
                                result = Lorentz.lorentz(loaded, params)
                            self.assertEqual(vars(params), original_params)
                            np.testing.assert_array_equal(loaded.sed_avg, original_sed)
                            self.assertEqual(set(Path(directory).rglob("*")), original_files)
                            np.testing.assert_allclose(
                                result.popt, [[2.0, amplitude, 0.08]], rtol=1e-5,
                            )
                            plotted_lines = []

                            def check_figure(*args, **kwargs):
                                self.assertIn(unit, plt.gca().get_ylabel())
                                plotted_lines.append(len(plt.gca().lines))

                            with patch("pySED.postprocess.Plot_SED.plt.savefig") as savefig:
                                savefig.side_effect = check_figure
                                plot_slice(loaded, params, fit_result=result)
                                plot_slice(loaded, params)
                            self.assertEqual(plotted_lines, [2, 1])
                            self.assertEqual(vars(params), original_params)
            finally:
                os.chdir(previous_directory)

    def test_legacy_unit_parser_defaults_and_boolean_spellings(self):
        """Legacy output is opt-in and accepts boolean words as well as 0/1."""
        with TemporaryDirectory() as directory:
            path = Path(directory) / "input.in"
            for value, expected in (
                (None, False), ("True", True), ("true", True), ("TRUE", True),
                ("False", False), ("false", False), ("FALSE", False),
                ("1", True), ("0", False),
            ):
                path.write_text(
                    "action = plot\n" + (
                        f"legacy_unit = {value}\n" if value is not None else ""
                    ), encoding="utf-8",
                )
                self.assertIs(read_input(path).legacy_unit, expected)
            path.write_text("action = plot\nlegacy_unit = yes\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "True or False"):
                read_input(path)

    def test_unmarked_sed_uses_input_scale_without_rescaling(self):
        """Headers from early 3.0 runs override input; unmarked 2.x files need True."""
        with TemporaryDirectory() as directory:
            path = Path(directory) / "old.SED"
            values = np.arange(1.0, 4.0)[:, None]
            for header in ("", "pySED spectral energy density; unit: eV/THz;"):
                np.savetxt(path, values, header=header)
                for fallback in (False, True):
                    loaded, legacy = FileIO._read_sed(path, fallback)
                    np.testing.assert_array_equal(loaded, values)
                    self.assertEqual(legacy, fallback if not header else False)

    def test_partial_sed_rejects_mixed_output_scales(self):
        """Element sums must not combine legacy and physical intensities."""
        with TemporaryDirectory() as directory:
            prefix = str(Path(directory) / "mixed")
            params = SimpleNamespace(
                out_files_name=prefix, plot_partial_SED=True,
                plot_partial_element="Si", plot_partial_dir=None,
            )
            partial_dir = Path(prefix + "_partial_SED")
            partial_dir.mkdir()
            for direction in ("x", "y", "z"):
                np.savetxt(
                    partial_dir / f"mixed.SED_Si_{direction}", [[1.0], [2.0]],
                    header="legacy_unit: true" if direction == "x" else FileIO.SED_HEADER,
                )
            with self.assertRaisesRegex(ValueError, "inconsistent output scales"):
                FileIO.load_data(params)

    def test_sed_reload_preserves_single_frequency_row_orientation(self):
        """A one-row file must stay (frequency, q-point), even for one q."""
        with TemporaryDirectory() as directory:
            params = SimpleNamespace(out_files_name=str(Path(directory) / "short"))
            for num_qpoints in (1, 3):
                with self.subTest(qpoints=num_qpoints):
                    phonons = SimpleNamespace(
                        sed_avg=np.arange(1, num_qpoints + 1)[None, :],
                        freq_fft=np.array([0.0]),
                    )
                    lattice = SimpleNamespace(
                        reduced_qpoints=np.zeros((num_qpoints, 3)),
                        q_distances=np.arange(num_qpoints), q_labels=((0.0, "G"),),
                    )
                    FileIO.write_output(phonons, params, lattice)
                    loaded = FileIO.load_data(params)
                    self.assertEqual(loaded.sed_avg.shape, (1, num_qpoints))
                    self.assertEqual(loaded.freq_fft.shape, (1,))
                    np.testing.assert_array_equal(loaded.sed_avg, phonons.sed_avg)
                    np.testing.assert_array_equal(loaded.freq_fft, phonons.freq_fft)
                    np.testing.assert_array_equal(loaded.qpoints, lattice.reduced_qpoints)

    def test_explicit_compute_runs_when_sed_already_exists(self) -> None:
        with TemporaryDirectory() as directory:
            output = Path(directory) / "result"
            output.with_suffix(".SED").write_text("old", encoding="utf-8")
            params = SimpleNamespace(
                out_files_name=str(output),
                total_num_steps=4,
                output_data_stride=1,
                num_blocks=2,
                trajectory_prefetch=False,
            )
            bz_info = object()
            sed = Mock()
            first_block = SimpleNamespace(
                data={"cells": np.eye(3)[None, :, :]},
            )
            second_block = SimpleNamespace(
                data={"cells": np.eye(3)[None, :, :]},
            )
            raw_blocks = iter([first_block, second_block])
            trajectory = SimpleNamespace(
                info=object(),
                iter_blocks=Mock(return_value=raw_blocks),
            )
            output = io.StringIO()

            with (
                patch("pySED.pipeline.step_prepare_trajectory"),
                patch(
                    "pySED.pipeline.trajectory_block_source",
                    return_value=trajectory,
                ),
                patch(
                    "pySED.pipeline.construct_BZ.BZ_methods",
                    return_value=bz_info,
                ),
                patch(
                    "pySED.pipeline.Phonon.spectral_energy_density",
                    return_value=sed,
                ),
                patch("pySED.pipeline.FileIO.write_output") as write_output,
                patch(
                    "pySED.pipeline.time.perf_counter",
                    side_effect=(10.0, 12.0, 20.0, 23.0),
                ),
                redirect_stdout(output),
            ):
                step_compute_sed(params)

            sed.compute_sed.assert_called_once_with(
                params,
                bz_info,
                first_block,
                raw_blocks,
            )
            write_output.assert_called_once_with(sed, params, bz_info)
            message = output.getvalue()
            self.assertIn("SED compute done (13.0 s)", message)
            self.assertIn("Initial read/setup : 2.0 s", message)
            self.assertIn("SED calculation    : 8.0 s", message)
            self.assertIn("Output/finalization: 3.0 s", message)

    def test_dump_group_supercell_overrides_the_trajectory_cell_for_bz(self) -> None:
        params = SimpleNamespace(
            out_files_name="result",
            total_num_steps=4,
            output_data_stride=1,
            num_blocks=2,
            trajectory_prefetch=False,
            dump_group_supercell=np.eye(3) * 4.0,
        )
        bz_info = object()
        sed = Mock()
        first_block = SimpleNamespace(
            data={"cells": np.eye(3)[None, :, :]},
        )
        second_block = SimpleNamespace(
            data={"cells": np.eye(3)[None, :, :]},
        )
        trajectory = SimpleNamespace(
            info=object(),
            iter_blocks=Mock(return_value=iter([first_block, second_block])),
        )

        with (
            patch("pySED.pipeline.step_prepare_trajectory"),
            patch(
                "pySED.pipeline.trajectory_block_source",
                return_value=trajectory,
            ),
            patch(
                "pySED.pipeline.construct_BZ.BZ_methods",
                return_value=bz_info,
            ) as make_bz,
            patch(
                "pySED.pipeline.Phonon.spectral_energy_density",
                return_value=sed,
            ),
            patch("pySED.pipeline.FileIO.write_output"),
        ):
            step_compute_sed(params)

        make_bz.assert_called_once()
        called_params, called_cell = make_bz.call_args.args
        self.assertIs(called_params, params)
        np.testing.assert_allclose(called_cell, np.eye(3) * 4.0)

    def test_single_trajectory_configuration(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            basis_file = root / "basis.in"
            basis_file.write_text(
                "test structure\n"
                "atoms_ids unitcell_index basis_index mass_types\n"
                "1 1 1 15.9994\n",
                encoding="utf-8",
            )
            input_file = root / "input.in"
            input_file.write_text(
                "action = plot\n"
                "trajectory_file = movie.nc\n"
                "lammps_unit = metal\n"
                "netcdf_compression_level = 0\n"
                "netcdf_batch_size = 64\n"
                "output_partial = 1\n"
                f"basis_lattice_file = {basis_file}\n"
                "plot_partial_SED = O y\n",
                encoding="utf-8",
            )

            params = read_input(str(input_file))

            self.assertEqual(
                params.trajectory_file,
                str((root / "movie.nc").resolve()),
            )
            self.assertEqual(params.lammps_unit, "metal")
            self.assertEqual(params.netcdf_compression_level, 0)
            self.assertEqual(params.netcdf_batch_size, 64)
            self.assertFalse(hasattr(params, "trajectory_path"))
            self.assertFalse(hasattr(params, "source_format"))
            self.assertTrue(params.output_partial)
            self.assertEqual(params.plot_partial_element, "O")
            self.assertEqual(params.plot_partial_dir, "y")
            self.assertFalse(hasattr(params, "plot_partial_type"))
            self.assertEqual(params.trajectory_read_mode, "cache")
            self.assertTrue(params.trajectory_prefetch)
            self.assertEqual(params.device, "cpu")
            self.assertFalse(hasattr(params, "backend"))

            with input_file.open("a", encoding="utf-8") as stream:
                stream.write("device = GPU\n")
            self.assertEqual(read_input(input_file).device, "gpu")

    def test_dump_group_supercell_is_parsed_and_requires_fixed_primitive(self) -> None:
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text(
                "action = plot\n"
                "rescale_prim = 0\n"
                "dump_group_supercell = 4 0 0  0 4 0  0 0 4\n",
                encoding="utf-8",
            )

            params = read_input(str(input_file))
            np.testing.assert_allclose(
                params.dump_group_supercell, np.eye(3) * 4
            )

            input_file.write_text(
                "action = plot\n"
                "dump_group_supercell = 4 0 0  0 4 0  0 0 4\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError, "dump_group_supercell is for a local group trajectory"
            ):
                read_input(str(input_file))

    def test_dump_group_supercell_supports_a_non_diagonal_repetition_matrix(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            basis_file = root / "basis.in"
            basis_file.write_text(
                "test structure\n"
                "atoms_ids unitcell_index basis_index mass_types\n"
                + "".join(
                    f"{atom_id} {atom_id} 1 58.6934\n"
                    for atom_id in range(1, 6145)
                ),
                encoding="utf-8",
            )
            input_file = root / "input.in"
            input_file.write_text(
                "action = compute\n"
                "num_atoms = 6144\n"
                "total_num_steps = 10\n"
                "time_step = 1\n"
                "output_data_stride = 1\n"
                "num_blocks = 1\n"
                f"basis_lattice_file = {basis_file}\n"
                "prim_unitcell = 2.522249888492 0 0 "
                "1.261124944246 2.184332478127 0 "
                "-1.261124944246 -0.728110826042 2.059408410199\n"
                "dump_group_supercell = 80.711996431757 0 0 "
                "40.355998215879 69.898639300061 0 "
                "0 0 12.356450461196\n"
                "rescale_prim = 0\n"
                "num_qpaths = 2\n"
                "q_path_name = GKX\n"
                "q_path = 0 0 0  3/4 3/8 -3/8  1 1/2 -1/2\n",
                encoding="utf-8",
            )

            params = read_input(str(input_file))
            bz = BZ_methods(
                params,
                params.dump_group_supercell,
                box_info=False,
                qpoints_info=False,
            )

            self.assertEqual(bz.num_qpoints, 18)
            np.testing.assert_allclose(
                bz.qpoints[[0, 12, 17]],
                [[0, 0, 0], [1.86832756, 0, 0], [2.49110341, 0, 0]],
                atol=1.0e-8,
            )

            params.dump_group_supercell = np.array(
                [[16, 0, 0], [0, 32, 0], [2, 2, 6]]
            ) @ params.prim_unitcell
            with self.assertRaisesRegex(
                ValueError, "The local cell does not match basis.in"
            ):
                BZ_methods(
                    params,
                    params.dump_group_supercell,
                    box_info=False,
                    qpoints_info=False,
                )

            params.dump_group_supercell = np.eye(3)
            with self.assertRaisesRegex(
                ValueError,
                "The local cell cannot be matched to prim_unitcell",
            ):
                BZ_methods(
                    params,
                    params.dump_group_supercell,
                    box_info=False,
                    qpoints_info=False,
                )

    def test_netcdf_tuning_parameters_are_validated(self) -> None:
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            for setting, message in (
                ("netcdf_compression_level = 10", "between 0 and 9"),
                ("netcdf_batch_size = 0", "must be positive"),
            ):
                input_file.write_text(
                    f"action = plot\n{setting}\n",
                    encoding="utf-8",
                )
                with (
                    self.subTest(setting=setting),
                    self.assertRaisesRegex(ValueError, message),
                ):
                    read_input(str(input_file))

    def test_partial_sed_rejects_unknown_element(self) -> None:
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text(
                "plot_partial_SED = NotAnElement y\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "was not found in the periodic table",
            ):
                read_input(str(input_file))

    def test_partial_sed_element_must_exist_in_structure(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            basis_file = root / "basis.in"
            basis_file.write_text(
                "test structure\n"
                "atoms_ids unitcell_index basis_index mass_types\n"
                "1 1 1 28.0855\n",
                encoding="utf-8",
            )
            input_file = root / "input.in"
            input_file.write_text(
                "action = plot\n"
                f"basis_lattice_file = {basis_file}\n"
                "plot_partial_SED = O y\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "was not found in the input structure",
            ):
                read_input(str(input_file))

    def test_malformed_inputs_are_rejected(self) -> None:
        cases = {
            "missing_equals": (
                "action plot\n",
                "expected 'key = value'",
            ),
            "unknown_parameter": (
                "action = plot\nnum_atom = 10\n",
                "parameter 'num_atom' is not valid",
            ),
            "extra_scalar_value": (
                "action = plot\nnum_atoms = 10 extra\n",
                "requires exactly 1 value",
            ),
            "duplicate_parameter": (
                "action = plot\naction = fit\n",
                "parameter 'action' is repeated",
            ),
            "invalid_boolean": (
                "action = plot\noutput_partial = 2\n",
                "must be 0 or 1",
            ),
            "invalid_device": (
                "action = plot\ndevice = unknown\n",
                "device 'unknown' is not available",
            ),
            "library_name_is_not_a_device": (
                "action = plot\ndevice = numpy\n",
                "device 'numpy' is not available",
            ),
            "old_backend_name": (
                "action = plot\nbackend = numpy\n",
                "parameter 'backend' is not valid",
            ),
            "removed_parameter": (
                "action = plot\ncompress = 1\n",
                "parameter 'compress' was removed",
            ),
            "negative_qpoint_index": (
                "action = plot\nqpoint_slice_index = -1\n",
                "qpoint_slice_index must be non-negative",
            ),
            "invalid_trajectory_read_mode": (
                "action = plot\ntrajectory_read_mode = parallel\n",
                "trajectory_read_mode must be 'cache' or 'direct'",
            ),
        }
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for name, (text, message) in cases.items():
                with self.subTest(name=name):
                    input_file = root / f"{name}.in"
                    input_file.write_text(text, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        read_input(input_file)

    def test_q_path_labels_are_validated(self) -> None:
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text(
                "action = plot\n"
                "num_qpaths = 2\n"
                "q_path_name = GX\n"
                "q_path = 0 0 0  0.5 0 0  0.5 0.5 0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "q_path_name requires 3 labels",
            ):
                read_input(input_file)

    def test_compute_requires_positive_sampling_values(self) -> None:
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text(
                "action = compute\n"
                "time_step = 0\n"
                "output_data_stride = 1\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "time_step must be positive",
            ):
                read_input(input_file)

    @unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
    def test_pipeline_prepares_text_for_compute(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            trajectory_directory = root / "trajectory"
            trajectory_directory.mkdir()
            trajectory = trajectory_directory / "dump.xyz"
            trajectory.write_text(GPUMD_TEXT, encoding="utf-8")
            input_file = root / "input.in"
            input_file.write_text(
                "action = plot\n"
                f"trajectory_file = {trajectory}\n"
                "netcdf_compression_level = 0\n"
                "netcdf_batch_size = 64\n",
                encoding="utf-8",
            )

            params = read_input(str(input_file))
            step_prepare_trajectory(params)

            self.assertEqual(params.source_format, "gpumd_xyz")
            self.assertTrue(Path(params.trajectory_path).is_file())
            self.assertEqual(
                Path(params.trajectory_path).name,
                "dump.xyz.pysed.nc",
            )
            self.assertEqual(Path(params.trajectory_path).parent, root)
            self.assertFalse(
                (trajectory_directory / "dump.xyz.pysed.nc").exists()
            )
            from netCDF4 import Dataset
            with Dataset(params.trajectory_path, "r") as dataset:
                self.assertFalse(
                    dataset.variables["coordinates"].filters()["zlib"]
                )
                self.assertEqual(dataset.pysed_compression_level, 0)

            cached_params = read_input(str(input_file))
            cached_params.trajectory_file = params.trajectory_path
            step_prepare_trajectory(cached_params)
            self.assertEqual(cached_params.source_format, "gpumd_xyz")
            self.assertEqual(
                cached_params.trajectory_path,
                params.trajectory_path,
            )

    def test_pipeline_can_read_text_directly_without_creating_a_cache(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            trajectory = root / "dump.xyz"
            trajectory.write_text(GPUMD_TEXT, encoding="utf-8")
            input_file = root / "input.in"
            input_file.write_text(
                "action = plot\n"
                f"trajectory_file = {trajectory}\n"
                "trajectory_read_mode = direct\n"
                "trajectory_prefetch = 1\n",
                encoding="utf-8",
            )

            params = read_input(str(input_file))
            step_prepare_trajectory(params)

            self.assertEqual(params.source_format, "gpumd_xyz")
            self.assertEqual(params.trajectory_path, str(trajectory))
            self.assertTrue(params.trajectory_prefetch)
            self.assertFalse((root / "dump.xyz.pysed.nc").exists())

    @unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
    def test_pipeline_reads_gpumd_netcdf_directly(self) -> None:
        from netCDF4 import Dataset

        with TemporaryDirectory() as directory:
            root = Path(directory)
            trajectory = root / "movie.nc"
            with Dataset(
                str(trajectory),
                "w",
                format="NETCDF3_64BIT_OFFSET",
            ) as dataset:
                dataset.program = "GPUMD"
                dataset.createDimension("frame", 1)
                dataset.createDimension("atom", 1)
                dataset.createDimension("spatial", 3)
                coordinates = dataset.createVariable(
                    "coordinates",
                    "f4",
                    ("frame", "atom", "spatial"),
                )
                coordinates.units = "angstrom"
                coordinates[:] = [[[0, 0, 0]]]

            input_file = root / "input.in"
            input_file.write_text(
                "action = plot\n"
                f"trajectory_file = {trajectory}\n"
                "trajectory_read_mode = direct\n",
                encoding="utf-8",
            )

            params = read_input(str(input_file))
            step_prepare_trajectory(params)

            self.assertEqual(params.source_format, "gpumd_netcdf")
            self.assertEqual(params.trajectory_path, str(trajectory))
            self.assertFalse((root / "movie.nc.pysed.nc").exists())


if __name__ == "__main__":
    unittest.main()
