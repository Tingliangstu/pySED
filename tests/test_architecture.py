"""Small contracts that keep the SED-only architecture understandable."""

import ast
import os
import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np

from pySED import Lorentz, OutputPaths, Plot_SED, pipeline
from pySED.parser import LEGACY_PARAMETER_HINTS, read_input
from pySED.parameters import PARAMETERS


class ArchitectureTests(unittest.TestCase):
    def _spectrum(self):
        frequency = np.linspace(1.0, 3.0, 401)
        data = SimpleNamespace(
            freq_fft=frequency,
            sed_avg=Lorentz._lorentzian(frequency, 2.0, 1.0, 0.08)[:, None],
            qpoints=np.array([[0.0, 0.0, 0.0]]),
        )
        params = SimpleNamespace(
            qpoint_slice_index=0, lorentz_fit_freq_min=1.0,
            lorentz_fit_freq_max=3.0, fitting_function="lorentz",
            peak_min_significance=4.0, modulate_factor=0,
            initial_guess_hwhm=0.05, peak_max_hwhm=1.0,
            plot_partial_SED=False, if_show_figures=False,
        )
        return data, params

    def test_fitting_returns_results_without_files_or_parameter_mutation(self):
        data, params = self._spectrum()
        original_params = vars(params).copy()
        original_sed = data.sed_avg.copy()
        with TemporaryDirectory() as directory:
            previous = os.getcwd()
            try:
                os.chdir(directory)
                with redirect_stdout(StringIO()):
                    result = Lorentz.lorentz(data, params)
                self.assertEqual(list(Path(directory).iterdir()), [])
            finally:
                os.chdir(previous)
        self.assertEqual(vars(params), original_params)
        np.testing.assert_array_equal(data.sed_avg, original_sed)
        np.testing.assert_allclose(result.popt[0], [2.0, 1.0, 0.08], rtol=1e-5)
        self.assertEqual(len(result.peak_fits), 1)

    def test_plain_plot_does_not_inherit_a_previous_fit_overlay(self):
        data, params = self._spectrum()
        with redirect_stdout(StringIO()):
            result = Lorentz.lorentz(data, params)
        original_params = vars(params).copy()
        with TemporaryDirectory() as directory:
            previous = os.getcwd()
            figures = []
            try:
                os.chdir(directory)
                with (
                    patch("pySED.Plot_SED.plt.savefig") as savefig,
                    patch("pySED.Plot_SED.plt.close"),
                ):
                    Plot_SED.plot_slice(data, params, fit_result=result)
                    figures.append(plt.gcf())
                    self.assertEqual(
                        savefig.call_args.args[0],
                        str(OutputPaths.fitting_qpoint_figure(0)),
                    )
                    self.assertEqual(len(figures[0].axes[0].lines), 2)
                    Plot_SED.plot_slice(data, params)
                    figures.append(plt.gcf())
                    self.assertEqual(savefig.call_args.args[0], "SED-0-qpoint.png")
                    self.assertEqual(len(figures[1].axes[0].lines), 1)
            finally:
                os.chdir(previous)
                for figure in figures:
                    plt.close(figure)
        self.assertEqual(vars(params), original_params)

    def test_pipeline_explicitly_orders_fit_write_and_plot(self):
        data, params = self._spectrum()
        result = SimpleNamespace(peak_fits=[])
        calls = Mock()
        with (
            patch("pySED.pipeline.Lorentz.lorentz", return_value=result) as fit,
            patch("pySED.pipeline.FileIO.write_phonon_lifetime", return_value="fit.dat") as write,
            patch("pySED.pipeline.Plot_SED.plot_slice") as plot,
            redirect_stdout(StringIO()),
        ):
            calls.attach_mock(fit, "fit")
            calls.attach_mock(write, "write")
            calls.attach_mock(plot, "plot")
            self.assertIs(pipeline._fit_one_qpoint(data, params), result)
        self.assertEqual([call[0] for call in calls.mock_calls], ["fit", "write", "plot"])
        plot.assert_called_once_with(data, params, fit_result=result)

    def test_readiness_includes_q_path_metadata(self):
        with TemporaryDirectory() as directory:
            prefix = str(Path(directory) / "sample")
            params = SimpleNamespace(out_files_name=prefix)
            paths = OutputPaths.sed_data_files(prefix)
            for path in paths[:3]:
                path.touch()
            self.assertFalse(pipeline._plot_data_ready(params))
            with self.assertRaisesRegex(FileNotFoundError, "Q_distances_and_labels"):
                pipeline.step_plot_sed(params)
            paths[3].touch()
            self.assertTrue(pipeline._plot_data_ready(params))

    def test_migration_hints_are_not_runtime_parameters(self):
        self.assertFalse(set(PARAMETERS) & set(LEGACY_PARAMETER_HINTS))
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / "input.in"
            input_file.write_text("action = plot\n", encoding="utf-8")
            params = read_input(input_file)
        self.assertFalse(hasattr(params, "with_eigs"))
        self.assertFalse(hasattr(params, "plot_lorentz"))

    def test_fitting_import_does_not_load_plotting_or_file_writers(self):
        subprocess.run(
            [sys.executable, "-c",
             "import sys; from pySED import Lorentz; "
             "assert 'matplotlib' not in sys.modules; "
             "assert 'pySED.FileIO' not in sys.modules; "
             "assert 'pySED.Plot_SED' not in sys.modules"],
            cwd=Path(__file__).resolve().parents[1],
            check=True, capture_output=True, text=True,
        )

    def test_numerics_and_backends_do_not_import_workflow_layers(self):
        root = Path(__file__).resolve().parents[1] / "pySED"
        boundaries = {
            "sed_kernel.py": ("pySED.backends", "pySED.io", "pySED.pipeline",
                              "pySED.FileIO", "pySED.Plot_SED"),
            "Lorentz.py": ("pySED.backends", "pySED.io", "pySED.pipeline",
                           "pySED.FileIO", "pySED.Plot_SED"),
            "backends/cpu.py": ("pySED.Phonon", "pySED.pipeline", "pySED.io"),
            "backends/gpu.py": ("pySED.Phonon", "pySED.pipeline", "pySED.io"),
        }
        for filename, forbidden in boundaries.items():
            tree = ast.parse((root / filename).read_text(encoding="utf-8"))
            imports = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports.extend(alias.name for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imports.extend(node.module + "." + alias.name for alias in node.names)
            for imported in imports:
                self.assertFalse(
                    imported.startswith(forbidden),
                    f"{filename} must not depend on {imported}",
                )


if __name__ == "__main__":
    unittest.main()
