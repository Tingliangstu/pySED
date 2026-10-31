"""Command-line help routing and input-file selection behavior."""

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from pySED import main as cli
from pySED.parameters import PARAMETERS


class CommandLineTests(unittest.TestCase):
    def invoke(self, *arguments):
        stdout, stderr = StringIO(), StringIO()
        with (
            patch.object(sys, "argv", ["pysed", *arguments]),
            redirect_stdout(stdout),
            redirect_stderr(stderr),
        ):
            try:
                status = cli.main()
            except SystemExit as error:
                status = error.code
        return status, stdout.getvalue(), stderr.getvalue()

    def test_global_help_is_short_and_points_to_parameter_help(self):
        with patch.object(cli, "read_input") as read, patch.object(cli, "run_pipeline") as run:
            status, output, error = self.invoke("-h")
        self.assertEqual(status, 0)
        self.assertFalse(error)
        self.assertIn("v-3.0.0", output)
        self.assertIn("https://pysed.readthedocs.io/", output)
        self.assertIn("pysed <parameter> -h", output)
        self.assertIn("pysed total_num_steps -h", output)
        self.assertNotIn("[ Workflow ]", output)
        self.assertLess(len(output.splitlines()), 55)
        read.assert_not_called()
        run.assert_not_called()

    def test_every_current_parameter_has_offline_help_without_running(self):
        with (
            patch.object(cli, "read_input") as read,
            patch.object(cli, "run_pipeline") as run,
            patch.object(cli, "_print_banner") as banner,
        ):
            for name in PARAMETERS:
                with self.subTest(parameter=name):
                    status, output, error = self.invoke(name, "-h")
                    self.assertEqual(status, 0, error)
                    self.assertFalse(error)
                    self.assertIn(f"PARAMETER: {name}", output)
                    self.assertIn("SYNTAX:", output)
                    self.assertIn("DEFAULT:", output)
                    self.assertNotIn(":doc:", output)
                    self.assertNotIn(":math:", output)
                    self.assertLessEqual(max(map(len, output.splitlines())), 96)
        read.assert_not_called()
        run.assert_not_called()
        banner.assert_not_called()

    def test_double_dash_help_and_mixed_case_parameter_names(self):
        for name in ("plot_partial_SED", "PLOT_PARTIAL_SED"):
            with self.subTest(name=name):
                status, output, error = self.invoke(name, "--help")
                self.assertEqual(status, 0, error)
                self.assertIn("PARAMETER: plot_partial_SED", output)

    def test_unknown_parameter_is_an_error_not_an_input_file(self):
        with patch.object(cli, "read_input") as read, patch.object(cli, "run_pipeline") as run:
            status, output, error = self.invoke("total_num_step", "-h")
        self.assertNotEqual(status, 0)
        self.assertFalse(output)
        self.assertIn("total_num_steps", error)
        read.assert_not_called()
        run.assert_not_called()

    def test_removed_names_are_not_help_aliases(self):
        for name in ("backend", "method"):
            with self.subTest(name=name):
                status, output, error = self.invoke(name, "-h")
                self.assertNotEqual(status, 0)
                self.assertFalse(output)
                self.assertIn(name, error)

    def test_single_argument_is_still_an_input_filename(self):
        for path in ("action", "device", "folder with spaces/custom.in"):
            with self.subTest(path=path), patch.object(sys, "argv", ["pysed", path]):
                self.assertEqual(cli._detect_input_file(), path)
        with patch.object(sys, "argv", ["pysed"]):
            self.assertEqual(cli._detect_input_file(), "input_SED.in")

    def test_malformed_invocation_does_not_start_calculation(self):
        with patch.object(cli, "run_pipeline") as run:
            status, output, error = self.invoke("action", "compute")
        self.assertNotEqual(status, 0)
        self.assertIn("pysed <parameter> -h", error)
        run.assert_not_called()

    def test_run_summary_uses_public_device_values(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "custom.in"
            for device, label in (("cpu", "CPU"), ("gpu", "GPU")):
                with self.subTest(device=device):
                    path.write_text(f"action = plot\ndevice = {device}\n", encoding="utf-8")
                    with patch.object(cli, "run_pipeline") as run:
                        status, output, error = self.invoke(str(path))
                    self.assertEqual(status, 0, error)
                    self.assertIn(f"Device  : {label}", output)
                    self.assertNotIn("Backend", output)
                    params = run.call_args.args[0]
                    self.assertEqual(params.device, device)
                    self.assertFalse(hasattr(params, "backend"))


if __name__ == "__main__":
    unittest.main()
