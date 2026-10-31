"""Command-line help routing and input-file selection behavior."""

from contextlib import contextmanager, redirect_stderr, redirect_stdout
from io import StringIO
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from pySED import main as cli
from pySED.parameters import PARAMETERS


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


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
        self.assertIn("pysed init", output)
        self.assertNotIn("[ Workflow ]", output)
        self.assertLess(len(output.splitlines()), 55)
        read.assert_not_called()
        run.assert_not_called()

    def test_init_creates_a_valid_basic_template_without_running(self):
        with TemporaryDirectory(prefix="pysed init ") as directory:
            calculation = Path(directory) / "SED"
            calculation.mkdir()
            with working_directory(calculation):
                with patch.object(cli, "read_input") as read, patch.object(cli, "run_pipeline") as run:
                    status, output, error = self.invoke("init")
                self.assertEqual(status, 0, error)
                self.assertIn("Created input_SED.in", output)
                self.assertEqual([path.name for path in Path.cwd().iterdir()], ["input_SED.in"])
                read.assert_not_called()
                run.assert_not_called()
                # Parsing needs the user's mapping file, but generation does not.
                structure = Path(directory) / "structure"
                structure.mkdir()
                (structure / "basis.in").touch()
                params = cli.read_input()
                self.assertEqual(params.action, "compute")
                self.assertEqual(params.device, "cpu")
                self.assertGreater(params.num_atoms, 0)
                self.assertEqual(params.q_path.size, (params.num_qpaths + 1) * 3)

    def test_init_preserves_an_existing_input_file_byte_for_byte(self):
        with TemporaryDirectory() as directory, working_directory(directory):
            path = Path("input_SED.in")
            original = b"# user input\r\n\xff\x00\n"
            path.write_bytes(original)
            with patch.object(cli, "run_pipeline") as run:
                status, output, error = self.invoke("init")
            self.assertNotEqual(status, 0)
            self.assertNotIn("Created", output)
            self.assertIn("already exists", error)
            self.assertEqual(path.read_bytes(), original)
            run.assert_not_called()

    def test_init_help_explains_generation_without_writing(self):
        with TemporaryDirectory() as directory, working_directory(directory):
            with patch.object(cli, "run_pipeline") as run:
                for flag in ("-h", "--help"):
                    with self.subTest(flag=flag):
                        status, output, error = self.invoke("init", flag)
                        self.assertEqual(status, 0, error)
                        self.assertIn("pysed init", output)
                        self.assertIn("never overwritten", output)
                        self.assertEqual(list(Path.cwd().iterdir()), [])
            run.assert_not_called()

    def test_init_extra_arguments_do_not_create_a_file(self):
        with TemporaryDirectory() as directory, working_directory(directory):
            status, output, error = self.invoke("init", "custom.in")
            self.assertNotEqual(status, 0)
            self.assertIn("pysed init", error)
            self.assertEqual(list(Path.cwd().iterdir()), [])

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
        for path in ("action", "device", "./init", "folder with spaces/custom.in"):
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
