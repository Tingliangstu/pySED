"""Tests for the compact command-line startup panel."""

import io
import unittest
from contextlib import redirect_stdout
from types import SimpleNamespace

from pySED.main import _print_run_summary, _show_help


class MainOutputTests(unittest.TestCase):
    def test_help_describes_the_sed_only_3_0_scope(self) -> None:
        output = io.StringIO()

        with redirect_stdout(output), self.assertRaises(SystemExit):
            _show_help()

        message = output.getvalue()
        self.assertIn("method = sed", message)
        self.assertIn("pySED", message)
        self.assertNotIn("planned extensions", message.lower())

    def test_run_summary_is_clear_and_aligned(self) -> None:
        params = SimpleNamespace(
            method="sed",
            action="thinking",
            backend="numpy",
            max_cores=8,
            out_files_name="SrTiO3",
        )
        output = io.StringIO()

        with redirect_stdout(output):
            _print_run_summary("input.in", params)

        message = output.getvalue()
        self.assertIn("Starting pySED task", message)
        self.assertIn("  Input   : input.in", message)
        self.assertIn("  Method  : SED", message)
        self.assertIn("  Action  : thinking", message)
        self.assertIn(
            "  Backend : NumPy (up to 8 CPU processes)",
            message,
        )
        self.assertIn("  Output  : SrTiO3", message)


if __name__ == "__main__":
    unittest.main()
