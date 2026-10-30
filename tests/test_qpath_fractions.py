"""Regression tests for exact commensurate q-path construction."""

from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from pySED.parser import read_input
from pySED.construct_BZ import BZPathHelper


class QPathFractionTests(unittest.TestCase):
    def test_explicit_fraction_is_preserved(self):
        params = self._parse_input(
            "num_qpaths = 1\n"
            "q_path_name = GQ\n"
            "q_path = 0 0 0  7/15 -17/31 0\n"
        )

        self.assertEqual(params.q_path.dtype, object)
        self.assertEqual(params.q_path[3], Fraction(7, 15))
        self.assertEqual(params.q_path[4], Fraction(-17, 31))

    def test_short_decimal_is_normalized(self):
        params = self._parse_input(
            "num_qpaths = 1\n"
            "q_path_name = GQ\n"
            "q_path = 0 0 0  0.33333 -0.70000000 0\n"
        )

        self.assertEqual(params.q_path[3], Fraction(1, 3))
        self.assertEqual(params.q_path[4], Fraction(-7, 10))

    def test_q_path_before_num_qpaths_remains_exact(self):
        params = self._parse_input(
            "q_path = 0 0 0  7/15 -7/15 -0.70000000\n"
            "num_qpaths = 1\n"
            "q_path_name = GQ\n"
        )

        reshaped = params.q_path.reshape(-1, 3)
        self.assertEqual(params.q_path.dtype, object)
        self.assertEqual(reshaped[1, 0], Fraction(7, 15))
        self.assertEqual(reshaped[1, 1], Fraction(-7, 15))
        self.assertEqual(reshaped[1, 2], Fraction(-7, 10))

    def test_15_by_15_by_10_path_contains_eight_qpoints(self):
        params = self._parse_input(
            "num_qpaths = 1\n"
            "q_path_name = GQ\n"
            "q_path = 0 0 0  7/15 -7/15 -0.70000000\n"
        )
        helper = BZPathHelper(np.eye(3), np.diag([15.0, 15.0, 10.0]))
        q_path = params.q_path.reshape(-1, 3)

        qpoints, distances = helper.build_commensurate_path(
            q_path[0],
            q_path[1],
        )
        reduced_qpoints = helper.cartesian_to_reduced(qpoints)
        expected = np.array(
            [[k / 15, -k / 15, -k / 10] for k in range(8)]
        )

        np.testing.assert_allclose(reduced_qpoints, expected)
        np.testing.assert_allclose(distances, np.arange(8) / 7)

    def test_programmatic_float_path_remains_supported(self):
        helper = BZPathHelper(np.eye(3), np.diag([15.0, 15.0, 10.0]))

        qpoints, _ = helper.build_commensurate_path(
            np.zeros(3),
            np.array([7 / 15, -7 / 15, -0.7]),
        )

        self.assertEqual(len(qpoints), 8)

    @staticmethod
    def _parse_input(contents):
        with TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.in"
            input_file.write_text(
                "action = plot\n" + contents,
                encoding="utf-8",
            )
            return read_input(input_file)


if __name__ == "__main__":
    unittest.main()
