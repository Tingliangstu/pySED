"""Regression tests for exact commensurate q-path construction."""

from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from pySED.parser import read_input
from pySED.sed.construct_BZ import BZPathHelper


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
        cases = [
            ([15, 15, 10], [7 / 15, -7 / 15, -0.7],
             [[k / 15, -k / 15, -k / 10] for k in range(8)]),
            ([300, 3, 3], [0.33333, 0, 0],
             [[k / 300, 0, 0] for k in range(101)]),
        ]
        for repeats, end, expected in cases:
            with self.subTest(repeats=repeats):
                helper = BZPathHelper(np.eye(3), np.diag(repeats))
                qpoints, _ = helper.build_commensurate_path(np.zeros(3), end)
                np.testing.assert_allclose(
                    helper.cartesian_to_reduced(qpoints), expected, atol=1e-14,
                )
                self._assert_periodic(helper, qpoints)

    def test_non_diagonal_paths_include_only_lattice_intersections(self):
        primitive = np.array([[2, 0, 0], [0.5, 3, 0], [0.2, 0.3, 4]])
        repeats = np.array([[4, 2, 0], [0, 6, 0], [0, 0, 2]])
        helper = BZPathHelper(primitive, repeats @ primitive)
        cases = [
            ([Fraction(1, 12), 0, 0], [Fraction(1, 12), 1, 0], [1 / 3, 5 / 6]),
            ([Fraction(1, 12), 1, 0], [Fraction(1, 12), 0, 0], [1 / 6, 2 / 3]),
            ([Fraction(1, 10), 0, 0], [Fraction(1, 10), 1, 0], []),
        ]
        for start, end, expected in cases:
            with self.subTest(start=start, end=end):
                qpoints, fractions = helper.build_commensurate_path(start, end)
                np.testing.assert_allclose(fractions, expected)
                start, end = np.asarray(start, float), np.asarray(end, float)
                reduced = start + np.asarray(expected)[:, None] * (end - start)
                np.testing.assert_allclose(
                    helper.cartesian_to_reduced(qpoints), reduced, atol=1e-14,
                )
                self._assert_periodic(helper, qpoints)

    def test_constant_and_short_paths_respect_periodicity(self):
        cases = [
            (4, 0.25, 0.25, [0.0]),
            (4, 0.1, 0.1, []),
            (1_000_000, 0.123456, 0.123457, [0.0, 1.0]),
        ]
        for repeats, start, end, expected in cases:
            with self.subTest(repeats=repeats, start=start, end=end):
                helper = BZPathHelper(np.eye(3), np.diag([repeats, 1, 1]))
                qpoints, fractions = helper.build_commensurate_path(
                    [start, 0, 0], [end, 0, 0],
                )
                np.testing.assert_allclose(fractions, expected)
                self.assertEqual(qpoints.shape, (len(expected), 3))
                self._assert_periodic(helper, qpoints)

    def test_invalid_cells_and_path_endpoints_are_rejected(self):
        for cell in [np.eye(2), np.zeros((3, 3)), np.full((3, 3), np.nan)]:
            for primitive, supercell in [(cell, np.eye(3)), (np.eye(3), cell)]:
                with self.subTest(primitive=primitive, supercell=supercell):
                    with self.assertRaises(ValueError):
                        BZPathHelper(primitive, supercell)
        for repeats in [[10000.05, 4, 4], [1e-5, 1, 1]]:
            with self.subTest(repeats=repeats):
                with self.assertRaises(ValueError):
                    BZPathHelper(np.eye(3), np.diag(repeats))

        helper = BZPathHelper(np.eye(3), np.eye(3))
        for endpoint in [[0, 0], [np.nan, 0, 0], [np.inf, 0, 0]]:
            with self.subTest(endpoint=endpoint):
                with self.assertRaises(ValueError):
                    helper.build_commensurate_path([0, 0, 0], endpoint)

        helper = BZPathHelper(np.eye(3), np.diag([4.00005, 4, 4]))
        with self.assertRaisesRegex(ValueError, "periodicity"):
            helper.build_commensurate_path([0, 0, 0], [3, 0, 0])

    def _assert_periodic(self, helper, qpoints):
        turns = qpoints @ helper.supercell.T / (2 * np.pi)
        np.testing.assert_allclose(turns, np.rint(turns), rtol=0, atol=1e-8)

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
