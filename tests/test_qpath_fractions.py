from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from pySED.My_Parsers import get_parse_input, parse_float_or_fraction
from pySED.construct_BZ import BZPathHelper


class QPathFractionTests(unittest.TestCase):
    def test_explicit_fraction_is_preserved(self):
        self.assertEqual(parse_float_or_fraction('7/15'), Fraction(7, 15))
        self.assertEqual(parse_float_or_fraction('-17/31'), Fraction(-17, 31))

    def test_short_decimal_is_normalized(self):
        self.assertEqual(parse_float_or_fraction('0.33333'), Fraction(1, 3))
        self.assertEqual(parse_float_or_fraction('-0.70000000'), Fraction(-7, 10))

    def test_q_path_before_num_qpaths_remains_exact(self):
        params = self._parse_input(
            """
            q_path = 0 0 0  7/15 -7/15 -0.70000000
            num_qpaths = 1
            q_path_name = 'GQ'
            """
        )

        self.assertEqual(params.q_path.dtype, object)
        self.assertEqual(params.q_path[1, 0], Fraction(7, 15))
        self.assertEqual(params.q_path[1, 1], Fraction(-7, 15))
        self.assertEqual(params.q_path[1, 2], Fraction(-7, 10))

    def test_15_by_15_by_10_path_contains_eight_qpoints(self):
        params = self._parse_input(
            """
            num_qpaths = 1
            q_path_name = 'GQ'
            q_path = 0 0 0  7/15 -7/15 -0.70000000
            """
        )
        helper = BZPathHelper(np.eye(3), np.diag([15.0, 15.0, 10.0]))

        qpoints, distances = helper.build_commensurate_path(
            params.q_path[0],
            params.q_path[1],
        )
        reduced_qpoints = helper.cartesian_to_reduced(qpoints)
        expected = np.array([
            [k / 15, -k / 15, -k / 10]
            for k in range(8)
        ])

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
            input_file = Path(tmpdir) / 'input_SED.in'
            input_file.write_text(contents, encoding='utf-8')
            return get_parse_input(str(input_file))


if __name__ == '__main__':
    unittest.main()
