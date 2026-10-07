"""PH3 exact part: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import ph3_two_readings_transport as y


class ExactTests(unittest.TestCase):
    def test_share_law_and_flat_readings(self):
        row = y.exact_controls()
        self.assertTrue(row['pure_readings_return_exactly'])
        self.assertEqual(len(row['points']), 3)

    def test_a_non_derivative_field_is_rejected(self):
        real = y.field

        def broken(x, y_):
            H, Hx, Hy, Hxy = real(x, y_)
            return H, Hx, Hy, [[F(2), F(5)], [F(4), F(0)]]
        with patch.object(y, 'field', broken):
            with self.assertRaises(ValueError):
                y.exact_controls()


if __name__ == '__main__':
    unittest.main()
