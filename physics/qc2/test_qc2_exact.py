"""QC2 exact part: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import qc2_unit_of_return as y

H = [[F(5, 2), F(-3, 4)], [F(-3, 4), F(7, 3)]]


class ExactTests(unittest.TestCase):
    def test_intensive_offset_is_the_scaled_derivative(self):
        row = y.derivative_control(H, F(2, 3), 5)
        self.assertEqual(row['identities'], 42)
        self.assertEqual(row['mean_pairing'], '4/3')

    def test_a_non_conjugate_offset_is_rejected(self):
        real = y.pmul
        with patch.object(y, 'pmul', lambda p, q: real({k: v*2 for k, v in p.items()}, q)):
            with self.assertRaises(ValueError):
                y.derivative_control(H, F(1), 3)

    def test_number_content_and_norms(self):
        row = y.wick_control(F(2, 3), 4)
        self.assertEqual(row['contents_by_number']['3'], [-3, -1, 1, 3])
        self.assertEqual(row['unit_per_quantum'], '4/3')

    def test_wrong_unit_breaks_the_wick_structure(self):
        real = y.wick
        with patch.object(y, 'wick', lambda a, b, k: real(a, b, k+F(1, 5))):
            with self.assertRaises(ValueError):
                y.wick_control(F(1), 3)

    def test_isotropic_moments(self):
        self.assertEqual(y.zmoment(2, 2, F(1, 2)), (F(2), F(0)))
        self.assertEqual(y.zmoment(3, 1, F(1, 2)), (F(0), F(0)))


if __name__ == '__main__':
    unittest.main()
