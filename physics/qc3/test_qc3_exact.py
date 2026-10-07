"""QC3: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import qc3_even_odd_arrow as y


class EvenOddTests(unittest.TestCase):
    def test_circular_contents(self):
        row = y.circular_control(F(1, 3), range(1, 6))
        self.assertEqual([r['limit'] for r in row['rows']], [1, 4, 9, 16, 25])
        self.assertEqual(row['rows'][0]['defect_ratio'], '1')

    def test_wrong_character_is_rejected(self):
        real = y.cheb_U
        with patch.object(y, 'cheb_U', lambda n, x: real(n, x)+(1 if n == 2 else 0)):
            with self.assertRaises(ValueError):
                y.circular_control(F(1, 3), range(1, 5))

    def test_rates_increase_to_the_square(self):
        self.assertEqual(y.circular_limit(3)['content'], 3)

    def test_su2_casimir_law(self):
        row = y.su2_control(F(1, 4))
        self.assertEqual(row['rate_content_half'], '64/257')
        self.assertEqual(F(row['ratio'])*(1+F(1, 16)**2), F(8, 3))

    def test_a_non_native_turn_is_rejected(self):
        real = y.turn

        def broken(axis, r):
            q = list(real(axis, r))
            q[0] = q[0]+F(1, 100)
            return tuple(q)
        with patch.object(y, 'turn', broken):
            with self.assertRaises(ValueError):
                y.su2_control(F(1, 4))

    def test_one_sided_record_has_no_even_odd_law(self):
        # dropping the -v turn: the 'mean' is U_v itself and E^2 - O^2 = I cannot hold with O = 0
        L = y.left_matrix(y.turn(0, F(1)))
        self.assertNotEqual(y.mat_mul(L, L), y.ident(4))


if __name__ == '__main__':
    unittest.main()
