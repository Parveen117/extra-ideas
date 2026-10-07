"""PR3: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import pr3_accelerated_frame as y

G = F(3, 2)


class AcceleratedFrameTests(unittest.TestCase):
    def test_split_term_present_without_half_density(self):
        self.assertTrue(y.boost_only_control(G)['needed'])

    def test_split_term_absorbed_by_half_density(self):
        self.assertEqual(y.half_density_control(G)['split_term'], '0')

    def test_frame_law_is_d_eta_minus_generator(self):
        for chi in y.TESTS:
            law = y.frame_law(chi, G)
            gen = y.G(chi, G)
            want = (y.rho(y.fadd(y.d_eta(chi[0]), gen[0], -1), -1), y.rho(y.fadd(y.d_eta(chi[1]), gen[1], -1), -1))
            self.assertEqual(law, want)

    def test_wrong_weight_is_rejected(self):
        real = y.mul
        with patch.object(y, 'HALF', F(1, 3)):
            with self.assertRaises(ValueError):
                y.half_density_control(G)
        self.assertIs(y.mul, real)

    def test_square_and_factorization(self):
        self.assertEqual(y.square_control(G)['well_depth'], '1/4 at g rho = 1/2')

    def test_wrong_generator_breaks_the_square(self):
        real = y.rho
        with patch.object(y, 'd_s', y.d_eta):
            with self.assertRaises(ValueError):
                y.square_control(G)
        self.assertIs(y.rho, real)

    def test_both_sheets_share_the_clock(self):
        self.assertEqual(y.antimass_control(G)['clock_factor'], 'rho for both sheets')

    def test_lapse_lemma(self):
        self.assertTrue(y.lapse_control()['lemma'])


if __name__ == '__main__':
    unittest.main()
