"""PR1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import pr1_sector_speed_and_mass as y


class SectorTests(unittest.TestCase):
    def test_walk_laws_for_several_coins(self):
        row = y.walk_control(F(3, 5), F(4, 5))
        self.assertEqual(row['cone_speed'], '3/5')
        self.assertEqual(row['curvature_over_two'], '4/5')
        self.assertEqual(y.walk_control(F(1), F(0))['curvature_over_two'], '0')

    def test_a_coin_that_is_not_a_turn_is_rejected(self):
        with self.assertRaises(ValueError):
            y.walk_control(F(3, 5), F(3, 5))

    def test_unconditioned_shift_breaks_the_walk_identity(self):
        with patch.object(y, 'SHIFT', [[{1: F(1)}, {}], [{}, {1: F(1)}]]):
            with self.assertRaises(ValueError):
                y.walk_control(F(3, 5), F(4, 5))

    def test_r28_coin(self):
        self.assertEqual(y.hadamard_control()['cone_speed_squared'], '1/2')

    def test_velocity_space(self):
        row = y.velocity_control()
        self.assertEqual(row['collinear_sum'], '11/17')
        self.assertEqual(row['rotation_tangent'], '26/127')

    def test_sector_boost_and_foreign_sector(self):
        self.assertEqual(y.sector_control(F(1), F(3, 4), F(5, 4), F(3, 4))['cosh_eta'], '17/8')
        with self.assertRaises(ValueError):
            y.sector_control(F(1), F(1), F(5, 4), F(1))

    def test_a_split_generator_is_not_a_mass_term(self):
        with patch.object(y, 'R', y.K):
            with self.assertRaises(ValueError):
                y.sector_control(F(1), F(3, 4), F(5, 4), F(3, 4))

    def test_two_parities(self):
        row = y.parity_control(F(1), F(4), F(3))
        self.assertEqual(row['frequency_squared'], '25')
        self.assertEqual(row['diffusion_times_two_flip_rate'], '1')


if __name__ == '__main__':
    unittest.main()
