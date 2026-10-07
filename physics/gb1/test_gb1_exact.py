"""GB1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import gb1_gravity_in_the_block as y


class GravityInTheBlockTests(unittest.TestCase):
    def test_three_factors_of_a_block(self):
        rows = y.factor_control()
        self.assertFalse(rows[0]['reading_changed'])
        self.assertEqual(rows[1]['invariant'], ['80/9', '720'])
        self.assertEqual(rows[2]['invariant'][0], rows[2]['invariant'][1])

    def test_scalar_only_gravity_bends_light_too_little(self):
        row = y.light_control()
        self.assertEqual(F(row['form_of_MO1']), 2*F(row['clock_factor_only']))
        self.assertEqual(row['common_factor'], '0')
        self.assertAlmostEqual(row['at_the_sun_arcsec']['clock_only'], 0.8756, places=3)

    def test_field_block_matches_gr1(self):
        rows = y.field_block_control()
        self.assertEqual((rows[0]['density'], rows[0]['proper_density']), ('5/4', '3/4'))

    def test_a_field_block_with_a_turn_is_not_self_dagger(self):
        real = y.vec
        with patch.object(y, 'vec', lambda v: y.cscale(y.c(0, 1), real(v))):
            with self.assertRaises(ValueError):
                y.field_block_control()

    def test_algebra_and_group_actions(self):
        row = y.action_control()
        self.assertEqual(row['light_change_by_content']['1'], ['-6', '-6', '14', '2'])
        self.assertEqual(row['light_change_by_content']['-1'], ['6', '6', '-14', '-2'])


if __name__ == '__main__':
    unittest.main()
