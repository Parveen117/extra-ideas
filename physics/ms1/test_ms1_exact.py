"""MS1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import ms1_speed_and_mass as y


class SpeedMassTests(unittest.TestCase):
    def test_state_law(self):
        rows = y.state_control()
        self.assertEqual(rows[0]['speed'], '35/37')
        self.assertEqual(rows[2]['memory'], '0')

    def test_flow_law(self):
        rows = y.flow_control()
        self.assertEqual(rows[2]['speed'], '8/17')
        self.assertEqual(rows[2]['proper_density'], '15/8')
        self.assertTrue(all(r['current'] == '1' for r in rows))

    def test_a_circular_flow_generator_is_rejected(self):
        with patch.object(y, 'K', y.R):
            with self.assertRaises(ValueError):
                y.flow_control()

    def test_combination_law(self):
        rows = y.combination_control()
        self.assertEqual(rows[0]['composite_mass'], '48/25')
        self.assertEqual(rows[0]['mass_defect'], '12/25')
        self.assertTrue(all(F(r['mass_defect']) > 0 for r in rows))

    def test_tower_and_character(self):
        t = y.tower_control(F(40, 41), F(9, 41))
        self.assertEqual(t['speed_reaches_zero_by'], 4)
        self.assertEqual(t['rows'][1]['binding'], '36/1681')

    def test_split_half_turn_breaks_the_character_bound(self):
        with self.assertRaises(ValueError):
            y.tower_control(F(5, 4), F(3, 4))


if __name__ == '__main__':
    unittest.main()
