"""TL1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import tl1_temperature_and_the_clock as y


class TemperatureClockTests(unittest.TestCase):
    def test_transport_of_a_turn(self):
        self.assertEqual(y.redshift_control()[0], ('3/5', '4/5', '7', '21/4'))

    def test_units_follow_the_clock(self):
        rows = y.unit_control()
        self.assertIn(('3/5', '4/5', '4', '3', True), rows)
        self.assertIn(('3/5', '4/5', '1', '1', False), rows)

    def test_equal_units_at_one_place_is_the_special_case(self):
        with patch.object(y, 'clock', lambda a, b: F(1)):
            with self.assertRaises(ValueError):
                y.unit_control()

    def test_no_net_flow_iff_equal_weight_per_static_quantum(self):
        rows = y.equilibrium_control()
        self.assertEqual(rows[0][2], '0')
        self.assertEqual(rows[2][2], '3751/22932')
        self.assertTrue(rows[3][2].startswith('-'))

    def test_weight_and_illustration(self):
        y.weight_control()
        ill = y.illustration()
        self.assertAlmostEqual(ill['earth_fraction_per_metre']*1e16, 1.0911, places=3)
        self.assertAlmostEqual(ill['sun_surface_to_far_fraction']*1e6, 2.1225, places=3)


if __name__ == '__main__':
    unittest.main()
