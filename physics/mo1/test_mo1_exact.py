"""MO1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import mo1_motion as y


class MotionTests(unittest.TestCase):
    def test_count_form_and_the_two_special_readings(self):
        rows = y.form_control()
        self.assertEqual(rows[0]['count_at_rest'], '4/5')
        self.assertEqual(rows[0]['count_with_the_frame'], '1')

    def test_wrong_clock_factor_is_rejected(self):
        with patch.object(y, 'PYTH', [(F(3, 5), F(3, 5))]):
            with self.assertRaises(ValueError):
                y.form_control()

    def test_radial_law(self):
        self.assertEqual(y.radial_law_control(), 12)

    def test_circles(self):
        rows = y.circular_control()
        self.assertEqual(rows[0]['count_factor_squared'], '5/6')
        self.assertTrue(rows[-1]['light_like'])

    def test_orbit_equation_and_advance(self):
        row = y.orbit_control()
        self.assertEqual(row['near_circular'][0]['advance_per_turn_over_2pi'], '1/4')

    def test_light(self):
        rows = y.light_control()
        self.assertEqual(rows[0]['deflection'], '1/50')
        self.assertEqual(F(rows[0]['deflection']), 2*F(rows[0]['clock_only_value']))

    def test_illustration_matches_the_known_values(self):
        ill = y.illustration()
        self.assertAlmostEqual(ill['mercury_advance_arcsec_per_century'], 42.98, places=1)
        self.assertAlmostEqual(ill['light_at_solar_limb_arcsec'], 1.75, places=2)


if __name__ == '__main__':
    unittest.main()
