"""DG1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import dg1_the_diagonal as y


class DiagonalTests(unittest.TestCase):
    def test_diagonal(self):
        res = y.run()
        self.assertEqual(res['orbit_energy_squared'], '1')
        self.assertEqual(res['orbit_speed_squared'], '1/2')
        self.assertEqual(res['outside'], 'bound')
        self.assertEqual(res['inside'], 'unbound')

    def test_a_different_orbit_law_moves_the_radius(self):
        with patch.object(y, 'orbit_speed2', lambda x: x/2):
            with self.assertRaises(ValueError):
                y.run()

    def test_last_circular_radius(self):
        with self.assertRaises(ZeroDivisionError):
            y.energy2(F(2, 3))


if __name__ == '__main__':
    unittest.main()
