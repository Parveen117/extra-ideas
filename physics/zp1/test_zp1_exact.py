"""ZP1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import zp1_floor_and_force as y


class FloorAndForceTests(unittest.TestCase):
    def test_all(self):
        res = y.run()
        self.assertEqual(res['plate_energy_coefficient'], '-1/720')
        self.assertEqual(res['plate_force_coefficient'], '-1/240')
        self.assertAlmostEqual(res['pressure_Pa']['100 nm'], 13.0, places=1)

    def test_constant_is_sector_blind(self):
        self.assertEqual(y.scale_sector()[2], y.phase_sector()[2])
        self.assertEqual(y.scale_sector()[0], -y.phase_sector()[0])

    def test_higher_series_terms(self):
        self.assertEqual(y.bernoulli_generating()[6], F(1, 30240))
        self.assertEqual(y.scale_sector()[6], F(-1, 6048))

    def test_a_wrong_series_is_rejected(self):
        with patch.object(y, 'bernoulli_generating', lambda: [F(1), F(-1, 2), F(1, 12), F(0), F(-1, 700)] + [F(0)]*5):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
