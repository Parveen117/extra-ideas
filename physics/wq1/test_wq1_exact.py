"""WQ1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import wq1_wave_count as y


class WaveCountTests(unittest.TestCase):
    def test_whole_areas_only(self):
        self.assertEqual(len(y.orbit_control()), 2)
        self.assertTrue(y.silent(1, y.action(F(6), F(3)), F(1, 2)))
        self.assertFalse(y.silent(1, y.action(F(1), F(3)), F(1, 2)))

    def test_ladder(self):
        res = y.ladder_control()
        self.assertEqual(res['lowest'], '7/10')
        self.assertEqual(res['step'], '7/5')

    def test_without_the_half_the_ladder_fails(self):
        def no_half(p, kappa, omega):
            zd = y.raise_(y.lower(p, kappa))
            return [omega*a for a in zd]
        with patch.object(y, 'hamiltonian', no_half):
            with self.assertRaises(ValueError):
                y.ladder_control()

    def test_count_is_frame_independent(self):
        self.assertEqual(y.frame_control()[0], ('2', '12', '9'))

    def test_superposition_and_numbers(self):
        self.assertEqual(y.superposition_control()['total_energy'], '33/4')
        self.assertAlmostEqual(y.illustration()['green_light_500nm_J']*1e19, 3.973, places=3)


if __name__ == '__main__':
    unittest.main()
