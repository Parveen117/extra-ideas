"""GR1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import gr1_memory_dictionary as y


class MemoryDictionaryTests(unittest.TestCase):
    def test_local_laws_and_no_memory_no_change(self):
        rows = y.local_control()
        self.assertEqual(rows[0]['dN'], '-480/289')
        self.assertEqual(rows[-1]['dN'], '0')

    def test_a_circular_static_generator_is_rejected(self):
        real = y.ms1.state
        # with psi -> g R psi instead of g K psi the current would not be constant
        def fake_state(psi):
            n, j, s = real(psi)
            return n, j+1, s
        with patch.object(y.ms1, 'state', fake_state):
            with self.assertRaises(ValueError):
                y.local_control()

    def test_dictionary(self):
        rows = y.dictionary_control()
        self.assertEqual(rows[0]['N'], '4/5')
        self.assertEqual(rows[0]['share'], '9/10')
        self.assertEqual(rows[-1]['share'], '1')

    def test_horizon_is_the_balanced_share_and_is_not_reached(self):
        row = y.horizon_control()
        self.assertGreater(F(row['last_share']), F(1, 2))
        self.assertGreater(F(row['last_memory_gap']), 0)

    def test_illustration_reproduces_surface_gravity(self):
        ill = y.illustration()
        self.assertAlmostEqual(ill['Earth surface']['acceleration_m_s2'], 9.82, places=2)
        self.assertAlmostEqual(ill['two_G_over_c2_m_per_kg']*1e27, 1.4852, places=3)


if __name__ == '__main__':
    unittest.main()
