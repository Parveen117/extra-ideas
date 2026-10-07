"""CL1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import cl1_clock_curvature_is_mass as y


class ClockTests(unittest.TestCase):
    def test_mass_is_clock_curvature(self):
        rows = y.run()['mass_is_clock_curvature']
        self.assertTrue(all(r['clock_curvature'] == r['rest_defect'] for r in rows))
        self.assertEqual(rows[0]['clock_curvature'], '4/5')

    def test_a_phase_that_is_not_a_turn_is_rejected(self):
        with self.assertRaises(ValueError):
            y.weyl_control(F(3, 5), F(3, 5))

    def test_tower(self):
        rows = y.tower_control()
        self.assertEqual([r['clock_curvature'] for r in rows[:2]], ['4', '2'])
        self.assertEqual(rows[1]['cone_speed'], '0')

    def test_varying_coin_exchange_term(self):
        self.assertEqual(y.varying_coin_control()['constant_coin'], 'no exchange term')

    def test_unconditioned_shift_breaks_the_walk_sum(self):
        real = y.walk_sum

        def broken(coins, psi):
            out = real(coins, psi)
            return {x: (v[0]+psi[x][1], v[1]) for x, v in out.items()}
        with patch.object(y, 'walk_sum', broken):
            with self.assertRaises(ValueError):
                y.varying_coin_control()

    def test_count_rotates_even_and_odd(self):
        self.assertEqual(y.count_control()['count_with_odd'], 'even part')

    def test_history_residue(self):
        rows = y.residue_control()
        self.assertEqual(rows[0]['closes_at'], [4, 8])
        self.assertEqual(rows[1]['closes_at'], [])


if __name__ == '__main__':
    unittest.main()
