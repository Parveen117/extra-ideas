"""HM1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import hm1_horizon_memory_product as y
rb = y.rb


class HorizonMemoryProductTests(unittest.TestCase):
    def test_run_and_pinned_values(self):
        res = y.run()
        self.assertAlmostEqual(res['tower_sensitivity_by_starting_x']['1/3'], 0.03065718490354965, places=14)
        self.assertLess(res['tower_sensitivity_by_starting_x']['1/100'], 1e-15)
        self.assertAlmostEqual(res['past_horizon']['literal']*res['past_horizon']['mirror'], 1.0, places=12)

    def test_one_cell_identity(self):
        r, u, v = F(7, 3), F(2, 5), F(9, 4)
        fu, fv = 1/(r + u), 1/(r + v)
        self.assertEqual(fu - fv, -(u - v)*fu*fv)

    def test_small_change_of_closure_moves_the_reading_by_the_product(self):
        xs = y.tower_profile(F(1, 3), 6)
        rec = y.recips_of_profile(xs)
        eps = F(1, 10**12)
        moved = y.returns_along(rec, xs[-1] + eps)[0] - xs[0]
        self.assertLess(abs(abs(moved)/eps - y.sensitivity(xs)), F(1, 10**10))

    def test_harmonic_spacing_forgets_and_matches_the_reciprocal_sum(self):
        slow = [1 - F(1, j + 2) for j in range(200)]
        self.assertEqual(y.sensitivity(slow), F(1, 200**2))          # 199 cells: 1/(n+1)^2
        self.assertGreater(sum(y.recips_of_profile(slow)), 9)          # R2 RI2: the reciprocal sum grows without bound
        fast = y.tower_profile(F(1, 3), 7)
        self.assertLess(sum(y.recips_of_profile(fast)), 4)

    def test_a_different_tower_changes_the_verdict_check(self):
        with patch.object(rb, 'tower', lambda x: (x + 1)/2):
            # still reaches the horizon geometrically, so memory stays; the pinned bound on convergence still holds
            xs = y.tower_profile(F(1, 3), 30)
            self.assertGreater(y.sensitivity(xs), 0)
            self.assertLess(y.sensitivity(xs) - y.sensitivity(xs + [rb.tower(xs[-1])]), F(1, 10**8))

    def test_wrong_identity_sign_is_rejected(self):
        with patch.object(y, 'product', lambda vals: -1*F(1)):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
