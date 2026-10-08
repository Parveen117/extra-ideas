"""EN1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import en1_native_energy_bound as y


class NativeEnergyBoundTests(unittest.TestCase):
    def test_run_and_pinned_values(self):
        res = y.run()
        self.assertEqual(res['boost_bound_along_tower'],
                         [('1/3', '2'), ('3/5', '4'), ('15/17', '16'), ('255/257', '256')])
        self.assertEqual(res['turn'], {'gain': '1', 'mass_bound': '49/25'})

    def test_other_seed(self):
        y.run(seed=2026, trials=120)

    def test_weighted_cut_square_identity(self):
        xs, ds = [F(3), F(-2, 5), F(7, 3)], [F(1, 2), F(4), F(5, 7)]
        lhs = sum(ds)*sum(x*x/d for x, d in zip(xs, ds)) - sum(xs)**2
        rhs = sum((ds[j]*xs[i] - ds[i]*xs[j])**2/(ds[i]*ds[j]) for i in range(3) for j in range(i + 1, 3))
        self.assertEqual(lhs, rhs)

    def test_bound_is_attained_only_by_the_co_moving_null_reading(self):
        b = F(3, 5)
        k2 = (1 + b)/(1 - b)
        n2 = 1 - b*b
        self.assertEqual(y.gain(b, (F(2), F(2)), 1)/n2, k2)
        self.assertEqual(y.gain(b, (F(2), F(-2)), 1)/n2, 1/k2)
        self.assertLess(y.gain(b, (F(2), F(1)), 1)/n2, k2)

    def test_own_square_of_the_boost_sector_can_be_created_from_zero(self):
        self.assertEqual(y.N_own((F(1), F(1)), 1), 0)
        self.assertEqual(y.N_own((F(2), F(0)), 1), 4)

    def test_a_larger_gauge_claim_is_rejected(self):
        with patch.object(y, 'D', lambda z: max(abs(z[0]), abs(z[1]))):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
