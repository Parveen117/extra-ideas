"""MC1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import mc1_minimum_cost as y


class MinimumCostTests(unittest.TestCase):
    def test_power_law_in_each_dimension(self):
        for d in (1, 3, 4, 5):
            self.assertEqual(y.shell_control(F(2), d, 6, F(1), F(0))['exponent'], -(d-2))
        self.assertEqual(y.shell_control(F(2), 3, 6, F(1), F(0))['profile'][:3], ['1', '31/63', '5/21'])

    def test_two_dimensions_is_linear_in_scale(self):
        self.assertEqual(y.shell_control(F(2), 2, 6, F(1), F(0))['profile'][:3], ['1', '5/6', '2/3'])

    def test_wrong_weights_do_not_give_the_power_law(self):
        with patch.object(y, 'weights', lambda b, d, n: [b**((d-1)*k) for k in range(n)]):
            with self.assertRaises(ValueError):
                y.shell_control(F(2), 3, 6, F(1), F(0))

    def test_a_non_minimizer_fails_stationarity(self):
        w = y.weights(F(2), 3, 4)
        f, _ = y.minimizer(w, F(1), F(0))
        f[2] += F(1, 10)
        with self.assertRaises(ValueError):
            y.tridiagonal_check(f, w)

    def test_gravity_case(self):
        g = y.gravity_control()
        self.assertEqual({r['flux'] for r in g['chains']}, {'-8/25'})
        self.assertEqual(g['least_cost_N_vs_memory']['difference_in_N_squared'], '81/2500')

    def test_four_dimensional_tail(self):
        self.assertEqual(y.tail_control()['exponent'], -2)


if __name__ == '__main__':
    unittest.main()
