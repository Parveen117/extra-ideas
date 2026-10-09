"""Tests for GC1.  python -m unittest test_gc1"""
from fractions import Fraction as F
from math import exp, factorial, pi, sqrt
import json
import os
import unittest

import gc1_core_gap_count as gc


class TestGC1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = gc.run()

    def test_all(self):
        bad = [k for k, v in self.out['detail'].items() if not v]
        self.assertEqual(bad, [])
        self.assertEqual(self.out['checks'], 42)

    def test_series_by_hand(self):
        # v = s - lam s^3/6 + s^4/12 + lam^2 s^5/120 - ...
        lam = F(7, 3)
        a, tail, slope = gc.series(lam, 4)
        c = [F(x, factorial(n)*lam.denominator**n) for n, x in enumerate(a[:6])]
        self.assertEqual(c, [0, 1, 0, -lam/6, F(1, 12), lam*lam/120])
        self.assertLess(tail, F(1, 10**20))

    def test_solution_against_a_plain_integration(self):
        lam = float(gc.E0)

        def f(s, y):
            return (y[1], (s - lam)*y[0])
        h, y, s = 1e-3, (0.0, 1.0), 0.0
        values = gc.grid_values(gc.E0, 6, 8)
        for k in range(1, 6001):
            k1 = f(s, y)
            k2 = f(s + h/2, (y[0] + h/2*k1[0], y[1] + h/2*k1[1]))
            k3 = f(s + h/2, (y[0] + h/2*k2[0], y[1] + h/2*k2[1]))
            k4 = f(s + h, (y[0] + h*k3[0], y[1] + h*k3[1]))
            y = (y[0] + h/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]), y[1] + h/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]))
            s += h
            if k % 1000 == 0:
                (vl, vu), (dl, du) = values[8*k//1000]
                self.assertAlmostEqual(y[0], float(vl), places=7)
                self.assertAlmostEqual(y[1], float(dl), places=7)
                self.assertLess(vu - vl, F(1, 10**11))

    def test_levels_bracket_the_older_values(self):
        # in older terms the two levels are minus the first two zeros of the Airy function: 2.338107.., 4.087949..
        self.assertTrue(gc.level_one(F(2338, 1000)))
        self.assertFalse(gc.level_one(F(23382, 10000)))
        ok, node = gc.level_two(F(40878, 10000))
        self.assertTrue(ok)
        self.assertLess(abs(float(node[0]) - (4.0878 - 2.338107)), 0.01)
        self.assertFalse(gc.level_two(F(4088, 1000))[0])
        # the bound on a step needs |s - lam| h^2 < 1 : a grid of whole steps is refused
        with self.assertRaises(ValueError):
            gc.level_one(gc.E0, step=1)
        # the certificate does not rest on the fine grid
        self.assertTrue(gc.level_one(gc.E0, step=16))

    def test_zero_average_local_rate(self):
        k = sqrt(2)/3

        def v(s):
            return s*s*exp(-k*s**1.5)
        low = 1e9
        for i in range(1, 400):
            s, h = i/40, 1e-4
            rate = -(v(s + h) - 2*v(s) + v(s - h))/(h*h)/v(s) + 2/(s*s) + s
            self.assertAlmostEqual(rate, 27*k/4/sqrt(s) + (1 - 9*k*k/4)*s, places=4)
            low = min(low, rate)
        self.assertGreater(low, float(gc.P0))
        self.assertLess(low, 3.25)

    def test_floor_of_lowest_on_a_small_case(self):
        # rates 1, 3, 5 ; a reading with weights 99/100, 1/200, 1/200
        w = (F(99, 100), F(1, 200), F(1, 200))
        eta = sum(x*r for x, r in zip(w, (1, 3, 5)))
        square = sum(x*r*r for x, r in zip(w, (1, 3, 5)))
        low = gc.floor_of_lowest(eta, square, F(3))
        self.assertLessEqual(low, 1)
        self.assertGreater(low, F(9, 10))
        self.assertEqual(gc.floor_of_lowest(F(1), F(1), F(3)), 1)
        with self.assertRaises(AssertionError):
            gc.floor_of_lowest(eta, square, eta)

    def test_comparison_on_gaussian_readings(self):
        # mean of h is above the mean of the d one-turn operators in every Gaussian reading exp(-w |C|^2/2)
        for d in (3, 4, 5):
            beta = sqrt((d - 1)/2)
            for i in range(1, 60):
                w = i/10
                mean_h = 3*d*w/4 + 3*d*(d - 1)/(4*w*w)
                mean_k = 3*d*w/8 + beta*d*2/sqrt(pi*w)
                self.assertGreater(mean_h, mean_k)

    def test_numbers(self):
        three, four = self.num['d = 3'], self.num['d = 4']
        self.assertEqual(three['lowest_rate'], [5.1865, 5.1868])
        self.assertEqual(three['gap'][0], 0.3342)
        self.assertEqual(four['gap'][0], 0.0031)
        # older values (small-volume, unit differing by 2^(1/3)): 4.1167 for the lowest rate of three turns
        self.assertTrue(three['lowest_rate'][0] <= 4.1167*2**(1/3) <= three['lowest_rate'][1] + 1e-4)
        self.assertEqual(self.num['rho']['d >= 4'], [1.3857, 1.4175])
        with open(os.path.join(os.path.dirname(os.path.abspath(gc.__file__)), 'GC1_RESULT.json')) as f:
            self.assertTrue(json.load(f)['all_pass'])


if __name__ == '__main__':
    unittest.main()
