"""Tests for DR1.  python -m unittest test_dr1"""
from fractions import Fraction as F
from math import exp, gamma, sqrt
import json
import os
import unittest

import dr1_core_by_rows as dr


class TestDR1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = dr.run()

    def test_all(self):
        bad = [k for k, v in self.out['detail'].items() if not v]
        self.assertEqual(bad, [])
        self.assertEqual(self.out['checks'], 20)

    def test_series_by_hand(self):
        lam = F(7, 2)
        a, tail, slope = dr.series(lam, 2, 6)                 # three directions: w = v/s of GC1
        self.assertEqual(a[:4], [1, 0, -lam/6, F(1, 12)])
        a, tail, slope = dr.series(lam, 4, 6)                 # five directions: w'' + (4/s) w' = (s - lam) w
        self.assertEqual(a[:4], [1, 0, -lam/10, F(1, 18)])
        self.assertLess(tail, F(1, 10**19))

    def test_series_solves_the_equation_in_floats(self):
        lam, big_p = F(29, 10), 3
        a, tail, slope = dr.series(lam, big_p, 8)
        c = [float(x) for x in a]

        def w(s, order=0):
            if order == 0:
                return sum(x*s**n for n, x in enumerate(c))
            if order == 1:
                return sum(n*x*s**(n - 1) for n, x in enumerate(c) if n)
            return sum(n*(n - 1)*x*s**(n - 2) for n, x in enumerate(c) if n > 1)
        for s in (0.3, 1.1, 2.7, 4.0):
            self.assertAlmostEqual(w(s, 2) + big_p/s*w(s, 1), (s - float(lam))*w(s), places=9)

    def test_float_levels_against_older_values(self):
        one, two = dr.float_levels(2)                          # three directions: 2.338107, 4.087949
        self.assertAlmostEqual(one, 2.338107, places=3)
        self.assertAlmostEqual(two, 4.087949, places=3)
        self.assertAlmostEqual(dr.float_levels(4)[0], 3.361255, places=3)      # l = 1 in three directions
        rows = self.num['floors of one row: m, first level, second level (l = 0)']
        for m, f1, f2 in rows[1:]:
            lv = dr.float_levels(m - 1)
            self.assertTrue(lv[0] - 0.004 < f1 < lv[0])
            if f2 is not None:
                self.assertTrue(lv[1] - 0.004 < f2 < lv[1])

    def test_certificates_do_not_rest_on_the_grid(self):
        self.assertTrue(dr.level_one(F(287, 100), 3, 9, step=32))
        self.assertTrue(dr.level_two(F(4491, 1000), 3, 11, step=128))
        with self.assertRaises(ValueError):
            dr.level_one(F(287, 100), 3, 9, step=1)

    def test_closed_floor_local_rate(self):
        k = sqrt(2)/3
        for p in (1.5, 2.5, 4.0):
            def v(s):
                return s**p*exp(-k*s**1.5)
            low = 1e9
            for i in range(10, 600):
                s, h = i/40, 1e-4
                rate = -(v(s + h) - 2*v(s) + v(s - h))/(h*h)/v(s) + p*(p - 1)/(s*s) + s
                self.assertLess(abs(rate - (3*k*(4*p + 1)/4/sqrt(s) + (1 - 9*k*k/4)*s)), 1e-4*rate)
                low = min(low, rate)
            self.assertGreater(low, 0.75*(4*p + 1)**(2/3) - 1e-6)
            self.assertLess(low, 0.75*(4*p + 1)**(2/3) + 1e-3)

    def test_comparison_by_rows_on_gaussian_readings(self):
        for d in range(3, 12):
            for i in range(1, 60):
                w = i/10
                mean_h = 3*d*w/4 + 3*d*(d - 1)/(4*w*w)
                mean_k = 3*(d*w/8 + (d - 1)/2*gamma((d + 1)/2)/(gamma(d/2)*sqrt(w)))
                self.assertGreater(mean_h, mean_k)

    def test_numbers(self):
        rows = {r[0]: r for r in self.num['d, line, lowest rate (low, high), gap floor, rho^3 (low, high)']}
        self.assertEqual(rows[3][1:5], [5.521, 5.1865, 5.1868, 0.3342])
        self.assertGreater(rows[4][4], 0.44)
        self.assertTrue(all(rows[d][4] < rows[d + 1][4] for d in range(3, 10)))
        self.assertTrue(all(rows[d][6] < 729/256 for d in rows))
        # older value for four turns (small-volume, unit differing by 2^(1/3)) is not assumed here; three turns is
        self.assertTrue(rows[3][2] <= 4.1167*2**(1/3) <= rows[3][3] + 1e-4)
        with open(os.path.join(os.path.dirname(os.path.abspath(dr.__file__)), 'DR1_RESULT.json')) as f:
            self.assertTrue(json.load(f)['all_pass'])


if __name__ == '__main__':
    unittest.main()
