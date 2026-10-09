"""SL1: stdlib only."""
from fractions import Fraction as F
from math import exp, pi, sqrt
import unittest

import sl1_seam_law_from_the_walk as m


def sums(t, c=pi, top=3000):
    a = 1 + 2*sum(exp(-c*t*n*n) for n in range(1, top))
    half = 2*sum(exp(-c*t*(n + .5)**2) for n in range(top))
    b = 1 + 2*sum((-1)**n*exp(-c*t*n*n) for n in range(1, top))
    return a*a, b*b, half*half


def mean(x, y):
    for _ in range(60):
        x, y = (x + y)/2, sqrt(x*y)
    return x


class SeamLawTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()

    def test_all(self):
        self.assertGreaterEqual(len(self.out), 16)
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_series_by_hand(self):
        a = m.series(lambda n: 4*n*n)
        self.assertEqual({k: v for k, v in a.items() if k <= 40}, {0: 1, 4: 2, 16: 2, 36: 2})     # 1 + 2q + 2q^4 + 2q^9
        s = m.smul(a, a)
        self.assertEqual([s.get(4*k, 0) for k in range(6)], [1, 4, 4, 0, 4, 8])                  # counts of norms 0..5
        b = m.series(lambda n: 4*n*n, lambda n: -1 if n % 2 else 1)
        self.assertEqual(b[4], -2)

    def test_chain_against_floats(self):
        for u in (0.5, 1.0, 2.0):
            s, f, l = sums(u)
            big = mean(s, l)
            self.assertAlmostEqual(u*big, 1, places=12)
            self.assertAlmostEqual(mean(s, f), 1, places=12)
            for n in (1, 2, 3):
                sn, _, ln = sums(u/2**n)
                self.assertLessEqual(ln/2**n, big*(1 + 1e-12))
                self.assertLessEqual(big, sn/2**n*(1 + 1e-12))
            s2, _, l2 = sums(2*u)
            self.assertAlmostEqual(s, s2 + l2, places=12)                                         # S = S' + L'
            self.assertAlmostEqual(l*l, 4*s2*l2, places=12)

    def test_weak_end_count_against_floats(self):
        for t in (0.3, 0.05, 0.01):
            s, _, l = sums(t)
            h = sqrt(t)
            self.assertTrue(1 - h <= h*sqrt(s) <= 1 + h)
            self.assertTrue(1 - h <= h*sqrt(l) <= 1 + 2*h)

    def test_another_constant(self):
        # with Exp(-c u) the same walk gives u M(S, L) = pi/c
        for c in (2.0, 3.0, 5.0):
            s, f, l = sums(1.0, c)
            self.assertAlmostEqual(mean(s, l), pi/c, places=12)
            self.assertAlmostEqual(mean(s, f), 1, places=12)

    def test_seam_law_against_floats(self):
        for u in (0.4, 1.7, 3.0):
            s, f, l = sums(u)
            sv, fv, lv = sums(1/u)
            self.assertAlmostEqual(sv/(u*s), 1, places=11)
            self.assertAlmostEqual(fv/(u*l), 1, places=10)
            self.assertAlmostEqual(lv/(u*f), 1, places=10)


if __name__ == '__main__':
    unittest.main()
