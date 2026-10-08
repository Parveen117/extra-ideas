"""AG1: stdlib only."""
from fractions import Fraction as F
import math
import unittest

import ag1_diagonal_walk as m


class DiagonalWalkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()
        cls.pi = m.pm1.pi_interval()

    def test_all(self):
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_first_counts(self):
        s, f, g = m.series(10)
        self.assertEqual(s, [1, 4, 4, 0, 4, 8, 0, 0, 4, 4, 8])
        self.assertEqual(f, [1, -4, 4, 0, 4, -8, 0, 0, 4, -4, 8])

    def test_values_against_floats(self):
        # cross-check only: plain floating sums
        def sums(t):
            q = math.exp(-math.pi*t)
            r = range(-40, 41)
            s = sum(q**(a*a + b*b) for a in r for b in r)
            f = sum((-1)**(a + b)*q**(a*a + b*b) for a in r for b in r)
            l = sum(q**((a + .5)**2 + (b + .5)**2) for a in r for b in r)
            return s, f, l
        for t in (0.5, 1.0, 2.0):
            s, f, l = sums(t)
            se, fe, le = m.sfl(F(t).limit_denominator(8), self.pi)
            self.assertAlmostEqual(float(se[0]), s, delta=1e-12)
            self.assertAlmostEqual(float(fe[0]), f, delta=1e-12)
            self.assertAlmostEqual(float(le[0]), l, delta=1e-12)
            self.assertAlmostEqual(s*s, f*f + l*l, delta=1e-10)
        s, f, l = sums(1.0)
        self.assertAlmostEqual(f/s, 1/math.sqrt(2), delta=1e-12)

    def test_turn_block_against_floats(self):
        for u in (1.0, 0.25):
            kp = sum(mm*mm*math.exp(-math.pi*u*mm*mm) for mm in range(1, 60))
            km = sum((-1)**(mm - 1)*mm*mm*math.exp(-math.pi*u*mm*mm) for mm in range(1, 60))
            a, b = m.k_contents(F(u).limit_denominator(8), self.pi)
            self.assertAlmostEqual(float(a[0]), kp, delta=1e-12)
            self.assertAlmostEqual(float(b[0]), km, delta=1e-12)

    def test_a_wrong_walk_is_separated(self):
        # the flip of the doubled cell is not the arithmetic mean of anything here: F(2t) against (S + F)/2
        s_a, f_a, _ = m.sfl(F(1, 2), self.pi)
        s_b, f_b, _ = m.sfl(F(1), self.pi)
        self.assertFalse(m.overlap(f_b, m.iscale(m.iadd(s_a, f_a), F(1, 2))))
        self.assertTrue(m.overlap(s_b, m.iscale(m.iadd(s_a, f_a), F(1, 2))))

    def test_interval_helpers(self):
        self.assertEqual(m.iscale((F(1), F(2)), -3), (F(-6), F(-3)))
        lo, hi = m.isqrt_iv((F(2), F(2)))
        self.assertTrue(lo*lo <= 2 <= hi*hi and hi - lo < F(1, 10**100))
        with self.assertRaises(AssertionError):
            m.overlap((F(1), F(0)), (F(0), F(1)))


if __name__ == '__main__':
    unittest.main()
