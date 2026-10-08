"""RW1: stdlib only."""
from fractions import Fraction as F
import unittest

import rw1_returned_winding_count as m


class ReturnedWindingTests(unittest.TestCase):
    def test_all(self):
        out, _ = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)

    def test_first_coefficients_of_the_mean(self):
        self.assertEqual(m.mean_series(0, 3), [0, 1, F(-1, 2), F(1, 3)])
        self.assertEqual(m.mean_series(2, 2), [0, F(1, 3), F(-1, 36)])

    def test_enclosure_is_two_sided_and_tight(self):
        lo, hi = m.mean_interval(1, F(4))
        self.assertLess(lo, hi)
        self.assertLess(hi - lo, F(1, 10**30))

    def test_a_wrong_spread_constant_fails(self):
        low, up = m.spread_series(1, 40, lower=F(3, 5), upper=F(9, 10))
        self.assertTrue(any(x < 0 for x in low))
        self.assertTrue(any(x < 0 for x in up))

    def test_centre_coupling_against_the_log_form(self):
        # |k_c - (1/2) ln(4 kappa/3)| < 1/(2 (kappa - 2)) for kappa >= 3   (floats, as a cross-check only)
        import math
        for kappa in (F(3), F(10), F(100), F(1000)):
            lo, hi = m.mean_interval(1, kappa*kappa/4)
            r = float(2*lo/kappa)
            kc = math.atanh(r)
            self.assertLess(abs(kc - 0.5*math.log(4*float(kappa)/3)), 1/(2*(float(kappa) - 2)))
            self.assertTrue(float(kappa)/2 < math.sinh(2*kc) < 2*float(kappa)/3)

    def test_the_inequality_can_fail_for_a_wrong_sector(self):
        # reading the nu = 1 sector with the nu = 3 bounds must fail
        lo, hi = m.mean_interval(1, F(4))
        self.assertFalse(hi*hi + F(7, 2)*hi < 4 < lo*lo + 4*lo)


if __name__ == '__main__':
    unittest.main()
