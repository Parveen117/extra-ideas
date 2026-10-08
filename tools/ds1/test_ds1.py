"""DS1: stdlib only."""
from fractions import Fraction as F
import unittest

import ds1_one_diagonal as m

G = m.G


class OneDiagonalTests(unittest.TestCase):
    def test_all(self):
        out, _ = m.run()
        for k, v in out.items():
            self.assertTrue(v, k)

    def test_rank_one_by_hand(self):
        # R = 1, D = v v*, v = (3/5, 4/5) t :  lost/seen = t^2 ;  count 0 below 1, kernel at 1, count 1 above
        for t, expect in ((F(1, 2), (2, 0, 0)), (F(1), (1, 0, 1)), (F(2), (1, 1, 0))):
            v = [[G(F(3, 5)*t)], [G(F(4, 5)*t)]]
            f = m.madd(m.eye(2), m.mmul(v, m.dagger(v)), -1)
            self.assertEqual(m.inertia(f), expect)
            b = m.mmul(m.dagger(v), v)
            self.assertEqual(b, [[G(t*t)]])

    def test_determinant_of_the_small_matrix(self):
        b = [[G(2), G(0, 1)], [G(0, -1), G(F(1, 2))]]
        self.assertEqual(m.charpoly_one_minus_zb(b), [F(1), F(-5, 2), F(0)])      # det B = 1 - 1 = 0, tr B = 5/2
        self.assertEqual(m.roots_with_multiplicity([F(1), F(-5, 2)], F(0), F(1)), 1)
        self.assertEqual(m.inertia(m.madd(b, m.eye(2), -1))[0], 1)

    def test_inertia_refuses_a_form_that_is_not_self_dagger(self):
        with self.assertRaises(AssertionError):
            m.inertia([[G(1), G(2)], [G(3), G(1)]])

    def test_diagonal_reading(self):
        z = G(1, 1)
        self.assertEqual(z/z.dag(), m.IOTA)
        self.assertEqual((G(2, 1)/G(2, 1).dag()), G(F(3, 5), F(4, 5)))       # PT1's u_5

    def test_odd_sector_margin_is_a_split_of_one(self):
        # RKF theorum/34-39 (cited, not re-run here): lost/seen <= 0.8290856201657449, margin 0.1709143798342551
        self.assertEqual(F(8290856201657449, 10**16) + F(1709143798342551, 10**16), 1)
        self.assertLess(F(8290856201657449, 10**16), 1)


if __name__ == '__main__':
    unittest.main()
