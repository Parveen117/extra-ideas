"""CR1: stdlib only."""
from fractions import Fraction as F
import unittest

import cr1_core_rates as m


class CoreRatesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()

    def test_all(self):
        self.assertGreaterEqual(len(self.out), 19)
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_generator_by_hand(self):
        # L(e1^2) = 2 e1 L e1 + 2 Gamma(e1, e1) = 6d e1 + 4 e1 ; L(e1 e2) = e2 3d + e1 2(d-1) e1 + 2 . 4 e2
        self.assertEqual(m.gen({(2, 0, 0): F(1)}, 3), {(1, 0, 0): 22})
        self.assertEqual(m.gen({(1, 1, 0): F(1)}, 3), {(0, 1, 0): 17, (2, 0, 0): 4})
        self.assertEqual(m.gen({(0, 0, 2): F(1)}, 4), {(0, 1, 1): 8})       # 2 e3 (d-2) e2 + 2 . 2 e2 e3 = 2d e2 e3
        self.assertEqual(m.degree((1, 2, 3)), 14)

    def test_means_by_hand(self):
        mean = m.free_means(3, F(2))
        self.assertEqual(mean({(1, 0, 0): F(1)}), F(9, 4))                  # nine entries of variance 1/4
        self.assertEqual(mean({(2, 0, 0): F(1)}), F(9, 4)*F(11, 4))         # <r^4> = <r^2>(N/2 + 1)/w
        self.assertEqual(mean({(0, 1, 0): F(1)}), 3*F(3, 2)/4)              # three pairs, (3/2)/w^2 each

    def test_matrix_element_by_a_second_formula(self):
        # phi = e1 exp(-w e1/2): |grad phi|^2 = e1 (2 - w e1)^2 exp(-w e1), so <phi, h phi> = mean(e1 (2 - w e1)^2/2 + e1^2 e2)
        for d in (3, 4):
            w = m.RATE[d]
            basis, s, h = m.ladder(d, w, 2)
            i = basis.index((1, 0, 0))
            mean = m.free_means(d, w)
            direct = mean({(1, 0, 0): F(2), (2, 0, 0): -2*w, (3, 0, 0): w*w/2, (2, 1, 0): F(1)})
            self.assertEqual(h[i][i], direct)
            self.assertEqual(s[i][i], mean({(2, 0, 0): F(1)}))

    def test_two_routes_of_the_count(self):
        basis, s, h = m.ladder(3, m.RATE[3], 6)
        for mu in (F(4), F(5187, 1000), F(26, 5), F(7), F(81, 10), F(12)):
            self.assertEqual(m.below(s, h, mu), m.below_by_minors(s, h, mu))
        self.assertEqual(m.below(s, h, F(4)), 0)
        self.assertEqual(m.below(s, h, F(81, 10)), 2)

    def test_numbers(self):
        three = self.num['d = 3, degree 10: lowest and second rate of the ladder, and their difference']
        self.assertAlmostEqual(three[0], 5.1867, places=3)
        self.assertAlmostEqual(three[1], 8.0467, places=3)
        # cross-check only (general knowledge, the usual normalization differs by 2^(1/3)): 4.1167 and 6.386
        self.assertAlmostEqual(three[0]/2**(1/3), 4.1167, places=3)
        self.assertAlmostEqual(three[1]/2**(1/3), 6.3866, places=3)
        self.assertAlmostEqual((18225/256)**(1/3), 4.1447, places=3)
        self.assertAlmostEqual((2025/8)**(1/3), 6.3257, places=3)


if __name__ == '__main__':
    unittest.main()
