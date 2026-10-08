"""PM1: stdlib only."""
from fractions import Fraction as F
import math
import unittest

import pm1_prime_turn_series as m


class PrimeTurnSeriesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()
        cls.sh = m.shells(200)
        cls.pi = m.pi_interval()

    def test_all(self):
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_first_values_by_hand(self):
        # norm 5: eight numbers, (1 + 2 iota)^4 = -7 - 24 iota  ->  (1/4)(4(-7) + 4(-7)) = -14
        self.assertEqual(m.gpow((1, 2), 4), (-7, -24))
        self.assertEqual(self.num['A_1_first_values'][:5], [1, -4, 0, 16, -14])
        self.assertEqual(m.from_primes(65, 1), (-14)*(-238))

    def test_exponential_enclosure_against_floats(self):
        for q in (F(0), F(1, 3), F(7, 3), F(40)):
            lo, hi = m.exp_neg(q)
            self.assertLessEqual(lo, hi)
            self.assertLess(hi - lo, F(1, 10**100))
            self.assertAlmostEqual(float(lo), math.exp(-float(q)), delta=1e-15)

    def test_pi_enclosure_against_floats(self):
        self.assertLess(self.pi[0], self.pi[1])
        self.assertAlmostEqual(float(self.pi[0]), math.pi, delta=1e-15)

    def test_seam_law_in_plain_floats(self):
        def theta(t, k):
            return sum((complex(a, b)**(4*k)).real*math.exp(-math.pi*n*t)
                       for n in range(1, 200) for a, b in self.sh[n])
        self.assertAlmostEqual(theta(0.5, 1), 2.0**5*theta(2.0, 1), delta=1e-12)
        self.assertAlmostEqual(theta(0.5, 2)/(2.0**9*theta(2.0, 2)), 1.0, delta=1e-10)

    def test_a_changed_coefficient_breaks_the_seam_law(self):
        coefs = [0] + [m.shell_sum(self.sh[n], 4)[0] for n in range(1, 110)]
        coefs[1] += 1
        left = m.weighted_sum(coefs, m.exp_neg_pi(F(1, 2), self.pi), F(1), 3)
        right = m.weighted_sum(coefs, m.exp_neg_pi(F(2), self.pi), F(32), 3)
        self.assertFalse(m.overlap(left, right))

    def test_negative_weight_gives_a_proper_enclosure(self):
        coefs = [0] + [m.shell_sum(self.sh[n], 4)[0] for n in range(1, 110)]
        e = m.exp_neg_pi(F(2), self.pi)
        pos = m.weighted_sum(coefs, e, F(32), 3)
        neg = m.weighted_sum(coefs, e, F(-32), 3)
        self.assertLessEqual(neg[0], neg[1])
        self.assertTrue(m.overlap(neg, (-pos[1], -pos[0])))          # the same number: must overlap
        self.assertFalse(m.overlap(neg, pos))                        # Theta_1(2) is not zero
        with self.assertRaises(AssertionError):
            m.overlap((F(1), F(0)), pos)                             # a reversed pair is refused

    def test_character_mod_four(self):
        self.assertEqual([m.chi4(n) for n in range(1, 9)], [1, 0, -1, 0, 1, 0, -1, 0])


if __name__ == '__main__':
    unittest.main()
