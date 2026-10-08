"""PT1: stdlib only."""
from fractions import Fraction as F
import unittest

import pt1_prime_turns as y


class PrimeTurnTests(unittest.TestCase):
    def test_all(self):
        res = y.run()
        self.assertEqual(res['prime_turns'][5], ('3/5', '4/5'))
        self.assertEqual(res['factored_turns'], 48)
        self.assertEqual(res['rational_cosine_orders'], {-2: 2, -1: 3, 0: 4, 1: 6, 2: 1})

    def test_unique_factorisation_of_a_composite_turn(self):
        u = y.tmul(y.b2((2, 1)), y.b2((3, 2)))            # 5 and 13
        self.assertEqual(y.factor_turn(u), (0, {5: 1, 13: 1}))
        v = y.tmul(y.b2((2, 1)), y.tpow(y.b2((3, 2)), -1))
        self.assertEqual(y.factor_turn(v), (0, {5: 1, 13: -1}))
        self.assertNotEqual(u, v)

    def test_not_a_turn_is_rejected(self):
        with self.assertRaises(ValueError):
            y.factor_turn((F(1, 2), F(1, 2)))

    def test_prime_three_mod_four_has_no_turn(self):
        self.assertIsNone(y.two_squares(7))
        self.assertIsNone(y.two_squares(11))
        self.assertEqual(y.two_squares(13), (3, 2))


if __name__ == '__main__':
    unittest.main()
