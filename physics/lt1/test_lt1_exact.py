"""LT1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import lt1_lambda_tower_inward as y


class LambdaTowerTests(unittest.TestCase):
    def test_tower_moves_inward_to_the_horizon(self):
        res = y.run()
        self.assertEqual(res['inward'][0], ('2', '9/8'))
        self.assertEqual(res['two_references'][0], ('1', '25/16', '3/5'))
        self.assertEqual(res['invariant_per_generation'][:2], ['16/25', '64/289'])

    def test_horizon_and_infinity_are_fixed(self):
        self.assertEqual(y.radius_map(F(1), F(1)), F(1))
        self.assertEqual(y.square_m(F(0)), 0)
        self.assertEqual(y.square_m(F(1)), 1)

    def test_a_different_generation_law_is_rejected(self):
        with patch.object(y, 'square_m', lambda m: 2*m):
            with self.assertRaises(ValueError):
                y.run()

    def test_lambda_end_points(self):
        q = F(3, 2)
        self.assertEqual(y.ratio_after(q, F(0)), q**2)
        self.assertLess(y.ratio_after(q, F(10**6)), q**4)


if __name__ == '__main__':
    unittest.main()
