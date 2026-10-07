"""LN1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import ln1_lost_and_returned as y


class LostAndReturnedTests(unittest.TestCase):
    def test_all(self):
        res = y.run()
        self.assertEqual(res['two_cuts'][0], ('5', '3', '4', '11'))
        self.assertEqual(res['three_cuts']['squares'], ['4', '1', '9'])
        self.assertEqual(res['gravity'][0], ('2', '9/16', '9/16', '9/8'))

    def test_lost_is_the_coupling_squared(self):
        H = y.sym([F(7), F(1), F(2)], {(0, 1): F(3), (1, 2): F(-4)})
        self.assertEqual(y.lost(H), [F(9), F(25), F(16)])

    def test_a_reading_with_no_coupling_loses_nothing(self):
        self.assertEqual(y.lost(y.sym([F(3), F(8)], {})), [F(0), F(0)])

    def test_a_layer_that_is_not_the_square_fails(self):
        with patch.object(y, 'mul', lambda a, b: a):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
