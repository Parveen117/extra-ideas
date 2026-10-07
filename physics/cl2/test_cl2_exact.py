"""CL2: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import cl2_two_histories as y

G = F(3, 2)


class HistoryTests(unittest.TestCase):
    def test_count_of_a_leg(self):
        self.assertEqual(y.leg_control(G)[0]['count'], '6')

    def test_a_spacelike_leg_is_rejected(self):
        with patch.object(y, 'LEGS', [(F(5), F(0), F(3), F(4))]):
            with self.assertRaises(ValueError):
                y.leg_control(G)

    def test_corner_cost_and_collinear_equality(self):
        row = y.corner_control(G)
        self.assertEqual(row['pairs'], 36)
        self.assertEqual(row['strict'], 35)
        self.assertEqual(row['example']['straight_squared'], '100')
        self.assertEqual(row['example']['turned_count_squared'], '64')

    def test_split_character(self):
        rows = y.hyperbola_control(F(5, 4), F(3, 4), F(2), levels=3)
        self.assertEqual(rows[-1]['ratio'], '1')
        self.assertEqual(rows[-2]['legs'], 2)
        # two legs: ratio = 2 / (2 cosh b) = 1 / cosh b, b = two base steps: cosh = 17/8
        self.assertEqual(rows[-2]['ratio'], '8/17')

    def test_a_non_boost_is_rejected(self):
        with self.assertRaises(ValueError):
            y.hyperbola_control(F(5, 4), F(1), F(1))

    def test_both_sheets(self):
        row = y.sheet_control(G)
        self.assertEqual(F(row['mass']), -F(row['antimass']))


if __name__ == '__main__':
    unittest.main()
