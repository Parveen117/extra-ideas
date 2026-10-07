"""BC1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import bc1_boundary_and_count as y


class BoundaryCountTests(unittest.TestCase):
    def test_all(self):
        res = y.run()
        self.assertTrue(res['boundary_never_returns'])
        self.assertEqual(res['clock'][0], ('3/5', '4/5', '12/5'))
        self.assertEqual(res['loop_areas'][0], '6')

    def test_interior_does_return(self):
        H = [[F(5), F(3)], [F(3), F(5)]]
        RH = y.mul(y.R, H)
        self.assertEqual(y.turn(F(1), F(0), F(4), RH), y.ID)
        quarter = y.turn(F(0), F(1), F(4), RH)
        self.assertEqual(y.mul(y.mul(quarter, quarter), y.mul(quarter, quarter)), y.ID)

    def test_wrong_generator_is_rejected(self):
        with patch.object(y, 'R', [[F(0), F(1)], [F(1), F(0)]]):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
