"""QC5: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import qc5_moving_share as y

X = (F(1, 2), F(-1, 3), F(2, 5), F(1))


class MovingShareTests(unittest.TestCase):
    def test_laws_at_a_point(self):
        rows = y.point_control(X, F(4, 9))
        self.assertTrue(rows['logistic']['self_dual'])
        self.assertFalse(rows['steeper']['self_dual'])
        self.assertEqual(rows['equal']['action_density_r4'], '3/4')
        self.assertEqual(rows['winding only']['action_density_r4'], '0')

    def test_scale_is_free(self):
        for l2 in (F(1, 100), F(7)):
            self.assertTrue(y.point_control(X, l2)['logistic']['self_dual'])

    def test_dropping_the_variance_term_is_detected(self):
        real = y.field

        def broken(x, q, dq):
            return real(x, q, dq) if q == 0 else {k: y.lin((1, v), (F(1, 7), y.UNIT[1])) for k, v in real(x, q, dq).items()}
        with patch.object(y, 'field', broken):
            with self.assertRaises(ValueError):
                y.point_control(X, F(4, 9))

    def test_wrong_share_law_is_not_self_dual(self):
        real = dict(y.SHARES)
        with patch.dict(y.SHARES, {'logistic': real['steeper']}):
            with self.assertRaises(ValueError):
                y.point_control(X, F(4, 9))

    def test_memory_form_of_the_cut(self):
        row = y.lift_control(X, F(3))
        self.assertTrue(row['action_density'])
        with patch.object(y, 'conj', lambda a: a):
            with self.assertRaises(ValueError):
                y.lift_control(X, F(3))

    def test_reduced_action(self):
        self.assertEqual(y.reduced_action()['minimum_reduced_action'], '2')


if __name__ == '__main__':
    unittest.main()
