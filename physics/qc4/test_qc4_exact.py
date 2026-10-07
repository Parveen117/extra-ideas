"""QC4: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import qc4_return_is_memory as y


class ReturnIsMemoryTests(unittest.TestCase):
    def test_curvature_is_minus_variance_and_compass_law(self):
        rows = y.readings_control()
        self.assertEqual(len(rows), 3)
        for row in rows:
            self.assertEqual(row['shares']['U:G = 3:1'], '3/4')
            self.assertEqual(row['shares']['F:Hc = 1:1'], '0')
            self.assertEqual(row['shares']['four corners equal'], '1/2')

    def test_weights_that_do_not_sum_to_one_are_rejected(self):
        read, _ = y.readings(F(1, 2), F(1, 3))
        with self.assertRaises(ValueError):
            y.mixture_curvature(read, {'U': F(1, 2), 'G': F(1, 3)})

    def test_a_curved_reading_is_rejected(self):
        real = y.ph3.field

        def broken(x, y_):
            H, Hx, Hy, Hxy, Hyx = real(x, y_)
            return H, Hx, Hy, Hxy, [[F(2), F(5)], [F(5), F(0)]]
        with patch.object(y.ph3, 'field', broken):
            with self.assertRaises(ValueError):
                y.readings(F(1, 2), F(1, 3))

    def test_memory_is_variance(self):
        self.assertEqual(y.memory_control()[0]['memory'], '576/625')

    def test_loop_arrow(self):
        row = y.loop_control(F(1, 4))
        self.assertLess(F(row['memory_over_twice_defect']), 1)
        self.assertGreater(F(row['memory_over_twice_defect']), F(99, 100))
        self.assertEqual(y.matrix_loop_control()['mean_defect'], '380/501')

    def test_a_non_isometric_link_is_rejected(self):
        real = y.qc3.turn

        def broken(axis, r):
            q = list(real(axis, r))
            q[0] = q[0]+F(1, 50)
            return tuple(q)
        with patch.object(y.qc3, 'turn', broken):
            with self.assertRaises(ValueError):
                y.loop_control(F(1, 4))


if __name__ == '__main__':
    unittest.main()
