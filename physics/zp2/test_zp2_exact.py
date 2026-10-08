"""ZP2: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import zp2_floor_outward_certificate as y


class FloorCertificateTests(unittest.TestCase):
    def test_run_and_pinned_values(self):
        res = y.run()
        self.assertEqual(res['line'][2], {'2': '1/12', '0': '-1/12'})
        self.assertEqual(res['line'][1], {'2': '1/6', '0': '-1/6'})
        self.assertEqual(res['plates'][3], {'1': '1/140', '-3': '-1/360', '-5': '1/252'})
        self.assertEqual(res['plates'][2], {'1': '1/60'})
        self.assertEqual(res['line_tail_e_m']['2'], '0')

    def test_exact_at_every_integer_cutoff_for_m2(self):
        for yy in (1, 2, 3, 17, 100):
            self.assertEqual(y.line_sum(2, yy), F(yy*yy - 1, 12))

    def test_sharp_and_kinked_cutoffs_give_other_constants(self):
        self.assertEqual(y.line_sum(0, 9), F(45))                # y^2/2 + y/2, constant 0
        self.assertEqual(y.line_sum(1, 9), F(81 - 1, 6))         # constant -1/6
        self.assertEqual(y.plate_sum(2, 11), F(11, 60))          # no d^-3 term at all

    def test_force_on_a_wall_inside_a_fixed_length(self):
        # units pi*kappa*c/2 = 1, lattice of cut-off cells: energy(L) = line_sum(2, y)/L with y = W L
        W, total = 12, 10
        def energy(a):
            return y.line_sum(2, W*a)/a + y.line_sum(2, W*(total - a))/(total - a)
        for a in range(1, total):
            self.assertEqual(energy(a), F(W*W*total, 12) - F(1, 12)*(F(1, a) + F(1, total - a)))
        self.assertLess(energy(1), energy(5))                    # pulled toward the nearer end

    def test_interpolation_refuses_a_wrong_degree(self):
        with self.assertRaises(ValueError):
            y.laurent(y.line_sum, 5, 5, 4)

    def test_a_changed_floor_weight_is_rejected(self):
        with patch.object(y, 'line_sum', lambda m, yy: sum(n*n*(1 - F(n, yy))**m for n in range(yy + 1))):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
