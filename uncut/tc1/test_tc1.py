"""TC1: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import tc1_three_cuts_turn_part as y


class ThreeCutTurnPartTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['determinant'], 'det rho = (h.h - k.k) + 2 iota h.k')
        self.assertIn('h = (4; 0,0,0), k = (0; 0,0,-3)', res['witness'])

    def test_light_field_is_the_traceless_case(self):
        # EM1: F = E + iota B, F.F = (E.E - B.B) + 2 iota E.B
        E, B = (1, 2, 2), (2, -1, 3)
        F = y.vec(0, E) + sp.I*y.vec(0, B)
        EE, BB, EB = sum(x*x for x in E), sum(x*x for x in B), sum(x*z for x, z in zip(E, B))
        self.assertEqual(sp.simplify(F.det() + (EE - BB) + 2*sp.I*EB), 0)
        self.assertEqual(sp.simplify(F*F - ((EE - BB) + 2*sp.I*EB)*y.I2), sp.zeros(2, 2))

    def test_opposite_directions_cancel(self):
        a = y.vec(4, (0, 0, 0)) + sp.I*y.vec(0, (1, 2, 2))
        b = y.vec(4, (0, 0, 0)) + sp.I*y.vec(0, (-1, -2, -2))
        self.assertEqual(a.det(), 25)
        self.assertEqual((a + b).det(), 64)
        self.assertEqual(y.split(a + b)[1], (0, [0, 0, 0]))

    def test_a_self_dagger_reading_has_no_second_half(self):
        rho = y.vec(5, (3, 0, 4))
        self.assertEqual(y.split(rho)[1], (0, [0, 0, 0]))
        self.assertEqual(rho.det(), 0)                               # a single reading: null (IN1-T2)

    def test_dagger_without_conjugation_is_rejected(self):
        with patch.object(y, 'dag', lambda M: M.T):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
