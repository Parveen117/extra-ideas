"""UP4: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import up4_cost_of_commuting_cuts as y

a, b = y.a, y.b


class CostOfCommutingCutsTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['witness']['chi'], '19/16')
        self.assertEqual(res['witness']['symmetric_part'], 'positive')

    def test_commuting_cuts_cost_nothing(self):
        U = sp.Function('U')(a, b)
        one, nil = sp.Integer(1), sp.Integer(0)
        self.assertEqual(y.defect(one, nil), (0, 0))
        T, P, A, B1, B2, C = y.response(U, one, nil)
        self.assertEqual(sp.simplify(B2 - B1), 0)

    def test_a_defect_along_the_other_operation(self):
        U = sp.Function('U')(a, b)
        ff, gg = sp.Integer(1), a                      # D_V = d_b + a d_a : [D_S, D_V] = D_S
        al, be = y.defect(ff, gg)
        self.assertEqual((sp.simplify(al), sp.simplify(be)), (1, 0))
        T, P, A, B1, B2, C = y.response(U, ff, gg)
        self.assertEqual(sp.simplify(B2 - B1 - T), 0)

    def test_small_defect_keeps_the_ratio_below_one(self):
        ff, gg = sp.exp(a/4), sp.Integer(0)            # beta = 1/4
        U = a**2 + b**2 - a*b/2 + b/2
        T, P, A, B1, B2, C = [sp.simplify(q.subs({a: 0, b: 0})) for q in y.response(U, ff, gg)]
        chi = (A*C - B1*B2)/(A*C)
        self.assertTrue(0 < chi < 1)
        self.assertNotEqual(B1, B2)

    def test_ignoring_the_defect_is_rejected(self):
        with patch.object(y, 'defect', lambda *args: (sp.Integer(0), sp.Integer(0))):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
