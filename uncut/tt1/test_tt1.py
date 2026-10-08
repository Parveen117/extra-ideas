"""TT1: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import tt1_time_in_the_thermo_diagram as y


class TimeInTheThermoDiagramTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['numbers_w_3_4']['clock_factor'], '8/17 = 1/(1 + 2 w^2)')
        self.assertEqual(res['gas']['closed_form'], 'tanh^2(eta) = (R c)^2 / (S^2 (S + c)^2 + (R c)^2)')

    def test_self_dagger_cut_cycle_is_the_identity(self):
        kap = sp.Rational(3, 5)*y.K + sp.Rational(4, 5)*y.S
        self.assertEqual((y.R*kap)**2, y.I2)

    def test_cycles_compose_by_adding_rapidity(self):
        # tanh(eta) = 3/5 ; one cycle 15/17 ; two cycles: add 15/17 to 15/17
        kap = sp.Rational(5, 4)*y.K + sp.Rational(3, 4)*y.R
        one = (y.R*kap)**2
        two = one*one
        a1, u1, v1, w1 = y.parts(one)
        a2, u2, v2, w2 = y.parts(two)
        x1 = sp.sqrt(u1*u1 + v1*v1)/a1
        x2 = sp.sqrt(u2*u2 + v2*v2)/a2
        self.assertEqual(x1, sp.Rational(15, 17))
        self.assertEqual(sp.simplify(x2 - 2*x1/(1 + x1*x1)), 0)
        self.assertEqual(w1, 0)
        self.assertEqual(w2, 0)

    def test_a_response_inside_the_cone_has_no_cut(self):
        L = sp.Matrix([[1, -2], [2, 1]])
        X = L - L.trace()/2*y.I2
        self.assertEqual(X*X, -4*y.I2)                                 # a turn, not a cut
        self.assertEqual(sorted(sp.im(e) for e in L.eigenvals()), [-2, 2])

    def test_reciprocal_response_has_no_time_part(self):
        L = sp.Matrix([[3, 1], [1, 2]])
        self.assertEqual(y.parts(L)[3], 0)

    def test_a_wrong_cycle_law_is_rejected(self):
        with patch.object(y, 'Exp_cut', lambda eta, n: sp.cosh(eta)*y.I2 - sp.sinh(eta)*n):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
