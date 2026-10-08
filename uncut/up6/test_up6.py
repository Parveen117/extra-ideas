"""UP6: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import up6_eight_scale_operations as y

a, b = y.a, y.b


class EightScaleOperationsTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['inside_a_corner'], 'commute')
        self.assertEqual(res['witnesses_defects_at_(1/2,1/2)']['not_a_power'], ['5/18', '-10/27', '1290275/1753182'])

    def test_gas_with_constant_capacities(self):
        # U = V^(-R/c) exp(S/c): P V = R T.  Forms (A) and (B) commute; the U and F corners do not.
        c, R = sp.Rational(3, 2), sp.Integer(1)
        U = b**(-R/c)*sp.exp(a/c)
        T, P, ops = y.setup(U)
        self.assertEqual(sp.simplify(P*b - R*T), 0)
        eT = sp.simplify(ops['F_V'](P)/P)
        self.assertEqual(eT, -1)
        self.assertEqual(sp.simplify(ops['F_V'](sp.log(-eT))), 0)
        m = sp.simplify(ops['F_V'](a)/a)
        self.assertEqual(sp.simplify(m - R/a), 0)
        self.assertEqual(sp.simplify(ops['U_S'](m) + R/a), 0)          # not zero: the cost is -R/S times S T = -R T

    def test_form_a_times_form_b(self):
        U = a**2/2 + a*b/3 + b**2/2 + a**2*b/5 + b**3/7
        T, P, ops = y.setup(U)
        A = -ops['G_P'](b)
        B = ops['F_V'](P)
        self.assertEqual(sp.simplify(A*B + P*b), 0)

    def test_one_corner_only_costs_nothing(self):
        U = sp.Function('U')(a, b)
        T, P, ops = y.setup(U)
        self.assertEqual(sp.simplify(ops['U_S'](ops['U_V'](U)) - ops['U_V'](ops['U_S'](U))), 0)

    def test_a_wrong_fixed_reading_is_rejected(self):
        real = y.setup
        def other(U):
            T, P, ops = real(U)
            ops['F_V'] = y.E(y.V, P)                                   # hold P instead of T
            return T, P, ops
        with patch.object(y, 'setup', other):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
