"""UP3: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import up3_sequence_of_diagrams as y

a, b = y.a, y.b


class SequenceOfDiagramsTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['power_form']['chi_1'], '1/4')
        self.assertEqual(res['not_a_power']['above'][3], '11636647/11469655')
        self.assertEqual(res['not_a_power']['below'][3], '-70384/167315')

    def test_the_rule_on_a_gas_with_constant_capacities(self):
        cv, R = sp.symbols('c_v R', positive=True)
        T, V = a, b
        S = cv*sp.log(T) + R*sp.log(V)
        P = R*T/V
        T2, S2, P2, V2 = [sp.simplify(q) for q in y.step(T, S, P, V)]
        self.assertEqual(sp.simplify(P2 - S*T/(cv + R)), 0)          # the owner's P_2 = S_1 T_1 / C_{P,1}
        self.assertEqual(sp.simplify(V2 - S*T/cv), 0)
        self.assertEqual(sp.simplify(T2 + P), 0)
        self.assertEqual(sp.simplify(S2 + (cv + R)/cv*P), 0)
        self.assertEqual(sp.simplify(y.br(T2, S2)), 0)               # the two ends of the axis are one reading

    def test_four_readings_with_no_potential(self):
        T, S, P, V = a, b + a**2, a*b + 1, a + b**3
        self.assertNotEqual(sp.simplify(y.br(T, S) - y.br(P, V)), 0)
        self.assertNotEqual(sp.simplify(y.br(T, S) + y.br(P, V)), 0)
        T2, S2, P2, V2 = y.step(T, S, P, V)
        at = {a: sp.Rational(2), b: sp.Rational(3)}
        self.assertEqual(sp.simplify((P2*S2 - V2*T2).subs(at)), 0)

    def test_exchanging_the_ends_of_one_axis_inverts_the_ratio(self):
        T, S, P, V = a, b + a**2, a*b + 1, a + b**3
        self.assertEqual(sp.simplify(y.cross_ratio(S, T, P, V)*y.cross_ratio(T, S, P, V) - 1), 0)
        self.assertEqual(sp.simplify(y.cross_ratio(T, S, V, P)*y.cross_ratio(T, S, P, V) - 1), 0)

    def test_a_rule_with_the_wrong_fixed_reading_is_rejected(self):
        def wrong(T, S, P, V):
            return (V*y.br(P, T)/y.br(V, T), V*y.br(P, S)/y.br(V, S), S*y.br(T, V)/y.br(S, V), S*y.br(T, P)/y.br(S, P))
        with patch.object(y, 'step', wrong):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
