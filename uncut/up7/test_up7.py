"""UP7: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import up7_turn_and_cut as y


class TurnAndCutTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['group']['order'], 8)
        self.assertEqual(res['group']['images_of_T_V_S_P_and_F_U_H_G']['thermal'], {'ends': 'SVTP', 'corners': 'UFGH'})
        self.assertEqual(res['open_walk']['eigenvalues_of_iota_kappa'], ['-2', '1/2'])

    def test_no_cut_on_the_carrier_commutes_with_the_turn(self):
        a, b = sp.symbols('a b', real=True)
        M = a*y.I2 + b*y.R                                   # everything that commutes with the turn
        sols = sp.solve(list(M*M - y.I2), [a, b], dict=True)
        self.assertEqual(sorted((s[a], s[b]) for s in sols), [(-1, 0), (1, 0)])

    def test_any_self_dagger_cut_gives_the_same_group(self):
        kap = y.cut(sp.Rational(3, 5), sp.Rational(4, 5), 0)
        self.assertEqual(kap*kap, y.I2)
        self.assertEqual(y.dagger(kap), kap)
        self.assertEqual(len(y.group_generated([y.R, kap])), 8)
        self.assertEqual((y.R*kap)**2, y.I2)

    def test_on_a_doubled_carrier_a_commuting_cut_needs_eight_steps(self):
        Z = sp.zeros(2, 2)
        turn = sp.Matrix(sp.BlockMatrix([[y.R, Z], [Z, y.R]]))
        swap = sp.Matrix(sp.BlockMatrix([[Z, y.I2], [y.I2, Z]]))
        self.assertEqual(turn*swap, swap*turn)
        chi = turn*swap
        self.assertEqual(chi**2, -sp.eye(4))
        self.assertEqual(chi**4, sp.eye(4))

    def test_miss_after_four_steps(self):
        w = sp.Rational(3, 4)
        kap = y.cut(sp.Rational(5, 4), 0, w)
        chi = y.R*kap
        T = sp.Matrix([sp.Rational(2), sp.Rational(7)])
        self.assertEqual(chi*chi*T - T, -2*w*chi*T)

    def test_dropping_the_turn_part_of_a_cut_is_rejected(self):
        with patch.object(y, 'cut', lambda u, v, w: u*y.K + v*y.S):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
