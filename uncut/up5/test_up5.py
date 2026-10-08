"""UP5: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import up5_seed_hierarchy_checked as y

a, b = y.a, y.b


class SeedHierarchyTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertIn('one equation', res['four_relations'])
        self.assertIn('F U', res['defect'])

    def test_seed_layer_on_a_potential(self):
        # U(S, V) = S^3 / V on the chart (a, b) = (S, V): the seed's layer 1, exact
        S, V = a, b
        U = S**3/V
        T, P = sp.diff(U, S), -sp.diff(U, V)
        L = {k: sp.simplify(v) for k, v in y.seed_layer(T, V, S, P).items()}
        self.assertEqual(sp.simplify(L['lp']/L['lv']), sp.Rational(1, 4))
        self.assertEqual(sp.simplify(L['zt']/L['ls']), sp.Rational(1, 4))
        self.assertEqual(sp.simplify(L['lt']*L['zt'] + P*V), 0)
        self.assertEqual(sp.simplify(y.br(T, S) - y.br(P, V)), 0)          # layer 0 has its one equation

    def test_layer_one_does_not_inherit_the_equation(self):
        S, V = a, b
        U = S**2/2 + S*V/3 + V**2/2 + S**2*V/5 + V**3/7
        T, P = sp.diff(U, S), -sp.diff(U, V)
        L = y.seed_layer(T, V, S, P)
        t1, v1, s1, p1 = L['zt'], L['lv'], L['ls'], L['lp']               # the closed layer (z_t in the t place)
        gap = sp.simplify((y.br(t1, s1) - y.br(p1, v1)).subs({a: sp.Rational(1, 2), b: sp.Rational(1, 2)}))
        gap2 = sp.simplify((y.br(t1, s1) + y.br(p1, v1)).subs({a: sp.Rational(1, 2), b: sp.Rational(1, 2)}))
        self.assertNotEqual(gap, 0)
        self.assertNotEqual(gap2, 0)

    def test_a_changed_definition_is_rejected(self):
        def other(T, V, S, P):
            L = y.seed_layer.__wrapped__(T, V, S, P) if hasattr(y.seed_layer, '__wrapped__') else None
            lp = -S*T/(T*y.at_fixed(S, T, V))                               # capacities exchanged
            lv = -S*T/(T*y.at_fixed(S, T, P))
            ls = V*y.at_fixed(P, V, S)
            lt = -P*y.at_fixed(V, P, T)
            return dict(lp=lp, lv=lv, ls=ls, lt=lt, zs=-P*V/ls, zt=-P*V/lt)
        with patch.object(y, 'seed_layer', other):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
