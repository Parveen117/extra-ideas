"""UP1: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import up1_degree_of_the_potential as y


class DegreeOfThePotentialTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['tower'], 'lambda and place enter only through lambda s^(k-2)')
        self.assertEqual(res['ends']['lambda m -> infinity'], 'rho -> 2 rho/(1 + rho^2), no lambda left')

    def test_degree_one_and_degree_two(self):
        a, b = sp.symbols('a b', positive=True)
        one = a**2/b                                   # degree 1: the far corner is silent
        self.assertEqual(sp.simplify(y.partial_cut(one, (a, b), [1, 1])), 0)
        self.assertEqual(sp.simplify(y.partial_cut(one, (a, b), [sp.Rational(1, 2)]*2) - one/2), 0)
        two = 3*a**2 - a*b + 5*b**2                    # degree 2: the centre is silent
        self.assertEqual(sp.simplify(y.partial_cut(two, (a, b), [sp.Rational(1, 2)]*2)), 0)
        self.assertEqual(sp.simplify(y.partial_cut(two, (a, b), [1, 1]) + two), 0)

    def test_cf_benchmark(self):
        # CF-4 (11): nu = 3/2, kappa = 2/3, r = t^4, W = (2/3) t^6
        t = sp.Symbol('t', positive=True)
        W = sp.Rational(2, 3)*(t**4)**sp.Rational(3, 2)
        nu = sp.Rational(3, 2)
        self.assertEqual(sp.simplify((1 - nu/2)*W - t**6/6), 0)
        self.assertEqual(sp.simplify((1 - nu)*W + t**6/3), 0)

    def test_tower_covariance_at_a_rational_point(self):
        a, b = sp.symbols('a b', positive=True)
        phi = a**4/b                                   # degree 3
        H = sp.hessian(phi, (a, b))
        at = lambda M, p, q: M.subs({a: p, b: q})
        lam = sp.Rational(1, 5)
        left = y.tower(at(H, 2, 6), lam)
        right = 2**(3 - 2)*y.tower(at(H, 1, 3), lam*2**(3 - 2))
        self.assertEqual(left, right)

    def test_a_mixed_potential_has_no_single_silent_depth(self):
        a, b, t = sp.symbols('a b t', positive=True)
        mixed = a**2 + a*b + a
        sol = sp.solve(sp.Eq(y.partial_cut(mixed, (a, b), [t, t]), 0), t)
        self.assertTrue(all(q.has(a) or q.has(b) for q in sol))

    def test_a_tower_of_another_shape_is_rejected(self):
        with patch.object(y, 'tower', lambda H, l: H + l*H):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
