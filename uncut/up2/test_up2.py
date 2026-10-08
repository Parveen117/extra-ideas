"""UP2: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import up2_two_point_potential as y


class TwoPointPotentialTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['one_pair']['depths'], '1/K + 1/K* = 1')
        self.assertEqual(res['two_pairs'], 'rest set of dimension 2 in 4; doubled response of rank 2')

    def test_duality_at_rational_points(self):
        for K, a, b in ((3, 4, 9), (2, 5, 7), (sp.Rational(3, 2), 4, 9), (5, 1, 32)):
            self.assertTrue(y.dual_check(K, a, b))

    def test_sign_follows_the_response(self):
        p, q = sp.symbols('p q', real=True)
        stable = (2*p*p + 2*p*q + 3*q*q)/2                       # positive response
        saddle = (2*p*p + 6*p*q + 3*q*q)/2                       # outside the disc
        Gs = y.two_point(stable, (p, q), (1, 1))
        Gd = y.two_point(saddle, (p, q), (1, 1))
        pts = [(0, 0), (3, -2), (-1, 4), (2, 1), (1, 0)]
        self.assertTrue(all(Gs.subs({p: a, q: b}) > 0 for a, b in pts))
        vals = [Gd.subs({p: a, q: b}) for a, b in pts]
        self.assertTrue(min(vals) < 0 < max(vals))

    def test_classical_degree_one_has_no_second_half(self):
        # degree 1: the first layer does not scale, so the y-half has weight 0 and no conjugate degree
        K = sp.Symbol('K', positive=True)
        self.assertEqual(sp.limit(K/(K - 1), K, 1, '+'), sp.oo)
        a, b = sp.symbols('a b', positive=True)
        one = a**2/b
        H = sp.hessian(one, (a, b))
        self.assertEqual(sp.simplify(H.det()), 0)                # its own response is already null

    def test_degree_two_is_symmetric(self):
        x, yv = sp.symbols('x y', positive=True)
        G = x**2/2 + yv**2/2 - x*yv
        self.assertEqual(sp.simplify(G - (x - yv)**2/2), 0)
        self.assertEqual(sp.simplify(G.subs({x: yv, yv: x}, simultaneous=True) - G), 0)

    def test_a_potential_that_keeps_the_first_layer_is_rejected(self):
        def wrong(phi, xs, base):
            return phi - phi.subs(dict(zip(xs, base)), simultaneous=True)
        with patch.object(y, 'two_point', wrong):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
