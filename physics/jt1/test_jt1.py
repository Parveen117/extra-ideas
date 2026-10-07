"""JT1 tests (sympy; exact rationals for the observer)."""
import unittest

import sympy as sp

import jt1_tower_of_the_frame_wave as y


class TowerTests(unittest.TestCase):
    def test_tower_observer(self):
        res = y.tower_observer()
        self.assertEqual(res['recognised_at_layer'], 3)

    def test_wave_tower(self):
        res = y.hyperbolic_form()
        self.assertEqual(res['layers'], {2: '1', 4: '1/3', 6: '2/45'})
        self.assertEqual(res['one_way'], '0')
        self.assertNotIn('phi', res['rigid_law'])

    def test_law_is_the_square_of_the_response_rate(self):
        M = sp.Matrix([[1 + y.T*y.Z, y.T], [y.T, 2 + y.Z**2]])
        self.assertEqual(sp.simplify(y.law_from_frame(M) - y.law_from_response(M)), 0)

    def test_volume_term_is_needed(self):
        M = sp.Matrix([[1 + y.T**2, 0], [0, 1 + y.Z]])
        G = M.T*M
        Gi = G.inv()
        bare = sp.Rational(1, 4)*((Gi*sp.diff(G, y.T))**2).trace() - sp.Rational(1, 4)*((Gi*sp.diff(G, y.Z))**2).trace()
        self.assertNotEqual(sp.simplify(y.law_from_frame(M) - bare), 0)


if __name__ == '__main__':
    unittest.main()
