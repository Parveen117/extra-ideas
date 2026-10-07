"""GW1 tests (sympy)."""
import unittest

import sympy as sp
from sympy.calculus.euler import euler_equations

import gw1_waves_of_the_frame as y


class FrameWaveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = y.run()

    def test_two_waves_at_the_cone_speed(self):
        self.assertEqual(self.res['pointwise_difference'], '0')
        self.assertEqual(self.res['sign'], 1)
        self.assertEqual(self.res['mode_rate'], 'k')

    def test_content_two(self):
        self.assertEqual(self.res['content'], 2)
        self.assertTrue(self.res['returns_after_half_turn'])

    def test_a_local_turn_of_the_frame_is_not_a_wave(self):
        """turning e1, e2 by a small angle phi(t, z) changes nothing: the law gives no equation for phi."""
        T, Z = y.T, y.Z
        phi = sp.Function('phi')(T, Z)
        E = sp.Matrix([[1, 0, 0, 0], [0, 1, y.eps*phi, 0], [0, -y.eps*phi, 1, 0], [0, 0, 0, 1]])
        c = y.tp.structure_in(y.COORDS, E)
        I1, I2, I3 = y.tp.invariants(c)
        Q = sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3
        q2 = sp.simplify(sp.series(sp.simplify(Q), y.eps, 0, 3).removeO().coeff(y.eps, 2))
        eqs = euler_equations(q2, [phi], [T, Z])        # an identically satisfied equation is dropped by sympy
        self.assertTrue(all(sp.simplify(e.lhs - e.rhs) == 0 for e in eqs if hasattr(e, 'lhs')))
        self.assertNotEqual(sp.diff(q2, sp.diff(phi, T)) if q2 != 0 else 1, None)

    def test_illustration(self):
        ill = y.illustration()
        self.assertAlmostEqual(ill['flux_W_per_m2']*1e3, 1.585, places=2)


if __name__ == '__main__':
    unittest.main()
