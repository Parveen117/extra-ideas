"""RV1: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import rv1_spin_and_strain_of_the_river as y


class SpinAndStrainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = y.run()

    def test_pinned_statements(self):
        self.assertEqual(self.res['rigid_turn'], {'I1': '4 Omega^2', 'I2': '-2 Omega^2', 'I3': '0', 'law': '0'})
        self.assertEqual(self.res['axial_river_exponents'], 'omega ~ r^n with n(n+3) = 0: n = 0 (rigid) or n = -3')
        self.assertAlmostEqual(self.res['numbers']['gyroscope_polar_orbit_7020km_mas_per_year'], 40.9, places=1)

    def test_other_combinations_of_the_invariants_do_see_the_spin(self):
        Om = sp.Symbol('Omega', positive=True)
        I1, I2, I3 = y.tp.invariants(y.tp.structure_in(y.COORDS, y.river_frame([-Om*y.y, Om*y.x, 0])))
        self.assertNotEqual(sp.simplify(I1 + I2 + I3), 0)
        self.assertEqual(sp.simplify(sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3), 0)

    def test_expansion_and_shear_are_seen(self):
        h = sp.Symbol('h', positive=True)
        self.assertEqual(sp.simplify(y.law_from_frame([h*y.x, h*y.y, h*y.z]) + 6*h**2), 0)       # CO1's value
        self.assertEqual(sp.simplify(y.law_from_frame([h*y.y, 0, 0]) - h**2/2), 0)               # simple shear = strain + spin

    def test_differential_turn_is_not_silent(self):
        rr = sp.sqrt(y.x**2 + y.y**2 + y.z**2)
        q = sp.simplify(y.law_from_strain([-y.y/rr**3, y.x/rr**3, 0]))
        self.assertNotEqual(q, 0)
        self.assertEqual(sp.simplify(q.subs({y.x: 0, y.y: 0})), 0)                                # vanishes on the axis

    def test_a_law_with_a_spin_term_is_rejected(self):
        def wrong(v):
            S, W = y.strain_spin(v)
            return sum(S[i, j]**2 + W[i, j]**2 for i in range(3) for j in range(3)) - S.trace()**2
        with patch.object(y, 'law_from_strain', wrong):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
