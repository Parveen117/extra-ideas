"""TP1 tests (sympy)."""
import unittest

import sympy as sp

import tp1_frame_defect_law as y


class FrameDefectLawTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = y.run()

    def test_equivalence_selects_one_combination(self):
        tf = self.res['turned_frame']
        self.assertTrue(tf['exact_zero_on_1_2_minus4'])
        self.assertEqual(tf['exact_rank_at_identity'], 2)
        self.assertAlmostEqual(tf['null_direction'][1], 2.0, places=9)
        self.assertAlmostEqual(tf['null_direction'][2], -4.0, places=9)

    def test_other_combinations_depend_on_the_frame(self):
        T, X, Y, Z = sp.symbols('T X Y Z', real=True)
        al = sp.Function('alpha')(T, X, Y, Z)
        turn = sp.Matrix([[1, 0, 0, 0], [0, sp.cos(al), sp.sin(al), 0], [0, -sp.sin(al), sp.cos(al), 0], [0, 0, 0, 1]])
        I1, I2, I3 = y.invariants(y.structure_in((T, X, Y, Z), turn))
        from sympy.calculus.euler import euler_equations
        eq = euler_equations(I1, [al], [T, X, Y, Z])[0]
        self.assertNotEqual(sp.simplify(eq.lhs - eq.rhs), 0)

    def test_own_connection_is_flat_and_law_is_the_curvature_law(self):
        ff = self.res['fall_frame']
        self.assertTrue(ff['own_connection_flat'])
        self.assertEqual(ff['Q_plus_R_minus_k_div'][-2], '0')
        self.assertNotEqual(ff['Q_plus_R_minus_k_div'][2], '0')

    def test_profile_plane(self):
        ff = self.res['fall_frame']
        self.assertEqual(ff['profile_condition'], '9*r*(2*a1 + a2 + a3)/2')
        self.assertEqual(ff['profile_condition_at_selected'], '0')


if __name__ == '__main__':
    unittest.main()
