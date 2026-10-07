"""CV1 tests (sympy)."""
import unittest

import sympy as sp

import cv1_frame_curvature as y


class FrameCurvatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = y.run()

    def test_contracted_curvature_and_vacuum(self):
        self.assertEqual(self.res['overall_sign'], -1)
        self.assertEqual(self.res['vacuum_tidal'], ['r_s/r**3', '-r_s/(2*r**3)', '-r_s/(2*r**3)'])
        self.assertEqual(self.res['curvature_square'], '12*r_s**2/r**6')

    def test_constant_memory_is_not_curvature_free(self):
        self.assertEqual(self.res['constant_memory_contracted'][2], 'B/r**2')
        self.assertEqual(self.res['mc1_operator_on_A_over_r_plus_B'], '0')

    def test_constant_curvature_case(self):
        self.assertEqual(self.res['constant_curvature_ratio'], '-L')

    def test_order_defect_of_the_frame_is_the_gradient_of_the_fall(self):
        c = self.res['structure_nonzero']
        self.assertEqual(c['011'], 'Derivative(beta(r), r)')
        self.assertEqual(c['022'], 'beta(r)/r')

    def test_a_frame_at_rest_has_no_curvature(self):
        E = y.frame(sp.Integer(0))
        c = y.structure(E)
        R = y.curvature(E, c, y.connection(c))
        self.assertTrue(all(sp.simplify(v) == 0 for v in R.values()))

    def test_a_uniform_fall_has_curvature(self):
        """a constant fall speed is not a change of frame: the angular defect beta/r remains."""
        b = sp.Rational(1, 3)
        E = y.frame(b)
        c = y.structure(E)
        R = y.curvature(E, c, y.connection(c))
        self.assertTrue(any(sp.simplify(v) != 0 for v in R.values()))


if __name__ == '__main__':
    unittest.main()
