"""NC1 tests (sympy)."""
import unittest

import sympy as sp

import nc1_native_curvature_of_the_field as y


class NativeCurvatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = y.run()

    def test_native_curvature_and_ratio(self):
        self.assertTrue(self.res['nt3_on_field'])
        self.assertEqual(self.res['ratio_of_field'], '-1/2')

    def test_variational_law_departs_at_second_order(self):
        self.assertEqual(self.res['kappa'], {'psi = 2 eta': '1/2', 'psi = eta': '-3/8'})
        self.assertNotEqual(self.res['variational_residual_of_ratio_field']['psi = 2 eta'], '0')

    def test_orbit_advance(self):
        m = self.res['mercury_arcsec_per_century']
        self.assertAlmostEqual(m['psi = 2 eta'], 57.3067, places=3)
        self.assertAlmostEqual(m['psi = eta'], 32.235, places=3)

    def test_another_ratio_is_another_field(self):
        """rho = -1 gives tanh(eta) = r0/r: memory falling as 1/r^2, not the observed 1/r."""
        r, r0 = sp.symbols('r r0', positive=True)
        t = r0/r
        deta = sp.diff(t, r)/(1 - t**2)
        self.assertEqual(sp.simplify(r*2*deta/(2*t/(1 - t**2))), -1)


if __name__ == '__main__':
    unittest.main()
