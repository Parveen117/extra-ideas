"""Controls for lost sheets, normalization and predecessor compatibility."""
from fractions import Fraction as F
import importlib.util
import unittest
import sympy as sp
import hc1_signed_horizon as h


class HC1Tests(unittest.TestCase):
    def test_response_radius_is_not_fall_speed_or_static_clock(self):
        u, v = h.signed_reading(F(1, 2))
        self.assertEqual((u, v), (F(4, 5), F(3, 5)))
        self.assertNotEqual(u, F(1, 2))
        self.assertNotEqual(1-u*u, F(3, 4))
        self.assertEqual(h.exterior_curvatures(4)['clock_squared'], F(3, 4))

    def test_sheet_loss_is_a_causal_counterexample(self):
        outside, inside = F(2, 3), F(3, 2)
        u1, v1 = h.signed_reading(outside)
        u2, v2 = h.signed_reading(inside)
        self.assertEqual(u1, u2)
        self.assertEqual(v1, -v2)
        self.assertGreater(h.radial_speeds(outside)[0], 0)
        self.assertLess(h.radial_speeds(inside)[0], 0)

    def test_clock_alone_does_not_determine_null_speed(self):
        N = F(3, 5)
        speeds = []
        for A in (N, 1/N):
            c = N/A
            self.assertEqual(N*N-A*A*c*c, 0)
            speeds.append(c)
        self.assertEqual(speeds, [1, F(9, 25)])

    def test_angular_motion_cannot_restore_an_outward_interior_ray(self):
        # Pythagorean local null velocities: radial w with transverse 1-w^2.
        for beta in (F(1, 2), F(1), F(3, 2)):
            for w in (F(-1), F(-3, 5), F(0), F(3, 5), F(1)):
                dr, angular_squared = w-beta, 1-w*w
                self.assertEqual((dr+beta)**2+angular_squared, 1)
                self.assertLessEqual(dr, h.radial_speeds(beta)[0])
                if beta > 1: self.assertLess(dr, 0)

    def test_nonunit_radius_normalization_and_weight(self):
        r, rs, scale = F(7, 3), F(2, 3), F(11, 5)
        a = h.exterior_curvatures(r, rs)
        b = h.exterior_curvatures(scale*r, scale*rs)
        for key in ('radial_angular', 'angular_angular', 'native_diagnostic', 'geometric_square'):
            self.assertEqual(b[key], a[key]/scale**4)
        self.assertEqual(a['native_diagnostic']*a['clock_squared']**4, a['geometric_square']/2)

    def test_positive_response_is_not_silently_extended_inside(self):
        for r, rs in ((1, 1), (F(1, 2), 1), (2, 0)):
            with self.assertRaises(ValueError): h.exterior_curvatures(r, rs)
        self.assertEqual(h.radial_speeds(F(2)), (-1, -3))

    def test_existing_lt1_generation_and_native_contraction(self):
        spec = importlib.util.spec_from_file_location('hc1_lt1',
            h.ROOT/'physics/lt1/lt1_lambda_tower_inward.py')
        lt = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(lt)
        for beta in (F(1, 7), F(2, 3), F(7, 5)):
            u, _ = h.signed_reading(beta)
            self.assertEqual(lt.square_m(beta*beta), u*u)
            self.assertEqual(lt.radius_map(1/beta**2, F(1)), 1/u**2)
        # Numerical-free substitution into NC1's derivative formula, outside.
        r, rs = F(5), F(4, 5)  # beta=2/5
        beta = F(2, 5)
        psi_prime, sinh = -beta/(r*(1-beta*beta)), 2*beta/(1-beta*beta)
        self.assertEqual(-(psi_prime*sinh)**2/(4*r*r), h.exterior_curvatures(r, rs)['radial_angular'])

    def test_horizon_is_regular_for_a_falling_orthonormal_frame(self):
        # Finite coframe and inverse at the horizon, unlike the static chart.
        coframe = sp.Matrix([[1, 0], [1, 1]])
        pairing = sp.diag(1, -1)
        metric = coframe.T*pairing*coframe
        self.assertEqual(metric.det(), -1)
        self.assertEqual(metric, sp.Matrix([[0, -1], [-1, -1]]))
        fall = sp.Matrix([1, -1])
        self.assertEqual(coframe*fall, sp.Matrix([1, 0]))
        self.assertEqual((fall.T*metric*fall)[0], 1)

    def test_field_squaring_is_not_radial_resampling(self):
        beta, r = F(1, 2), F(4)
        c = (1+beta*beta)/(1-beta*beta)
        original = h.exterior_curvatures(r)
        fixed_base_radial = 16*c*c*original['radial_angular']
        fixed_base_angular = 16*c**4*original['angular_angular']
        fixed_base_D = -(2*fixed_base_radial+fixed_base_angular)
        next_beta, _ = h.signed_reading(beta)
        resampled = h.exterior_curvatures(1/next_beta**2)
        self.assertGreater(fixed_base_D, 16*original['native_diagnostic'])
        self.assertNotEqual(fixed_base_D, resampled['native_diagnostic'])
        self.assertNotEqual(fixed_base_radial, resampled['radial_angular'])
        self.assertNotEqual(fixed_base_angular/fixed_base_radial, 4)


if __name__ == '__main__':
    unittest.main()
