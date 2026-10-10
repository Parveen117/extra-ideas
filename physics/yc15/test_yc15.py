"""Independent source, residual, channel and full-return controls for YC15."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import yc15_boundary_return_channels as y


def signed_axes():
    return [tuple(sign*int(i == j) for i in range(4))
            for j in range(4) for sign in (-1, 1)]


class YC15Tests(unittest.TestCase):
    def test_signed_bridge_cubature_keeps_norm_and_both_link_correlations(self):
        # Exact on degree<=2 in each independent quaternion bridge.
        a, b = (Q(3, 5), Q(4, 5), 0, 0), (Q(5, 13), 0, Q(12, 13), 0)
        aa, bb = (0, 0, 0, 1), (Q(4, 5), Q(3, 5), 0, 0)
        axes = signed_axes()
        mean = lambda fn: sum(fn(u, v) for u in axes for v in axes)/64
        self.assertEqual(mean(lambda u, v: y.wilson(a, b, u, v)), 0)
        self.assertEqual(mean(lambda u, v: y.wilson(a, b, u, v)**2), Q(1, 4))
        self.assertEqual(mean(lambda u, v: y.wilson(a, b, u, v)*y.wilson(aa, bb, u, v)),
                         y.dot(a, aa)*y.dot(b, bb)/4)
        self.assertEqual(mean(lambda u, v: y.wilson(a, b, u, v)*y.wilson(b, a, u, v)),
                         y.dot(a, b)**2/4)

    def test_conditional_residual_variance_has_no_cross_cube_term(self):
        a, b = (Q(3, 5), Q(4, 5), 0, 0), (Q(5, 13), 0, Q(12, 13), 0)
        da = y.quaternion_product((0, 2, 0, 0), a)
        db = y.quaternion_product((0, 0, Q(1, 3), 0), b)
        self.assertEqual(y.dot(a, da), 0)
        self.assertEqual(y.dot(b, db), 0)
        direct = sum((-2*(y.wilson(da, b, u, v)+y.wilson(a, db, u, v)))**2
                     for u in signed_axes() for v in signed_axes())/64
        self.assertEqual(direct, y.dot(da, da)+y.dot(db, db))
        mixed = sum(y.wilson(da, b, u, v)*y.wilson(a, db, u, v)
                    for u in signed_axes() for v in signed_axes())/64
        self.assertEqual(mixed, 0)

    def test_kinetic_bound_and_source_variance_cover_the_whole_window(self):
        previous = Q(0)
        for i in range(17):
            theta = Q(i, 16)
            value = y.edge_kinetic_bound(theta)
            self.assertGreaterEqual(value, previous)
            self.assertLessEqual(value, Q(1, 24))
            previous = value
            row = y.source_bounds(theta, 1-theta)
            self.assertLessEqual(row['variance_upper'], Q(1, 12))
            self.assertLessEqual(row['inverse_vector_norm_squared_upper'], y.INVERSE_NORM_CAP**2)

    def test_complete_inverse_residual_identity_not_a_source_span_closure(self):
        r = sp.Matrix([sp.Rational(1, 2), 0])
        for off in (sp.Rational(1, 8), sp.Rational(1, 3), sp.Rational(1, 2)):
            D = sp.Matrix([[12, off], [off, 12]])
            s = (D-12*sp.eye(2))*r
            actual = (r.T*D.inv()*r)[0]
            residual = (s.T*D.inv()*s)[0]/144
            self.assertEqual(r.dot(s), 0)
            self.assertEqual(actual, Q(1, 48)+residual)
            self.assertGreater(actual, Q(1, 48))
            self.assertNotEqual((D.inv()*r)[1], 0)
            self.assertLessEqual(actual, y.source_bounds(1, 1)['susceptibility_upper'])
            self.assertLessEqual((D.inv()*r).dot(D.inv()*r),
                                 y.source_bounds(1, 1)['inverse_vector_norm_squared_upper'])

    def test_residual_remainder_bound_is_quadratic_in_internal_coupling(self):
        zero = y.source_bounds(0, 0)
        self.assertEqual(zero['susceptibility_upper'], zero['susceptibility_lower'])
        self.assertEqual(zero['inverse_vector_norm_squared_upper'], Q(1, 576))
        self.assertEqual(zero['cube_return_difference_norm_squared_upper'], 0)
        half, full = y.source_bounds(Q(1, 2), Q(1, 2)), y.source_bounds(1, 1)
        self.assertEqual(full['variance_upper']/half['variance_upper'], 9)
        self.assertEqual(full['susceptibility_upper']-Q(1, 48), Q(1, 11232))

    def test_common_edge_product_resolves_into_degree_zero_and_two(self):
        # Direct sphere Laplacian on u_i W, keeping the other three links'
        # energy9. No use of the heat-response coefficient implementation.
        u = sp.symbols('u0:4')
        a = sp.symbols('a0:4')
        W = sum(x*z for x, z in zip(u, a))
        for i in range(4):
            polynomial = u[i]*W
            euler = lambda f: sum(u[j]*sp.diff(f, u[j]) for j in range(4))
            sphere_H = euler(euler(polynomial))+2*euler(polynomial)-sum(sp.diff(polynomial, z, 2) for z in u)
            actual = sp.expand(sphere_H+9*polynomial)
            expected = 9*a[i]/4+17*(polynomial-a[i]/4)
            self.assertEqual(sp.expand(actual-expected), 0)

    def test_laplace_and_stationary_derivatives_independently_match(self):
        t = sp.symbols('t', nonnegative=True)
        f = sp.exp(-3*t)/84-sp.exp(-9*t)/48+sp.exp(-17*t)/112
        exact_integral = sp.integrate(sp.exp(-9*t)*f/4, (t, 0, sp.oo))
        self.assertEqual(exact_integral, Q(1, 22464))
        self.assertEqual(exact_integral, y.first_channel_coefficient(True))
        self.assertEqual(y.first_channel_coefficient(False), 0)

    def test_omitting_reference_ground_change_gives_wrong_channel(self):
        raw = Q(1, 6)*(Q(1, 16*18)+Q(3, 16*26))
        self.assertNotEqual(raw, Q(1, 22464))
        self.assertEqual(raw-Q(1, 576), Q(1, 22464))
        # A disjoint face has no return after that same subtraction.
        self.assertEqual(Q(1, 6*4*24)-Q(1, 576), 0)

    def test_incident_face_response_positive_but_not_instantaneous(self):
        z = sp.symbols('z')
        polynomial = sp.div(4-7*z**3+3*z**7, (1-z)**2)[0]
        self.assertTrue(all(coefficient > 0 for coefficient in sp.Poly(polynomial, z).all_coeffs()))
        t = sp.symbols('t')
        f = sp.exp(-3*t)/84-sp.exp(-9*t)/48+sp.exp(-17*t)/112
        self.assertEqual(sp.series(f, t, 0, 3).removeO(), t*t/2)

    def test_each_cube_edge_has_two_first_return_faces(self):
        edges, faces, _, _ = y.y12.cube()
        for e in range(len(edges)):
            active = [face for face in faces if e in face]
            self.assertEqual(len(active), 2)
            self.assertEqual(active[0] & active[1], frozenset((e,)))

    def test_whole_interface_window_has_both_nonlinear_reserves(self):
        for i in range(13):
            row = y.join_bounds(y.ETA_MAX*i/12)
            self.assertGreater(row['mapping_reserve'], 0)
            self.assertLess(row['contraction'], 1)
            self.assertLess(row['relative_return'], 1)
            self.assertGreater(row['full_gap_lower'], Q(9, 40))
            self.assertEqual(row['full_gap_lower'], y.GAP-Q(1024, 5)*row['eta'])
            self.assertEqual(row['nonlinear_minimum_support'], 1)

    def test_both_source_and_charged_gap_refinements_are_needed(self):
        row = y.join_bounds(y.ETA_MAX)
        old_seed = Q(96, 13)*y.ETA_MAX
        self.assertGreater(old_seed+row['contraction']*y.RADIUS, y.RADIUS)
        old_q = Q(5184, 5)*y.ETA_MAX/Q(1, 4)
        self.assertGreater(row['seed']+old_q*y.RADIUS, y.RADIUS)
        self.assertLess(row['fixed_point_norm_upper'], y.RADIUS)

    def test_boundary_signature_rule_survives_periodic_seams(self):
        for shape in ((4, 4, 4), (4, 6, 8)):
            row = y.y14.geometry_certificate(shape)
            self.assertTrue(row['bridge_signatures_unique'])
            self.assertEqual(row['cube_profiles'], [(24, 0)])
            self.assertEqual(row['bridge_profiles'], [(2, 2)])

    def test_invalid_windows_rejected(self):
        for theta in (-1, Q(1001, 1000)):
            with self.assertRaises(ValueError):
                y.source_bounds(theta, 0)
        for eta in (-1, y.ETA_MAX+Q(1, 10**9)):
            with self.assertRaises(ValueError):
                y.join_bounds(eta)

    def test_frozen_input_pins_match(self):
        for path in ('physics/yc12/YC12_RESULT.json', 'physics/yc14/YC14_RESULT.json'):
            y.verify_predecessor(path)


if __name__ == '__main__':
    unittest.main()
