"""Independent source, jet, continuation and recognition controls for LS2."""
import hashlib
import json
import unittest

import sympy as sp

import ls2_response_tower as l


class LS2Tests(unittest.TestCase):
    def test_first_coefficients_against_direct_corner_response_derivatives(self):
        source = l.up8.family()['L'].subs({l.x: 2*l.t, l.y: l.t/3})
        coefficients = l.response_coefficients(2, sp.Rational(1, 3), 3)
        for n in range(4):
            derivative = source.diff(l.t, n).subs(l.t, 0)/sp.factorial(n)
            self.assertTrue(l.zero(derivative-coefficients[n]))

    def test_full_source_includes_derivatives_of_the_constrained_operator(self):
        # A potential affine in lambda has a nonterminating response tower.
        source = l.up8.family()['L'].subs({l.x: l.t, l.y: l.t})
        self.assertNotEqual(source[1, 1].diff(l.t, 6).subs(l.t, 0), 0)
        self.assertEqual(source[0, 0].diff(l.t, 2), 0)

    def test_mean_and_size_are_retained_in_the_lift(self):
        coefficients = [3*l.I+5*l.K, 7*l.I+2*l.K+3*l.J+l.R,
                        -2*l.I+l.K-l.J]
        means, sizes, cuts = l.normalize_jets(coefficients)
        self.assertEqual(means, [3, 7, -2])
        self.assertEqual(sizes[0], 5)
        for n in range(3):
            self.assertTrue(l.zero(means[n]*l.I+l.matrix_sum(sizes[j]*cuts[n-j]
                                                           for j in range(n+1))-coefficients[n]))

    def test_lift_accepts_a_cut_with_an_elliptic_first_derivative(self):
        coefficients = [l.K, l.R, l.K/2, sp.zeros(2), -l.K/8]
        means, sizes, cuts = l.normalize_jets(coefficients)
        self.assertEqual(means, [0]*5)
        self.assertEqual(sizes, [1, 0, 0, 0, 0])
        self.assertEqual(cuts, coefficients)
        self.assertEqual(cuts[1]**2, -l.I)

    def test_nonpositive_base_discriminant_is_rejected(self):
        for matrix in (l.I, l.R, l.K+l.R):
            with self.assertRaises(ValueError):
                l.normalize_jets([matrix])

    def test_exact_reconstruction_on_three_source_rays_beyond_first_disk(self):
        for alpha, beta, value in ((1, 1, 1), (2, 1, 2), (3, 2, 1)):
            coefficients = l.response_coefficients(alpha, beta, 4)
            rebuilt = l.reconstruct_from_five(coefficients, 2*alpha+beta)
            target = l.up8.at(l.up8.family()['L'], alpha*value, beta*value)
            self.assertTrue(l.zero(rebuilt.subs(l.t, value)-target))

    def test_known_denominator_is_a_necessary_model_input(self):
        coefficients = l.response_coefficients(1, 1, 4)
        wrong = l.reconstruct_from_five(coefficients, 4)
        true = l.up8.family()['L'].subs({l.x: l.t, l.y: l.t})
        # Both match the supplied finite packet, but not its continuation.
        self.assertTrue(l.zero(wrong.subs(l.t, 0)-true.subs(l.t, 0)))
        self.assertFalse(l.zero(wrong.subs(l.t, 1)-true.subs(l.t, 1)))

    def test_state_decoder_on_exact_positive_inputs(self):
        r1, phi1 = l.first_compass_jets()
        for ss, vv in ((1, 1), (sp.Rational(2, 3), 5), (7, sp.Rational(1, 4))):
            rr, pp = [z.subs({l.S: ss, l.V: vv}) for z in (r1, phi1)]
            self.assertEqual(-37*rr/4, ss)
            self.assertEqual(420/(1369*pp+652*ss), vv)

    def test_source_metric_leading_term_from_raw_response_differentials(self):
        # Independently differentiate w/sqrt(Delta) and arg(p+iq) in x,y.
        f = l.up8.family()
        origin = {l.x: 0, l.y: 0}
        p0, q0, d0 = [f[name].subs(origin) for name in ('p', 'q', 'delta')]
        dr = [f['w'].diff(z).subs(origin)/sp.sqrt(d0) for z in (l.x, l.y)]
        dp = [(p0*f['q'].diff(z)-q0*f['p'].diff(z)).subs(origin)/(p0*p0+q0*q0)
              for z in (l.x, l.y)]
        # dx/lambda=dS, dy/lambda=-dV/V^2 on the state slice.
        jac = sp.Matrix([dr, dp])*sp.diag(1, -1/l.V**2)
        expected = sp.Matrix(l.first_compass_jets()).jacobian((l.S, l.V))
        self.assertTrue(l.zero(jac.T*jac-expected.T*expected))

    def test_curvature_order_needs_constant_zero_angle(self):
        t, x, y = l.t, l.x, l.y
        depth, angle = t*x, y
        curvature = 2*depth*(depth.diff(x)*angle.diff(y)-depth.diff(y)*angle.diff(x))
        self.assertEqual(curvature, 2*t*t*x)
        # A varying zero angle changes the leading state-surface order to two.

    def test_one_parameter_restriction_cannot_have_two_form_curvature(self):
        x, y = l.x, l.y
        depth, angle = x*x, 3*x
        self.assertEqual(depth.diff(x)*angle.diff(y)-depth.diff(y)*angle.diff(x), 0)
        self.assertNotEqual(depth.diff(x)**2+angle.diff(x)**2, 0)

    def test_error_gate_with_exact_rational_perturbations(self):
        v = sp.Matrix([3, 4, 1])
        d = 24
        M = 6  # ||v||=sqrt(26)<6.
        for e in (sp.Matrix([sp.Rational(1, 10), 0, 0]),
                  sp.Matrix([0, 0, sp.Rational(1, 10)])):
            epsilon = sp.Rational(1, 10)
            lower = d-2*M*epsilon-epsilon**2
            new = v+e
            delta = new[0]**2+new[1]**2-new[2]**2
            self.assertGreater(lower, 0)
            self.assertGreaterEqual(delta, lower)

    def test_source_pins_and_claim_boundaries(self):
        record = json.loads((l.HERE/'LS2_RESULT.json').read_text())
        for path, digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256((l.ROOT/path).read_bytes()).hexdigest(), digest)
        self.assertTrue(all(record['checks'].values()))
        self.assertTrue(record['claims']['first_jet_repairs_two_state_directions'])
        self.assertFalse(record['claims']['lambda_is_physical_RG_scale'])
        self.assertFalse(record['claims']['new_yang_mills_gap'])


if __name__ == '__main__':
    unittest.main()
