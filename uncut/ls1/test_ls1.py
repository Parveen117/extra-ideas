"""Independent geometric, atlas, covariance and scope controls."""
import hashlib
import json
import unittest
import sympy as sp
import ls1_compass_space as l


class LS1Tests(unittest.TestCase):
    def test_zero_seam_is_not_the_return_point(self):
        self.assertNotEqual(l.K, l.J)
        self.assertEqual((l.R*l.K)**2, l.I)
        self.assertEqual((l.R*l.J)**2, l.I)
        self.assertEqual(l.chord_squared(l.K, l.J), 2)
        self.assertEqual(l.chord_squared(l.K, -l.K), 4)

    def test_rational_chart_positive_return_and_inverse(self):
        for b in (sp.Rational(1, 3), 1, 2, 5):
            for c, s in ((1, 0), (sp.Rational(3, 5), sp.Rational(4, 5)), (0, 1)):
                k = l.rational_cut(b, c, s)
                ret = (l.R*k)**2
                self.assertEqual(k*k, l.I)
                self.assertEqual(ret, ret.T)
                self.assertEqual(ret.det(), 1)
                self.assertGreater(ret[0, 0], 0)
                self.assertEqual(ret.trace(), sp.Rational(b)**2+sp.Rational(b)**-2)

    def test_pairing_covariance_requires_transporting_the_pairing(self):
        z = sp.Matrix([[2, 1], [0, sp.Rational(1, 2)]])
        x = sp.Matrix([[0, 2], [3, 0]])
        g = sp.diag(2, 5)
        xp, gp = z*x*z.inv(), z.inv().T*g*z.inv()
        original = sp.trace(g.inv()*x.T*g*x)/2
        transported = sp.trace(gp.inv()*xp.T*gp*xp)/2
        self.assertEqual(original, transported)
        self.assertNotEqual(sp.trace(x.T*x), sp.trace(xp.T*xp))

    def test_gaussian_curvature_via_independent_christoffel_calculation(self):
        r = l.r
        e, g = (1+2*r*r)/(1+r*r), 1+r*r
        # R^r_(phi r phi) from the three nonzero Christoffel entries.
        a = sp.diff(e, r)/(2*e)
        b = -sp.diff(g, r)/(2*e)
        c = sp.diff(g, r)/(2*g)
        riemann = sp.diff(b, r)+a*b-b*c
        self.assertEqual(sp.cancel(riemann/g-l.gauss_curvature(r)), 0)

    def test_metric_bounds_on_unit_strip(self):
        r = l.r
        e, g = l.metric()[0, 0], l.metric()[1, 1]
        self.assertEqual(sp.cancel(e-1-r*r/(1+r*r)), 0)
        self.assertEqual(sp.cancel(2-e-1/(1+r*r)), 0)
        self.assertEqual(g-1, r*r)
        for x in (0, sp.Rational(1, 4), sp.Rational(3, 4), 1):
            self.assertTrue(1 <= e.subs(r, x) <= 2)
            self.assertTrue(1 <= g.subs(r, x) <= 2)

    def test_curvature_pullback_can_vanish_on_rank_one_source(self):
        x, y = sp.symbols('x y', real=True)
        depth, angle = x*x, 3*x
        pulled = 2*depth*(sp.diff(depth, x)*sp.diff(angle, y)-sp.diff(depth, y)*sp.diff(angle, x))
        self.assertEqual(pulled, 0)
        angle = y
        pulled = 2*depth*(sp.diff(depth, x)*sp.diff(angle, y)-sp.diff(depth, y)*sp.diff(angle, x))
        self.assertEqual(pulled, 4*x**3)

    def test_annulus_curvature_matches_loop_return_difference(self):
        a, b, r, phi = sp.symbols('a b r phi', real=True)
        integral = sp.integrate(2*r, (r, a, b), (phi, 0, 2*sp.pi))
        self.assertEqual(sp.expand(-integral+2*sp.pi*(b*b-a*a)), 0)
        # At depth1 endpoint rotation is identity but unwrapped angle is -2pi.
        theta = -2*sp.pi
        rot = sp.cos(theta)*l.I+sp.sin(theta)*l.R
        self.assertEqual(rot, l.I)
        self.assertNotEqual(theta, 0)

    def test_zero_observer_reset_sign_after_one_winding(self):
        p = l.phi
        c = sp.cos(p/2)*l.I-sp.sin(p/2)*l.R
        self.assertTrue(l.zero(c*l.chart(0, p)*c.T-l.K))
        self.assertEqual(c.subs(p, 0), l.I)
        self.assertEqual(c.subs(p, 2*sp.pi), -l.I)

    def test_analytic_at_one_with_a_different_chart_expansion(self):
        x = sp.symbols('x', real=True)
        f = sp.sqrt(1+(1+x)**2)
        self.assertEqual(f.subs(x, 0), sp.sqrt(2))
        self.assertEqual(sp.diff(f, x).subs(x, 0), sp.sqrt(2)/2)
        series = sp.series(sp.sqrt(1+x*x), x, 0, 8).removeO()
        self.assertEqual(series, 1+x*x/2-x**4/8+x**6/16)
        self.assertEqual((1+2*x*x).subs(x, sp.I/sp.sqrt(2)), 0)

    def test_tower_level_labels_do_not_select_an_interpolating_path(self):
        x = sp.symbols('x', real=True)
        for n in range(6):
            self.assertEqual(sp.sin(sp.pi*x).subs(x, n), 0)
        self.assertEqual(sp.sin(sp.pi*x).subs(x, sp.Rational(1, 2)), 1)

    def test_scope_and_source_pins(self):
        record = json.loads((l.HERE/'LS1_RESULT.json').read_text())
        for path, digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256((l.ROOT/path).read_bytes()).hexdigest(), digest)
        self.assertTrue(all(record['checks'].values()))
        self.assertFalse(record['claims']['full_primitive_lambda_space_selected'])
        self.assertFalse(record['claims']['new_yang_mills_gap'])


if __name__ == '__main__':
    unittest.main()
