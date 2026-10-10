"""Frame covariance, independent matrix curvature and full-hidden controls."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import yc13_frame_flow as y


class YC13Tests(unittest.TestCase):
    def test_matrix_connection_independently_recovers_curvature_square(self):
        f = y.family()
        B = y.u8.cycle(f['L'])
        for xx, yy in ((Q(1, 8), Q(1, 32)), (Q(1, 4), Q(1, 4))):
            sub = {y.x: sp.Rational(xx), y.y: sp.Rational(yy)}
            b = B.subs(sub)
            X = [b.inv()*B.diff(t).subs(sub) for t in (y.x, y.y)]
            curvature = -(X[0]*X[1]-X[1]*X[0])/4
            expected = f['fxy'].subs(sub)
            self.assertEqual(sp.factor(-sp.trace(curvature*curvature)/2-expected**2), 0)

    def test_every_fixed_frame_is_flat_without_closing_its_cut(self):
        f = y.family()
        sub = {y.x: y.lam*y.V/y.S, y.y: y.lam*(y.V/y.S)**3}
        e = f['ell'].subs(sub)
        ph = {t: f['phix'].subs(sub)*sp.diff(sub[y.x], t)+f['phiy'].subs(sub)*sp.diff(sub[y.y], t)
              for t in (y.S, y.V)}
        self.assertTrue(y.zero(sp.diff(e, y.S)*ph[y.V]-sp.diff(e, y.V)*ph[y.S]))
        self.assertGreater(f['w'].subs({y.x: 1, y.y: 1}), 0)
        self.assertFalse(y.zero(y.u8.cycle(f['L']).subs({y.x: 1, y.y: 1})-sp.eye(2)))

    def test_interval_evaluation_contains_exact_point_curvatures(self):
        f = y.family()
        for lam, ratio in ((Q(0), Q(1, 2)), (Q(1, 4), Q(1)), (Q(1, 8), Q(3, 4))):
            exact = Q(2*sp.Rational(lam)*sp.Rational(ratio)**3*f['fxy'].subs(
                {y.x: sp.Rational(lam*ratio), y.y: sp.Rational(lam*ratio**3)}))
            interval = y.curvature_interval(y.I.exact(lam), y.I.exact(ratio))
            self.assertTrue(interval.contains(exact))

    def test_whole_loop_enclosures_overlap_and_exclude_zero(self):
        coarse, fine = y.loop_certificate(8), y.loop_certificate(16)
        self.assertLess(fine.hi, 0)
        self.assertLessEqual(max(coarse.lo, fine.lo), min(coarse.hi, fine.hi))
        self.assertLessEqual(fine.hi-fine.lo, coarse.hi-coarse.lo)

    def test_potential_is_not_its_corner_average(self):
        U = y.S**3/y.V+y.lam*(y.S*y.S+y.V*y.V)
        centre = U-(y.S*sp.diff(U, y.S)+y.V*sp.diff(U, y.V))/2
        self.assertTrue(y.zero(centre))
        self.assertFalse(y.zero(centre+y.lam*sp.diff(U, y.lam)/2))
        self.assertGreater(U.subs({y.S: 1, y.V: 1, y.lam: 0}), 0)

    def test_moving_observer_derivative_term_preserves_curvature(self):
        t, s = sp.symbols('t s', real=True)
        O = sp.Matrix([[1-t*t, -2*t], [2*t, 1-t*t]])/(1+t*t)*y.OBS
        a = [s*y.R, t*t*y.R]
        transformed = [O*a[i]*O.T-sp.diff(O, variable)*O.T for i, variable in enumerate((t, s))]
        F = sp.diff(transformed[1], t)-sp.diff(transformed[0], s)
        F += transformed[0]*transformed[1]-transformed[1]*transformed[0]
        self.assertTrue(y.zero(F-O*((2*t-1)*y.R)*O.T))
        self.assertTrue(y.zero(transformed[0]+s*y.R+2*y.R/(1+t*t)))

    def test_observation_preserves_finite_transition_composition(self):
        rot = lambda c, s: sp.Matrix([[c, -s], [s, c]])
        T1, T2 = rot(sp.Rational(3, 5), sp.Rational(4, 5)), rot(sp.Rational(5, 13), sp.Rational(12, 13))
        O0, O1, O2 = y.OBS, T1*y.OBS, T2*y.OBS
        observed1, observed2 = O1*T1*O0.T, O2*T2*O1.T
        self.assertEqual(observed2*observed1, O2*(T2*T1)*O0.T)
        self.assertEqual(observed1.T*observed1, sp.eye(2))

    def test_source_amplitude_null_direction_is_not_a_new_coordinate(self):
        k, h = y.lam*y.S, y.lam*y.V
        ell, phi = k*k+h, k*h*h
        variables = (y.lam, y.S, y.V)
        F = sp.Matrix([[sp.diff(ell, a)*sp.diff(phi, b)-sp.diff(ell, b)*sp.diff(phi, a)
                        for b in variables] for a in variables])
        self.assertTrue(y.zero(F*sp.Matrix([y.lam, -y.S, -y.V])))
        self.assertTrue(y.zero(F[0, 1]+y.V/y.lam*F[1, 2]))
        self.assertTrue(y.zero(F[0, 2]-y.S/y.lam*F[1, 2]))

    def test_full_hidden_step_retains_excursions_outside_source(self):
        D = sp.diag(19, 21, 26)
        T = sp.Matrix([[0, 0, 2], [0, 0, 3], [2, 3, 0]])
        source = sp.Matrix([[1, 0], [0, 1], [0, 0]])
        self.assertNotEqual(T*source, sp.zeros(3, 2))
        for zeta in (0, 11):
            sigma = lambda t: source.T*(D-sp.Rational(t)*T-zeta*sp.eye(3)).inv()*source
            for t1, t2 in ((Q(0), Q(1, 2)), (Q(1, 2), Q(1))):
                a, b = sigma(t1), sigma(t2)
                upper = sp.Rational(y.return_step(t1, t2, zeta))*a
                for difference in (b-a, upper-b):
                    self.assertEqual(y.y12.y4.inertia(difference.tolist())[0], 0)
                self.assertGreater(b[0, 0], a[0, 0])

    def test_parity_hypothesis_is_essential_for_even_return(self):
        theta = sp.symbols('theta')
        bad_return = 1/(20-2*theta)
        self.assertNotEqual(sp.diff(bad_return, theta).subs(theta, 0), 0)
        D = sp.diag(19, 26)
        T = sp.Matrix([[0, 2], [2, 0]])
        good = (D-theta*T).inv()[0, 0]
        self.assertEqual(sp.diff(good, theta).subs(theta, 0), 0)
        self.assertGreater(sp.diff(good, theta, 2).subs(theta, 0), 0)

    def test_full_window_face_certificate_interior_controls(self):
        for tmax, zeta in ((Q(1), Q(11)), (Q(3, 2), Q(42, 5)), (Q(2), Q(27, 5))):
            for k in range(5):
                theta = tmax*k/4
                pencil = y.y12.pencil(theta, zeta, True)
                self.assertEqual(y.y12.y4.inertia(pencil), (1, 0, 6))
                self.assertEqual(y.y12.y4.inertia([row[1:] for row in pencil[1:]]), (0, 0, 6))

    def test_bad_return_and_gap_reserves_are_rejected(self):
        for args in ((1, 0, 0), (0, 1, 12), (0, 1, 18), (-1, 1, 0)):
            with self.assertRaises(ValueError):
                y.return_step(*args)
        for args in ((-1, 1), (1, 12), (1, 0), (1, Q(23, 2))):
            with self.assertRaises(ValueError):
                y.gap_window(*args)

    def test_cube_gauge_floor_is_not_reused_for_charged_blocks(self):
        self.assertEqual(y.gap_window(2, Q(27, 5))['gauge_gap_lower'], Q(27, 5))
        self.assertEqual(y.y12.join_bounds(0)['g'], Q(1, 4))
        with self.assertRaises(ValueError):
            y.y12.cube_bounds(2)

    def test_frozen_actual_source_mixed_components_keep_full_errors(self):
        records = y.cube_mixed_certificates()
        for row in records:
            state = tuple(Q(v) for v in row['state'])
            a = tuple(Q(v) for v in row['mixed_lambda_S'])
            b = tuple(Q(v) for v in row['mixed_lambda_V'])
            self.assertLess(a[1], 0)
            self.assertGreater(b[0], 0)
            l, S, V = row['lambda_value'], row['S'], row['V']
            self.assertLessEqual(a[0], -V/l*state[1])
            self.assertGreaterEqual(a[1], -V/l*state[0])
            self.assertLessEqual(b[0], S/l*state[0])
            self.assertGreaterEqual(b[1], S/l*state[1])


if __name__ == '__main__':
    unittest.main()
