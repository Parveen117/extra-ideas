"""Independent source, truncation, transport and boundary controls for UP9."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import up9_haar_source_flow as u


class UP9Tests(unittest.TestCase):
    def test_moments_from_beta_integral_not_cumulant_recursion(self):
        for n in range(8):
            # (2/pi) integral_-1^1 c^(2n) sqrt(1-c^2) dc.
            beta_value = sp.simplify(2*sp.gamma(n+sp.Rational(1,2))
                *sp.gamma(sp.Rational(3,2))/(sp.pi*sp.gamma(n+2)))
            self.assertEqual(sp.Rational(u.haar(2*n)), beta_value)

    def test_source_centering_and_parity(self):
        self.assertEqual(u.mixed(0,1),0)
        for i in (1,3,5):
            for j in range(5):
                self.assertEqual(u.mixed(i,j),0)
        U=u.cumulant_potential()
        self.assertEqual(sp.expand(U.subs(u.S,-u.S)-U),0)
        self.assertNotEqual(sp.expand(U.subs(u.V,-u.V)-U),0)

    def test_gaussian_source_is_a_negative_control(self):
        # Its cumulants above degree two vanish; every member stays flat.
        S,V=u.S,u.V
        U=S*S/8+V*V/32
        m=sp.cancel(-V*sp.diff(U,S,V)/(S*sp.diff(U,S,2)))
        self.assertEqual(m,0)
        self.assertEqual(sp.diff(U,S,V),0)
        self.assertEqual(U-(S*sp.diff(U,S)+V*sp.diff(U,V))/2,0)

    def test_cut_flatness_does_not_mean_all_corner_frames_commute(self):
        # Selected cut defect begins at fourth order, but [V dV,D2]
        # already has a first-order coefficient V*dV(m)=-lambda V/4.
        m=u.leading_geometry()['m']
        self.assertEqual(sp.expand(u.V*sp.diff(m,u.V)).coeff(u.lam,1),-u.V/4)
        self.assertEqual(sp.expand(u.S*sp.diff(m,u.S)).coeff(u.lam,1),0)

    def test_sixth_cumulant_cannot_be_dropped(self):
        g=u.leading_geometry()
        bad=g['U']-g['U'].coeff(u.lam,4)*u.lam**4
        T=sp.diff(bad,u.S)
        bad_m=sp.series(-u.V*sp.diff(T,u.V)/(u.S*sp.diff(T,u.S)),u.lam,0,5).removeO()
        self.assertNotEqual(sp.factor(bad_m.coeff(u.lam,4)-g['m'].coeff(u.lam,4)),0)

    def test_cut_degeneracy_is_not_a_singular_source(self):
        g=u.leading_geometry()
        self.assertEqual(g['p0'].subs(u.V,2*u.S),0)
        h=sp.hessian(g['U'].subs(u.lam,0),(u.S,u.V))
        self.assertEqual(h.det(),sp.Rational(1,64))

    def test_outward_grid_and_interval_arithmetic(self):
        a=u.Interval.exact(Q(1,3));b=u.Interval.exact(Q(-7,11))
        for interval,answer in ((a+b,Q(1,3)-Q(7,11)),(a*b,-Q(7,33)),(a/b,-Q(11,21))):
            self.assertTrue(interval.contains(answer))
        self.assertLess(a.lo,Q(1,3));self.assertGreater(a.hi,Q(1,3))
        with self.assertRaises(ZeroDivisionError):
            u.Interval(Q(-1),Q(2)).reciprocal()

    def test_full_integral_tail_is_retained_and_shrinks(self):
        a=u.tilted_moment(0,0,Q(1,2),order=20)
        b=u.tilted_moment(0,0,Q(1,2),order=24)
        self.assertGreater(a.hi-a.lo,0)
        self.assertLessEqual(a.lo,b.lo);self.assertLessEqual(b.hi,a.hi)
        self.assertGreater(b.lo,0)
        with self.assertRaises(ValueError):u.tilted_moment(0,0,Q(1))
        with self.assertRaises(ValueError):u.haar(-1)

    def test_jet_transport_recovers_predecessor_exact_rational_witness(self):
        S,V=u.S,u.V
        potential=S**3/V+S**4/V+S**3/V**2
        jet={(i,j):u.Interval.exact(Q(sp.diff(potential,S,i,V,j).subs({S:1,V:1})
             /(sp.factorial(i)*sp.factorial(j))))
             for i in range(5) for j in range(5-i) if i+j}
        f=u.geometry_from_jet(jet)
        self.assertTrue(f['ell'].contains(Q(40000,45653449)))
        self.assertTrue(f['curvature'].contains(-Q(76489856000,2084237405595601)))

    def test_response_is_covariance_in_the_actual_tilted_measure(self):
        l=Q(1,4)
        z=u.tilted_moment(0,0,l)
        m=u.tilted_moment(1,0,l)/z
        variance=u.tilted_moment(2,0,l)/z-m*m
        got=u.finite_geometry(l)['uss']
        self.assertLessEqual(max(got.lo,variance.lo),min(got.hi,variance.hi))

    def test_source_selected_curvature_has_certified_printed_bounds(self):
        for l,low,high in ((Q(1,4),-Q(1229129,10**16),-Q(1229128,10**16)),
                          (Q(1,2),-Q(6081351,10**14),-Q(6081350,10**14))):
            f=u.finite_geometry(l)
            self.assertGreater(f['delta'].lo,0)
            self.assertLess(f['w'].hi,0)
            self.assertLessEqual(low,f['curvature'].lo)
            self.assertLessEqual(f['curvature'].hi,high)

    def test_aggregation_transports_curvature_as_a_two_form(self):
        S,V=u.S,u.V
        b=sp.symbols('b',positive=True)
        f9=u.leading_geometry()['f9']
        pull=sp.cancel(f9.subs({S:S/sp.sqrt(b),V:V/sp.sqrt(b)},simultaneous=True)/b)
        self.assertEqual(sp.simplify(pull-f9/b**sp.Rational(9,2)),0)
        # Missing the area Jacobian gives the wrong power by a factor b.
        self.assertNotEqual(sp.simplify(pull*b-f9/b**sp.Rational(9,2)),0)


if __name__=='__main__':
    unittest.main()
