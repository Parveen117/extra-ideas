"""Topology, normalization, discarded-channel and full-tail controls."""
from fractions import Fraction as Q
from itertools import combinations
from math import factorial
import unittest
import sympy as sp
import yc11_closed_cube_response as y


class YC11Tests(unittest.TestCase):
    def test_character_coefficients_from_direct_Haar_polynomials(self):
        c=sp.symbols('c')
        for n in range(1,6):
            for power in range(13):
                poly=sp.Poly(sp.expand(c**power*sp.chebyshevu(n-1,c)),c)
                direct=sum(Q(coeff)*y.u9.haar(k[0]) for k,coeff in poly.terms())
                self.assertEqual(direct/(n*factorial(power)),y.face_coefficient(n,power))

    def test_every_proper_face_subset_has_a_free_boundary(self):
        faces=y.cz.box_surface(1,1,1)
        edges={edge for f in faces for edge,sign in f}
        for count in range(1,6):
            for subset in combinations(range(6),count):
                hits=[sum(edge==ee for i in subset for ee,sgn in faces[i]) for edge in edges]
                self.assertIn(1,hits)
        self.assertTrue(all(sum(edge==ee for f in faces for ee,sgn in f)==2 for edge in edges))

    def test_marginals_are_three_independent_faces(self):
        for power in range(11):
            total=Q(0)
            for a in range(power+1):
                for b in range(power-a+1):
                    c=power-a-b
                    total+=Q(factorial(power),factorial(a)*factorial(b)*factorial(c))*y.u9.haar(a)*y.u9.haar(b)*y.u9.haar(c)
            self.assertEqual(y.moment(power,0),total)
            self.assertEqual(y.moment(0,power),total)

    def test_low_order_blindness_does_not_remove_the_joint_record(self):
        self.assertEqual(y.moment(1,1),0)
        self.assertEqual(y.moment(2,2),y.moment(2,0)*y.moment(0,2))
        self.assertEqual(y.moment(3,3),Q(9,256))
        self.assertEqual(y.moment(3,0)*y.moment(0,3),0)

    def test_complete_moments_obey_bounded_support(self):
        for a in range(11):
            for b in range(11-a):
                self.assertLessEqual(abs(y.moment(a,b)),3**(a+b))
                self.assertEqual(y.moment(a,b),y.moment(b,a))

    def test_no_closure_gives_a_flat_selected_compass(self):
        S,V=y.S,y.V
        U=y.potential()
        independent=U.subs(V,0)+U.subs(S,0)
        self.assertEqual(sp.diff(independent,S,V),0)
        self.assertNotEqual(sp.diff(U,S,V),0)
        self.assertEqual(sp.cancel(-V*sp.diff(independent,S,V)/(S*sp.diff(independent,S,2))),0)

    def test_angular_alignment_delays_the_curvature(self):
        g=y.leading_geometry()
        self.assertNotEqual(g['w4'],0)
        self.assertEqual(g['f12'],0)
        self.assertNotEqual(g['f14'],0)
        # A commuting scalar centre-weight algebra is a different object.
        self.assertEqual(sp.expand(g['ell8']*g['phi4']-g['phi4']*g['ell8']),0)

    def test_the_diagonal_cut_degeneracy_is_not_a_source_pole(self):
        g=y.leading_geometry()
        self.assertEqual(g['p0'].subs(y.V,y.S),0)
        self.assertEqual(sp.hessian(g['U'].subs(y.e,0),(y.S,y.V)).det(),sp.Rational(9,16))

    def test_boundary_tail_is_complete_and_decreases(self):
        k,h=Q(1),Q(1,2)
        self.assertGreater(y.interface_tail(k,h,6),0)
        self.assertLess(y.interface_tail(k,h,7),y.interface_tail(k,h,6))
        a,b=y.interface_interval(k,h,4),y.interface_interval(k,h,6)
        self.assertLessEqual(a.lo,b.lo);self.assertLessEqual(b.hi,a.hi)
        self.assertLess(b.hi,y.boundary_budget(k)*y.boundary_budget(h))

    def test_outward_log_and_denominator_guards(self):
        x=Q(1,2)
        got=y.log_one_plus(y.I.exact(x))
        lower=sum((Q((-1)**(n+1),n)*x**n for n in range(1,101)),Q(0))
        upper=lower+x**101/101
        self.assertLessEqual(got.lo,lower);self.assertGreaterEqual(got.hi,upper)
        with self.assertRaises(ValueError):y.log_one_plus(y.I.exact(1))
        with self.assertRaises(ValueError):y.boundary_budget(5)
        with self.assertRaises(ValueError):y.moment(-1,0)

    def test_full_integral_printed_curvature_and_centre_enclosures(self):
        cases=((Q(1,4),Q(1205337,10**19),Q(1205338,10**19),-Q(575901,10**10),-Q(575900,10**10)),
               (Q(1,2),Q(1283772,10**15),Q(1283773,10**15),-Q(778350,10**9),-Q(778349,10**9)))
        for l,lo,hi,mlo,mhi in cases:
            f=y.finite_geometry(l)['curvature']
            du,dm=y.centre_memory(l)
            self.assertLessEqual(lo,f.lo);self.assertLessEqual(f.hi,hi)
            self.assertGreater(du.lo,0)
            self.assertLessEqual(mlo,dm.lo);self.assertLessEqual(dm.hi,mhi)

    def test_full_source_and_commutator_readings_agree(self):
        for l in (Q(1,4),Q(1,2)):
            g=y.finite_geometry(l);a,b=g['w'],g['direct_w']
            self.assertLessEqual(max(a.lo,b.lo),min(a.hi,b.hi))

    def test_shared_kinetic_metric_survives_zero_Haar_covariance(self):
        # YC8's exact one-link formula: Gamma_e(Wp,Wq)=a.b-(a.u)(b.u).
        u=sp.Matrix([1,0,0,0]);a=sp.Matrix([0,1,0,0]);b=a
        self.assertEqual(a.dot(b)-a.dot(u)*b.dot(u),1)
        self.assertEqual(y.moment(1,1),0)

    def test_area_factor_in_source_dilation(self):
        b=sp.symbols('b',positive=True)
        f=y.leading_geometry()['f14']
        pulled=f.subs({y.S:y.S/sp.sqrt(b),y.V:y.V/sp.sqrt(b)},simultaneous=True)/b
        self.assertEqual(sp.simplify(pulled-f/b**7),0)


if __name__=='__main__':unittest.main()
