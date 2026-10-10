"""Negative controls and independent checks of UP8's source/transport map."""
import importlib.util
import unittest
import sympy as sp
import up8_deformed_compass as u


class UP8Tests(unittest.TestCase):
    def test_closed_limit_matches_existing_cut_carrier(self):
        spec=importlib.util.spec_from_file_location('up8_up7',u.ROOT/'uncut/up7/up7_turn_and_cut.py')
        old=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(old)
        self.assertEqual((u.R,u.K,u.J),(old.R,old.K,old.S))
        L=u.at(u.family()['L'],0,0)
        a,p,q,w=u.parts(L)
        kap=(L-a*u.I)/sp.sqrt(p*p+q*q-w*w)
        self.assertEqual(kap.T,kap)
        self.assertEqual(u.cycle(L),u.I)
        self.assertEqual(len(old.group_generated([u.R,kap])),8)

    def test_coordinate_maxwell_holds_while_cross_corner_defect_survives(self):
        S,V=sp.symbols('S V',positive=True)
        U=S**3/V+S**4/V+S**3/V**2
        self.assertEqual(sp.diff(U,S,V),sp.diff(U,V,S))
        hes=sp.hessian(U,(S,V)).subs({S:1,V:1})
        self.assertGreater(hes[0,0],0)
        self.assertGreater(hes.det(),0)
        L=u.at(u.family()['L'],1,1)
        self.assertNotEqual(L[1,0],L[0,1])
        self.assertNotEqual(hes,L)

    def test_symmetrizing_the_source_erases_a_real_return(self):
        L=u.at(u.family()['L'],1,1)
        self.assertNotEqual(u.cycle(L),u.I)
        self.assertEqual(u.cycle((L+L.T)/2),u.I)
        self.assertGreater(u.at(u.family()['fxy'],1,1),0)

    def test_common_energy_unit_cancels_but_source_is_retained(self):
        L=u.at(u.family()['L'],sp.Rational(1,2),sp.Rational(1,3))
        self.assertEqual(u.cycle(L),u.cycle(sp.Rational(7,3)*L))
        _,p,q,w=u.parts(L)
        _,pp,qq,ww=u.parts(7*L)
        self.assertEqual(w*w/(p*p+q*q-w*w),ww*ww/(pp*pp+qq*qq-ww*ww))

    def test_connection_share_is_a_real_hypothesis(self):
        t=sp.symbols('t',positive=True)
        B=sp.diag(t,1/t)
        dB=sp.diff(B,t)
        X=B.inv()*dB
        for share in (sp.Rational(0),sp.Rational(1,4),sp.Rational(1,2),sp.Rational(1)):
            A=share*X
            residual=sp.simplify(dB-A.T*B-B*A)
            self.assertEqual(residual==sp.zeros(2),share==sp.Rational(1,2))

    def test_one_scalar_response_can_be_open_but_connection_flat(self):
        f=u.family()
        # Dropping the third potential term makes all normalized data depend on x alone.
        ell=f['ell'].subs(u.y,0)
        phi_x=f['phix'].subs(u.y,0)
        self.assertGreater(ell.subs(u.x,1),0)
        self.assertNotEqual(phi_x.subs(u.x,1),0)
        self.assertEqual(sp.diff(ell,u.y),0)
        self.assertEqual(sp.diff(phi_x,u.y),0)
        # On this two-state surface a=ell(x) phi_x(x) dx is closed.
        self.assertEqual(sp.diff(ell*phi_x,u.y),0)

    def test_complete_polynomial_certificate_not_sampled_signs(self):
        f=u.family()
        for name,terms,constant in (('H',45,110889),('G',28,945)):
            record=u.polynomial_record(f[name])
            self.assertEqual(record['terms'],terms)
            self.assertEqual(record['constant'],str(constant))
            reconstructed=sum(sp.Rational(t['coefficient'])*u.x**t['powers'][0]*u.y**t['powers'][1]
                              for t in record['coefficients'])
            self.assertTrue(u.zero(reconstructed-f[name]))
            self.assertTrue(all(sp.Rational(t['coefficient'])>0 for t in record['coefficients']))

    def test_oriented_area_gives_closed_protocol_return(self):
        angle=sp.symbols('angle',real=True)
        low,high=sp.Rational(1,4),sp.Rational(9,4)
        width=sp.pi/3
        integral=sp.integrate(high,(angle,0,width))+sp.integrate(low,(angle,width,0))
        self.assertEqual(-integral,-2*sp.pi/3)
        # Reversal reverses the returned angle; identical lost weights cancel.
        self.assertEqual(integral,2*sp.pi/3)
        self.assertEqual(sp.integrate(low,(angle,0,width))+sp.integrate(low,(angle,width,0)),0)

    def test_ordered_return_rotates_without_individual_rotation(self):
        f=u.family()
        b1,b2=u.cycle(u.at(f['L'],1,1)),u.cycle(u.at(f['L'],2,2))
        for b in (b1,b2):
            self.assertEqual(b,b.T)
            self.assertEqual(b.det(),1)
            self.assertGreater(b.trace(),2)
            self.assertEqual(u.polar_tangent(b),0)
        self.assertLess(u.polar_tangent(b2*b1),0)
        self.assertEqual(u.polar_tangent(b1*b2),-u.polar_tangent(b2*b1))

    def test_physical_state_orientation_and_first_visible_order(self):
        lam,S,V=sp.symbols('lambda S V',positive=True)
        jac=sp.Matrix([[sp.diff(lam*S,S),sp.diff(lam*S,V)],
                       [sp.diff(lam/V,S),sp.diff(lam/V,V)]]).det()
        self.assertEqual(jac,-lam**2/V**2)
        lead=sp.Rational(13440,1874161)*u.x
        self.assertEqual(sp.expand(jac*lead.subs(u.x,lam*S)),
                         -sp.Rational(13440,1874161)*lam**3*S/V**2)
        self.assertEqual(u.at(u.family()['fxy'],0,0),0)


if __name__=='__main__':
    unittest.main()
