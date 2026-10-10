"""Focused exact tests of YC10's new support, norm and scaling inputs."""
from fractions import Fraction as F
from itertools import combinations
import unittest
import sympy as sp
import yc10_closed_return_and_scale as y


class YC10Tests(unittest.TestCase):
    def test_both_exact_windows_and_gap_floors(self):
        c3,c4=y.constants(3),y.constants(4)
        self.assertEqual(c3['theta_max'],F(1,64))
        self.assertEqual(c4['theta_max'],F(3,160))
        self.assertEqual(c3['contraction_at_endpoint'],F(39,40))
        self.assertEqual(c4['contraction_at_endpoint'],F(5,6))
        self.assertEqual(y.isotropic_gap(F(1,64),(3,4,5)),F(81,10))
        self.assertEqual(y.isotropic_gap(F(3,160),(4,5,6)),8)

    def test_the_marked_creator_cost_is_not_divided_by_support(self):
        self.assertEqual(y.constants(3)['map_coefficient'],F(39,10))
        self.assertEqual(y.constants(4)['map_coefficient'],F(10,3))
        self.assertEqual(y.constants(3)['lipschitz_coefficient'],F(78,5))
        self.assertEqual(y.constants(4)['lipschitz_coefficient'],F(100,9))

    def test_creator_support_bound_on_overlapping_hyperedges(self):
        for minimum in (3,4):
            weights={I:F(1,100*(1+sum(I))) for size in range(minimum,7)
                     for I in combinations(range(6),size)}
            r=max(sum(len(I)*v for I,v in weights.items() if e in I) for e in range(6))
            X={0,1,2,3}
            touching=sum(v for I,v in weights.items() if X.intersection(I))
            self.assertLessEqual(touching,4*r/minimum)

    def test_orthogonal_projection_gain_with_multidimensional_components(self):
        # Each component is a 2-vector (3a,4a), with norm 5|a|.
        a=[F((-1)**i*(i+1),19) for i in range(16)]
        sum_norm=sum(5*abs(z) for z in a)
        norm_squared=sum(25*z*z for z in a)
        self.assertLessEqual(sum_norm**2,16*norm_squared)
        self.assertEqual(F(16)**2,16*F(16))  # all unit component norms: sharp

    def test_length_three_adjoint_is_even_and_nonzero(self):
        q=sp.symbols('q',real=True)
        chi=4*q*q-1
        self.assertEqual(sp.expand(chi.subs(q,-q)-chi),0)
        self.assertNotEqual(chi,0)
        # SU(2) Haar moments: q^2=1/4, q^4=1/8.
        self.assertEqual(4*F(1,4)-1,0)
        self.assertEqual(16*F(1,8)-8*F(1,4)+1,1)
        self.assertEqual(y.profile_for_shape((3,9,8))['minimum_support'],3)
        self.assertEqual(y.profile_for_shape((4,9,8))['minimum_support'],4)

    def test_exponential_bounds_have_positive_rational_reserves(self):
        self.assertEqual(y.exponential_upper(F(2,3),5),F(404671,207765))
        self.assertLess(y.exponential_upper(F(2,3),5),F(39,20))
        self.assertEqual(y.exponential_upper(F(1,2),1),F(33,20))
        self.assertLess(y.exponential_upper(F(1,2),1),F(5,3))

    def test_full_tower_remainder_contracts_without_claiming_finite_closure(self):
        for s in (3,4):
            beta=y.constants(s)['beta_max']
            q=y.constants(s)['contraction_at_endpoint']
            errors=[y.iteration_error(beta,s,m) for m in range(4)]
            self.assertTrue(all(v>0 for v in errors))
            for a,b in zip(errors,errors[1:]): self.assertEqual(b,q*a)
        self.assertEqual(y.iteration_error(0,3,0),0)

    def test_new_bound_improves_yc9_on_its_original_window(self):
        theta=y.yc9.THETA_MAX
        old=y.yc9.isotropic_gap(theta)
        self.assertGreater(y.isotropic_gap(theta,(3,3,3)),old)
        self.assertGreater(y.isotropic_gap(theta,(4,4,4)),old)
        self.assertEqual(y.constants(3)['theta_max']/theta,10)
        self.assertEqual(y.constants(4)['theta_max']/theta,12)

    def test_action_and_spatial_energy_have_different_scale_weights(self):
        for eps in (F(1,2),F(1,11),F(2,3)):
            four=y.curvature_scaling(4,eps)
            three=y.curvature_scaling(3,eps)
            self.assertEqual(four['curvature_density']*four['measure'],1)
            self.assertEqual(three['curvature_density']*three['measure'],1/eps)
            self.assertGreater(four['curvature_density'],1)

    def test_symbolic_nonabelian_and_native_scale_controls(self):
        result=y.symbolic_checks()
        self.assertTrue(all(result.values()),[k for k,v in result.items() if not v])

    def test_large_scalar_energy_does_not_control_a_nontrivial_gap(self):
        eps,z=sp.symbols('eps z',positive=True)
        H=sp.Matrix([[eps**-4,1],[1,eps**-4+2]])
        centred=H-eps**-4*sp.eye(2)
        self.assertEqual(sp.expand((centred-z*sp.eye(2)).det()),z*z-2*z-1)
        self.assertEqual(set(centred.eigenvals()),{1-sp.sqrt(2),1+sp.sqrt(2)})

    def test_scope_and_geometry_guards(self):
        with self.assertRaises(ValueError): y.isotropic_gap(F(3,160),(3,4,4))
        with self.assertRaises(ValueError): y.isotropic_gap(F(1,10),(4,4,4))
        with self.assertRaises(ValueError): y.isotropic_gap(F(1,100),(2,4,4))
        with self.assertRaises(ValueError): y.constants(2)
        with self.assertRaises(ValueError): y.curvature_scaling(4,0)
        with self.assertRaises(ValueError): y.exponential_upper(4,1)
        with self.assertRaises(ValueError): y.iteration_error(0,4,-1)


if __name__=='__main__':
    unittest.main()
