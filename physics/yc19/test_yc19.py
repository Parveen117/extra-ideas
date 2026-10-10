"""Independent local-operator, metric, range and resolvent controls for YC19."""
from fractions import Fraction as Q
from itertools import product
import unittest
import sympy as sp
import yc19_local_excitation_operator as y


class YC19Tests(unittest.TestCase):
    def test_creation_lift_is_local_and_not_a_global_rank_one(self):
        vector=sp.Matrix([sp.Rational(1,3),sp.Rational(1,5),0,sp.Rational(1,7)])
        lifted=y.vacuum_lift(vector,(2,2))
        self.assertEqual(lifted[:,0],vector)
        self.assertNotEqual(lifted,vector*sp.eye(4)[0,:])
        A=sp.Matrix([[1,2,0,1],[3,4,1,0],[2,0,5,1],[1,2,3,6]])/20
        self.assertEqual(y.normal_order(A,(2,2))[:,0],sp.zeros(4,1))
        extended=sp.kronecker_product(A,sp.eye(3))
        self.assertEqual(y.normal_order(extended,(2,2,3)),
                         sp.kronecker_product(y.normal_order(A,(2,2)),sp.eye(3)))

    def test_all_excited_vectors_keep_the_complete_commutator_action(self):
        A=sp.Matrix([[1,2,0,1],[3,4,1,0],[2,0,5,1],[1,2,3,6]])/20
        u=sp.Matrix([0,sp.Rational(2,7),sp.Rational(-3,11),sp.Rational(5,13)])
        C=y.y18.creator(u)
        self.assertEqual(y.normal_order(A,(2,2))*u,(A*C-C*A)*sp.eye(4)[:,0])
        scalar_only=A-A[0,0]*sp.eye(4)
        self.assertNotEqual(scalar_only*u,y.normal_order(A,(2,2))*u)

    def test_projection_count_depends_on_original_face_not_excited_dimension(self):
        # Three qutrit factors, root X={middle}; two entangled creators.
        c=sp.zeros(9,1);c[4]=sp.Rational(1,10);c[8]=sp.Rational(2,10)
        d=sp.zeros(9,1);d[4]=sp.Rational(1,20);d[5]=sp.Rational(-1,20);d[7]=sp.Rational(1,20)
        vac=sp.eye(9)[:,0]
        left=sp.kronecker_product(c*vac.T,sp.eye(3))
        right=sp.kronecker_product(sp.eye(3),d*vac.T)
        V=sp.kronecker_product(sp.eye(3),sp.ones(3)/3,sp.eye(3))
        A=left*V*right
        source=A[:,0]
        patterns=set()
        for row,digits in enumerate(product(range(3),repeat=3)):
            if source[row]!=0:
                self.assertTrue(digits[0]!=0 and digits[2]!=0)
                patterns.add(tuple(v!=0 for v in digits))
        self.assertTrue(patterns)
        self.assertLessEqual(len(patterns),2)
        self.assertEqual(y.normal_order(A,(3,3,3))[:,0],sp.zeros(27,1))

    def test_complete_excitation_polynomial_matches_independent_similarity(self):
        energies,V=y.y18.chain_fixture()
        direct=y.excitation_coefficients(energies,V)
        independent=y.similarity_coefficients(energies,V)
        self.assertEqual(direct,independent)
        self.assertTrue(any(A!=A.T for A in direct))
        for A in direct:self.assertEqual(A[:,0],sp.zeros(8,1))

    def test_exact_operator_normal_order_includes_scalar_and_ground_response(self):
        row=y.exact_metric_control()
        self.assertEqual(row['dressed'],row['H0']+y.normal_order(row['W'],(2,2)))
        self.assertEqual(row['dressed'][:,0],sp.zeros(4,1))
        self.assertNotEqual(row['dressed'],row['dressed'].T)
        self.assertEqual(row['H'],row['H'].T)
        self.assertNotEqual(row['S'].T*row['S'],sp.eye(4))

    def test_full_and_quotient_metrics_both_travel_with_similarity(self):
        row=y.exact_metric_control()
        H,M,A,G=row['dressed'],row['metric'],row['A'],row['quotient']
        self.assertEqual(H.T*M,M*H)
        self.assertEqual(A.T*G,G*A)
        self.assertNotEqual(G,M[1:,1:])  # Projecting the metric is insufficient.
        for n in (1,2,3):self.assertGreater(G[:n,:n].det(),0)
        upper=H[0:1,1:]
        self.assertEqual(upper,-M[0:1,1:]*A/M[0,0])

    def test_local_root_budget_and_relative_budget_have_different_constants(self):
        self.assertEqual(5*24*(1+2*y.BALL)*Q(16,15),130)
        self.assertEqual(y.full_local_bound(y.RHO),Q(13,432))
        self.assertEqual(y.RELATIVE,y.y15.join_bounds(y.RHO)['relative_return'])
        self.assertNotEqual(y.LOCAL,y.RELATIVE)

    def test_tail_bound_controls_full_operators_in_both_norms(self):
        half=y.truncation_bounds(y.RHO/2)
        quarter=y.truncation_bounds(y.RHO/4)
        self.assertEqual(half['local_tail_upper'],Q(13,6912))
        self.assertEqual(half['relative_tail_upper'],Q(346112,31987775))
        self.assertEqual(quarter['local_tail_upper'],Q(13,331776))
        self.assertGreater(half['local_tail_upper'],quarter['local_tail_upper'])
        self.assertEqual(y.truncation_bounds(0)['relative_tail_upper'],0)

    def test_size_and_exponential_range_bounds_sum_every_order(self):
        x=Q(1,4)
        row=y.truncation_bounds(x*y.RHO,range_base=2)
        self.assertEqual(row['exponential_range_upper'],y.LOCAL)
        self.assertEqual(row['size_weighted_upper'],Q(13,243))
        self.assertLess(sum(4*n*y.LOCAL*x**n for n in range(1,24)),row['size_weighted_upper'])
        self.assertLess(sum(y.LOCAL*(2*x)**n for n in range(1,24)),row['exponential_range_upper'])

    def test_resolvent_certificate_uses_relative_error_and_can_decline(self):
        row=y.resolvent_transfer(y.RHO/2,y.GAP/2)
        self.assertTrue(row['certified'])
        self.assertEqual(row['error_product'],Q(692224,21604415))
        self.assertGreater(row['exact_relative_inverse_upper'],row['approximate_relative_inverse_upper'])
        # The approximate resolvent bound may exist while the full-error
        # transfer does not certify it; do not silently report a spectrum.
        b4=y.truncation_bounds(y.RHO/2)['truncated_relative_upper']
        delta=y.truncation_bounds(y.RHO/2)['relative_tail_upper']
        rejected=y.resolvent_transfer(y.RHO/2,y.GAP*(1-b4-delta/2))
        self.assertFalse(rejected['certified'])
        self.assertIsNone(rejected['exact_relative_inverse_upper'])

    def test_regrouping_preserves_terms_and_can_increase_incidence(self):
        terms=[({0,1},Q(1,5)),({1,2},Q(1,7)),({3,4},Q(1,11)),({4,5},Q(1,13))]
        partition={i:i//2 for i in range(6)}
        fine=y.incidence_budget(terms)
        coarse=y.incidence_budget(terms,partition)
        self.assertLessEqual(coarse,2*fine)
        separate=[({0},1),({1},1)]
        self.assertEqual(y.incidence_budget(separate,{0:0,1:0}),2*y.incidence_budget(separate))
        bounds=y.truncation_bounds(y.RHO/4,block_size=8)
        self.assertEqual(bounds['regrouped_tail_upper'],8*bounds['local_tail_upper'])

    def test_domain_and_scope_rejections(self):
        self.assertGreater(y.full_local_bound(y.RHO),0)
        with self.assertRaises(ValueError):y.truncation_bounds(y.RHO)
        with self.assertRaises(ValueError):y.truncation_bounds(y.RHO/2,range_base=2)
        with self.assertRaises(ValueError):y.truncation_bounds(y.RHO/4,block_size=0)
        with self.assertRaises(ValueError):y.full_local_bound(-1)
        with self.assertRaises(ValueError):y.vacuum_lift(sp.zeros(3,1),(2,2))
        with self.assertRaises(ValueError):y.nilpotent_exp(sp.eye(2))
        with self.assertRaises(ValueError):y.resolvent_transfer(y.RHO/2,y.GAP)

    def test_frozen_predecessors_are_pinned(self):
        y.y15.verify_predecessor('physics/yc18/YC18_RESULT.json')
        y.y15.verify_predecessor('physics/yc15/YC15_RESULT.json')


if __name__=='__main__':unittest.main()
