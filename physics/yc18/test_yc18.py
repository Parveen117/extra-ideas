"""Topology, full-complement bounds and independent series controls for YC18."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import yc18_connected_fourth_order as y


class YC18Tests(unittest.TestCase):
    def test_complete_pair_xor_search_includes_four_bridge_faces(self):
        for shape,expected in (((4,4,4),24),((4,6,8),72),((6,6,6),81)):
            row=y.lattice_audit(shape)
            self.assertEqual(row['third_order_distinct_returns'],0)
            self.assertEqual(row['fourth_order_distinct_returns'],expected)
            self.assertEqual(row['max_shared_bridges'],1)
            self.assertTrue(row['all_fourth_distinct_are_exactly_tubes'])
            self.assertTrue(row['all_fourth_distinct_have_only_two_bridge_faces'])

    def test_repetitions_leave_only_the_declared_odd_multiplicity_sets(self):
        # Two distinct signatures cannot cancel; the p^3 q and p^2 q r
        # partitions reduce to that same obstruction.
        masks=(3,6,12,9,0b11110000)
        for i,a in enumerate(masks):
            self.assertNotEqual(a^a^a,0)
            self.assertEqual(a^a^a^a,0)
            for j,b in enumerate(masks):
                if i==j:continue
                self.assertNotEqual(a^b,0)
                self.assertNotEqual(a^a^a^b,0)
                self.assertEqual(a^a^b^b,0)
        self.assertEqual(masks[0]^masks[1]^masks[2]^masks[3],0)

    def test_even_nonconstant_bridge_floor_is_eight(self):
        self.assertEqual(min(l*(l+2) for l in range(1,8) if l%2),3)
        self.assertEqual(min(l*(l+2) for l in range(1,8) if not l%2),8)
        self.assertEqual(y.self_four_bound(2),Q(1,4*6*8*6))
        self.assertEqual(y.self_four_bound(4),Q(1,4*12*8*12))

    def test_six_order_pair_sum_and_bridge_disjoint_zero_orders(self):
        for m,n in ((2,2),(2,4),(4,2),(4,4)):
            for shared in (0,1):
                value,count=y.paired_bound_by_words(m,n,shared)
                self.assertEqual(value,y.paired_bound(m,n,shared))
                self.assertEqual(value,y.paired_bound(n,m,shared))
                self.assertEqual(count,6 if shared else 4)
        self.assertLess(y.paired_bound(2,2,0),y.paired_bound(2,2,1))

    def test_spectator_spectral_shift_is_exact_and_not_identity_extension(self):
        ds=[sp.diag(6,9),sp.diag(6,11),sp.diag(6,13)]
        X=sp.Matrix([[0,1],[1,0]])
        Z=sp.diag(1,-1)
        S=sp.diag(0,2,5)
        inv=[(sp.kronecker_product(d,sp.eye(3))+sp.kronecker_product(sp.eye(2),S)).inv() for d in ds]
        global_word=inv[2]*sp.kronecker_product(X,sp.eye(3))*inv[1]*sp.kronecker_product(Z,sp.eye(3))*inv[0]/4
        def local(e):
            r=[(d+e*sp.eye(2)).inv() for d in ds]
            return r[2]*X*r[1]*Z*r[0]/4
        shifted=sp.zeros(6)
        for i,e in enumerate((0,2,5)):
            p=sp.zeros(3);p[i,i]=1
            shifted+=sp.kronecker_product(local(e),p)
        self.assertEqual(global_word,shifted)
        self.assertNotEqual(global_word,sp.kronecker_product(local(0),sp.eye(3)))
        embed=sp.kronecker_product(sp.eye(2),sp.Matrix([1,0,0]))
        self.assertEqual(global_word*embed,embed*local(0))
        # One ordered adjacent word has derivative cap (1/864)*(3/6).
        for e in (2,5):
            difference=local(e)-local(0)
            positive=(sp.Rational(e,1728))**2*sp.eye(2)-difference.T*difference
            self.assertGreaterEqual(positive[0,0],0)
            self.assertGreaterEqual(positive[1,1],0)
            self.assertGreaterEqual(positive.det(),0)

    def test_tube_spectator_lipschitz_constant_counts_all_three_positions(self):
        adjacent=16*Q(1,4*6**3)*Q(3,6)
        opposite=8*Q(1,4*6*12*6)*(Q(2,6)+Q(1,12))
        self.assertEqual(adjacent+opposite,Q(29,2592))

    def test_raw_disconnected_schur_return_cancels_at_the_actual_energy(self):
        row=y.disconnected_control()
        self.assertEqual(row['raw_fourth_cross'],Q(1,13824))
        self.assertEqual(row['energy_argument_cross'],row['raw_fourth_cross'])
        self.assertEqual(row['net_energy_cross'],0)
        self.assertTrue(row['all_mixed_creator_components_zero'])

    def test_general_recursion_matches_independent_rayleigh_logarithm(self):
        energies,V=y.chain_fixture()
        psi,energy=y.rayleigh_coefficients(energies,V)
        direct=y.commutator_coefficients(energies,V)
        self.assertEqual(direct,y.creation_log(psi))
        self.assertTrue(any(direct[n]!=psi[n+1] for n in range(4)))
        H,_,_=y.qubit_reference(energies)
        for n in range(1,5):
            rhs=V*psi[n-1]
            for k in range(1,n+1):rhs+=energy[k]*psi[n-k]
            self.assertEqual(sp.simplify(H*psi[n]-rhs),sp.zeros(8,1))

    def test_creators_commute_but_overlapping_supports_multiply_to_zero(self):
        vectors=[]
        for mask in (1,2,3):
            vector=sp.zeros(4,1);vector[mask]=1;vectors.append(vector)
        a,b,ab=[y.creator(v) for v in vectors]
        self.assertEqual(a*b,b*a)
        self.assertEqual(a*b,ab)
        self.assertEqual(a*ab,sp.zeros(4))
        self.assertEqual(ab*b,sp.zeros(4))

    def test_local_coefficients_ignore_disconnected_spectator_gap(self):
        X=sp.Matrix([[0,1],[1,0]])
        first=y.commutator_coefficients((2,3),y.factor_matrix(X,0,2))
        changed=y.commutator_coefficients((2,101),y.factor_matrix(X,0,2))
        self.assertEqual(first,changed)
        for c in first:self.assertEqual(list(c[2:]),[0,0])

    def test_cauchy_tail_is_uniform_profile_bound_not_endpoint_certificate(self):
        self.assertEqual(y.tail_bound(0),0)
        self.assertEqual(y.tail_bound(y.RHO/2),Q(1,2048))
        self.assertEqual(y.tail_bound(y.RHO/4),Q(1,98304))
        self.assertGreater(y.tail_bound(3*y.RHO/4),y.tail_bound(y.RHO/2))
        for x in (Q(1,4),Q(1,2),Q(3,4)):
            finite=sum(y.BALL*x**n for n in range(5,30))
            self.assertLess(finite,y.tail_bound(x*y.RHO))
        row=y.y15.join_bounds(y.RHO)
        self.assertLess(row['contraction'],1)
        self.assertGreater(row['mapping_reserve'],0)

    def test_free_tube_kinetic_source_needs_the_final_energy_inverse(self):
        returned=Q(y.y17.tube_bounds()['free_vacuum_face_pair_coefficient'])
        creator=returned/24
        self.assertEqual(creator,Q(11,3981312))
        self.assertEqual(creator/4,Q(11,15925248))
        self.assertEqual(2*creator/4,Q(11,7962624))
        self.assertNotEqual(creator,returned)

    def test_outside_declared_ranges_are_rejected(self):
        for eta in (-1,y.RHO,2*y.RHO):
            with self.assertRaises(ValueError):y.tail_bound(eta)
        with self.assertRaises(ValueError):y.paired_bound(2,4,2)
        with self.assertRaises(ValueError):y.self_four_bound(3)
        with self.assertRaises(ValueError):y.creator(sp.Matrix([1,0]))
        with self.assertRaises(ValueError):y.qubit_reference((0,3))

    def test_frozen_inputs_match(self):
        y.y15.verify_predecessor('physics/yc17/YC17_RESULT.json')
        y.y15.verify_predecessor('physics/yc15/YC15_RESULT.json')


if __name__=='__main__':unittest.main()
