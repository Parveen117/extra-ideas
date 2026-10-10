"""Independent tensor, topology, parity and scale controls for YC17."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations,product
import unittest
import sympy as sp
import yc17_four_face_bridge_return as y


class YC17Tests(unittest.TestCase):
    def test_direct_haar_quadrature_at_noncommuting_rational_links(self):
        a=[(Q(3,5),Q(4,5),0,0),(Q(5,13),0,Q(12,13),0),y.AXES[3],y.AXES[0]]
        b=[y.AXES[0],y.AXES[1],(Q(3,5),0,0,Q(4,5)),y.AXES[2]]
        self.assertNotEqual(y.mul(a[0],a[1]),y.mul(a[1],a[0]))
        total=Q(0)
        # Each bridge has homogeneous degree two, so its four unit-axis
        # second-moment quadrature is exact even without the negative axes.
        for us in product(y.AXES,repeat=4):
            term=Q(1)
            for i in range(4):term*=y.y15.wilson(a[i],b[i],us[(i+1)%4],us[i])
            total+=term/256
        self.assertEqual(total,Q(-3,325))
        self.assertEqual(y.static_tube(a,b),total)
        self.assertEqual(total,y.quaternion_chain(a)[0]*y.quaternion_chain(b)[0]/64)

    def test_tensor_identity_is_complete_not_only_identity_link_control(self):
        self.assertEqual(y.tensor_checks(),(True,True))
        self.assertEqual(y.static_tube([y.AXES[0]]*4,[y.AXES[0]]*4),Q(1,64))
        # Reversing a single physical link can change the sign; do not use
        # absolute traces or discard orientation before gluing.
        a=[y.AXES[1],y.AXES[1],y.AXES[0],y.AXES[0]]
        self.assertEqual(y.static_tube(a,[y.AXES[0]]*4),Q(-1,64))
        a[0]=y.conj(a[0])
        self.assertEqual(y.static_tube(a,[y.AXES[0]]*4),Q(1,64))

    def test_every_proper_subset_has_an_exposed_bridge(self):
        for size in (1,2,3):
            for faces in combinations(range(4),size):
                incidence=Counter(v for f in faces for v in (f,(f+1)%4))
                self.assertTrue(any(count%2 for count in incidence.values()))
        # The larger odd bridge floor is needed only for opposite faces.
        for first,second in combinations(range(4),2):
            mask=y.SIGNATURES[first]^y.SIGNATURES[second]
            self.assertEqual(mask.bit_count(),4 if (first-second)%2==0 else 2)

    def test_small_periodic_lattice_keeps_parallel_interfaces_distinct(self):
        tubes,faces,_=y.tube_geometry((4,4,4))
        pairs=Counter(frozenset(t['cubes']) for t in tubes)
        self.assertEqual(len(tubes),24)
        self.assertEqual(set(pairs.values()),{2})
        self.assertEqual(len({frozenset(t['bridges']) for t in tubes}),24)
        self.assertEqual(len(faces),96)

    def test_actual_bridge_graph_is_union_of_four_cycles_on_asymmetric_torus(self):
        tubes,faces,_=y.tube_geometry((4,6,8))
        counts=Counter(f for t in tubes for f in t['faces'])
        self.assertEqual(set(counts),faces)
        self.assertEqual(set(counts.values()),{1})
        for t in tubes:
            degree=Counter(e for face in t['faces'] for e in face if e in t['bridges'])
            self.assertEqual(set(degree.values()),{2})
            self.assertEqual(len(degree),4)

    def test_full_bound_agrees_with_independent_adjacent_opposite_count(self):
        for z in (Q(-2),Q(0),Q(1,2),Q(2),Q(299,100)):
            row=y.tube_bounds(z)
            self.assertEqual(row['full_operator_norm_upper'],4/(6-z)**3+2/((6-z)**2*(12-z)))
            exact=Q(1,64)*(16/((12-z)*(18-z)*(24-z))+8/((12-z)*(24-z)**2))
            self.assertEqual(row['free_vacuum_face_pair_coefficient'],exact)
            self.assertLess(row['free_vacuum_source_norm'],row['full_operator_norm_upper'])

    def test_completed_bridge_adjoint_has_no_later_insertion_to_return_it(self):
        # In every permutation a bridge occurs twice. After its second
        # occurrence, its constant projection commutes past all later faces.
        from itertools import permutations
        for order in permutations(range(4)):
            for bridge in range(4):
                hits=[step for step,face in enumerate(order) if bridge in (face,(face+1)%4)]
                self.assertEqual(len(hits),2)
                self.assertTrue(all(bridge not in (face,(face+1)%4) for face in order[hits[-1]+1:]))
        self.assertEqual([r['free_surviving_energy'] for r in y.prefix_data((0,1,2,3))],[12,18,24])
        self.assertEqual([r['free_surviving_energy'] for r in y.prefix_data((0,2,1,3))],[12,24,24])

    def test_static_moment_and_kinetic_coefficient_are_different_targets(self):
        self.assertEqual(Q(1,64)*Q(1,4)**2,Q(1,1024))
        row=y.tube_bounds()
        self.assertEqual(row['free_vacuum_face_pair_coefficient'],Q(11,165888))
        self.assertEqual(row['free_vacuum_source_norm'],Q(11,663552))
        self.assertEqual(row['incident_selected_coefficient_cost'],Q(5,36))
        self.assertEqual(y.y16.source_bounds(0,0)['connected_source_squared_upper'],0)

    def test_connection_commutator_without_derivative_is_not_general_curvature(self):
        self.assertTrue(all(y.curvature_controls().values()))
        # A constant off-diagonal connection may have no curvature at all.
        R=sp.Matrix([[0,-1],[1,0]])
        self.assertEqual(R*(2*R)-(2*R)*R,sp.zeros(2))

    def test_invalid_carrier_and_resolvent_domain_rejected(self):
        for shape in ((3,4,4),(4,4),(4,4,5)):
            with self.assertRaises(ValueError):y.tube_geometry(shape)
        for z in (3,4):
            with self.assertRaises(ValueError):y.tube_bounds(z)
        with self.assertRaises(ValueError):y.prefix_data((0,0,2,3))
        with self.assertRaises(ValueError):y.static_tube([y.AXES[0]]*3,[y.AXES[0]]*4)

    def test_predecessor_is_byte_pinned(self):
        y.y15.verify_predecessor('physics/yc16/YC16_RESULT.json')


if __name__=='__main__':unittest.main()
