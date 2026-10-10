"""Independent charge, full-inverse, actual-incidence and scope controls."""
from fractions import Fraction as Q
from itertools import combinations, product
import unittest
import sympy as sp
import yc22_gauge_protected_channels as y


class YC22Tests(unittest.TestCase):
    def test_single_link_normalization_and_shared_star_factor_two(self):
        q = tuple(sp.symbols('q0:4',real=True))
        self.assertEqual(y.sphere_casimir(q[0],q),3*q[0])
        # Degree-one endpoints each have Casimir3, but one link has energy3.
        self.assertEqual(y.charge_pattern_floor((1,1),(1,1)),3)
        self.assertGreater(sum(y.casimir(n) for n in (1,1)),3)
        self.assertEqual(y.charged_floor(1),3)

    def test_global_centre_kernel_and_nonuniform_degrees(self):
        ds = (3,4,4,3)
        for ns in product(range(4),repeat=4):
            if sum(ns)%2:
                with self.assertRaises(ValueError):
                    y.charge_pattern_floor(ns,ds)
            elif any(ns):
                self.assertGreaterEqual(y.charge_pattern_floor(ns,ds),Q(3,4))
            else:
                self.assertEqual(y.charge_pattern_floor(ns,ds),0)

    def test_every_link_derivative_is_one_adjoint_under_gauge_action(self):
        qs,W,G,C,_ = y.normalization_controls()
        for edge in range(3):
            for a in range(3):
                v = y.vector_field(W,qs[edge],a)
                for vertex in range(3):
                    self.assertEqual(C(v,vertex),(8*v if vertex==edge else 0))

    def test_force_identity_uses_the_true_ground_not_a_trial_floor(self):
        w,dw,t = sp.symbols('w dw t',real=True)
        # On the triangle, H0 W=9W and Omega=1+tW is positive for |t|<1.
        V = -9*t*w/(1+t*w)
        h_derivative = 9*t*dw+V*t*dw
        negative_force = -sp.diff(V,w)*dw*(1+t*w)
        self.assertEqual(sp.cancel(h_derivative-negative_force),0)

    def test_full_inverse_preserves_labels_but_leaves_its_source_span(self):
        D = sp.diag(sp.Matrix([[2,sp.Rational(1,3)],[sp.Rational(1,3),4]]),
                    sp.Matrix([[5,sp.Rational(1,2)],[sp.Rational(1,2),7]]))
        C = sp.diag(3,3,8,8)
        self.assertEqual(D*C,C*D)
        s,t = sp.eye(4)[:,0],sp.eye(4)[:,2]
        for z in (0,sp.Rational(1,2),-1):
            inverse = (D-z*sp.eye(4)).inv()
            self.assertEqual(s.T*inverse*t,sp.zeros(1))
            self.assertNotEqual((inverse*s)[1],0)

    def test_same_charge_initial_orthogonality_does_not_survive_mixing(self):
        D = sp.Matrix([[3,1],[1,4]])
        s,t = sp.eye(2)[:,0],sp.eye(2)[:,1]
        self.assertEqual(s.T*t,sp.zeros(1))
        self.assertEqual(D*(3*sp.eye(2)),(3*sp.eye(2))*D)
        self.assertNotEqual(s.T*D.inv()*t,sp.zeros(1))

    def test_charged_floor_does_not_bound_the_physical_excitation(self):
        for epsilon in (sp.Rational(1,10),sp.Rational(1,1000)):
            H = sp.diag(0,epsilon,3)
            C = sp.diag(0,0,6)
            self.assertTrue(y.y21.y20.positive_semidefinite(H-C/2))
            self.assertLess(epsilon,3)
            self.assertEqual(C[:2,:2],sp.zeros(2))

    def test_complete_bridge_flip_is_an_actual_vertex_gauge_transform(self):
        b = y.open_block()
        def endpoint(edge):
            p,a = edge
            q = list(p)
            q[a] += 1
            return tuple(q)
        flipped = {e for e in b['edges'] if (e[0][0]<2)!=(endpoint(e)[0]<2)}
        self.assertEqual(len(flipped),4)
        self.assertTrue(all(e[1]==0 and e[0][0]==1 for e in flipped))
        self.assertTrue(all(sum(e in flipped for e in f)%2==0 for f in b['faces']))

    def test_boundary_force_classes_have_exact_periodic_partners(self):
        records = {tuple(r['edge']):r for r in y.open_block()['records']}
        shape,size = (8,4,4),(4,2,2)
        def shift(p,a):
            q=list(p); q[a]=(q[a]+1)%shape[a]
            return tuple(q)
        def block(p):
            return tuple(p[i]//size[i] for i in range(3))
        found=0
        for p in product(*(range(n) for n in shape)):
            for a,b in combinations(range(3),2):
                face=((p,a),(shift(p,a),b),(shift(p,b),a),(p,b))
                internal=[e for e in face if block(e[0])==block(shift(*e))]
                if len(internal)==2:
                    kinds=[]
                    for x,c in internal:
                        local=(tuple(x[i]%size[i] for i in range(3)),c)
                        rec=records[local]
                        kinds.append((min(rec['degrees']),rec['faces']))
                    self.assertEqual(kinds[0],kinds[1])
                    found+=1
        self.assertEqual(found,48*(shape[0]*shape[1]*shape[2]//16)//2)

    def test_source_residual_caps_retain_the_second_moment(self):
        rows=((Q(9,8),Q(59,28800),Q(91,2000)),
              (Q(2),Q(11,4800),Q(6,125)),
              (Q(9,2),Q(43,14400),Q(11,200)))
        for variance,expected,cap in rows:
            got=y.inverse_square_bound(variance,Q(15,2))
            self.assertEqual(got,expected)
            self.assertGreater(got,Q(1,576))
            self.assertLess(got,cap*cap)

    def test_local_force_bound_uses_the_better_endpoint(self):
        self.assertEqual(y.edge_force_bound((3,4),2),Q(9,16))
        self.assertEqual(y.edge_force_bound((4,4),2),1)
        self.assertEqual(y.edge_force_bound((4,4),3),Q(9,4))
        self.assertEqual(y.edge_force_bound((4,3),0),0)

    def test_resolved_boundary_seed_is_not_a_global_quadrature_replacement(self):
        b=y.open_block()
        costs=[4*r['external_incidence']*r['cap'] for r in b['records']]
        self.assertEqual(sum(costs),Q(228,25))
        self.assertEqual(sum(b['boundary_incidence'].values()),48)
        self.assertGreater(sum(costs)**2,sum(c*c for c in costs))

    def test_physical_bridge_parity_allows_two_parallel_edges_but_not_one(self):
        # Coarse graph with two blocks and two distinct connecting bridges.
        endpoints=((0,1),(0,1))
        allowed=[]
        for parity in product((0,1),repeat=2):
            counts=[sum(parity[e] for e,p in enumerate(endpoints) if v in p) for v in (0,1)]
            if all(c%2==0 for c in counts):
                allowed.append(parity)
        self.assertEqual(allowed,[(0,0),(1,1)])
        self.assertEqual(3*sum(allowed[1]),6)

    def test_physical_gap_uses_delta_without_replacing_g_in_the_budget(self):
        for name in ('cube','tube','force'):
            row=y.join_bounds(name)
            self.assertLess(row['full_factor_floor'],row['physical_reference_floor'])
            self.assertEqual(row['physical_gap']/row['full_gap'],row['physical_reference_floor']/row['full_factor_floor'])
            self.assertGreater(row['ball_reserve'],0)
            self.assertLess(row['contraction'],1)
            self.assertEqual(row['nonlinear_minimum_support'],1)
            self.assertEqual(y.join_bounds(name,0)['physical_gap'],row['physical_reference_floor'])

    def test_invalid_or_unproved_inputs_are_rejected(self):
        for bad in (-1,0,Q(3,2)):
            with self.assertRaises(ValueError): y.charged_floor(bad)
        with self.assertRaises(ValueError): y.inverse_square_bound(1,0)
        with self.assertRaises(ValueError): y.edge_force_bound((3,4),-1)
        with self.assertRaises(ValueError): y.join_bounds('force',Q(1,3000))
        with self.assertRaises(ValueError): y.join_bounds('continuum')


if __name__ == '__main__':
    unittest.main()
