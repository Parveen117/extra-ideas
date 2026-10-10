"""Independent exact controls for the YC9 analytic gap proof."""
from fractions import Fraction as F
from itertools import product
import unittest
import sympy as sp
import yc9_shared_plane_gap as y


class YC9Tests(unittest.TestCase):
    def test_rational_contraction_and_reserve(self):
        self.assertEqual(y.MAP_COEFFICIENT*y.BETA_MAX, F(1,16))
        self.assertEqual(y.LIPSCHITZ_COEFFICIENT*y.BETA_MAX, F(11,18))
        self.assertEqual(1-y.RELATIVE_COEFFICIENT*y.BETA_MAX, F(8,9))
        self.assertEqual(y.gap_lower(y.BETA_MAX), F(32,3))

    def test_exponential_majorant_by_positive_series_tail(self):
        # For n>=2, (n+1)!/n! >=3; geometric tail begins at n=2.
        for n in range(2,20):
            self.assertGreaterEqual(int(sp.factorial(n)), 2*3**(n-2))
        bound = F(3,2)+F(1,8)/(1-F(1,6))
        self.assertEqual(bound,F(33,20))
        self.assertLess(bound,y.EXP_MAJORANT)

    def test_plane_budget_matches_direct_shared_edge_incidence(self):
        g = y.yc8.geometry((3,4,5))
        couplings = {frozenset((0,1)):F(1,2000), frozenset((1,2)):F(-1,3000),
                     frozenset((0,2)):F(1,2500)}
        face_t = [couplings[frozenset(g['edges'][e][1] for e in face)] for face in g['faces']]
        direct = max(sum(abs(face_t[p]) for p in ps) for ps in g['incident'])
        self.assertEqual(direct,y.plane_budget(F(1,2000),F(-1,3000),F(1,2500)))

    def test_continuous_plane_join_has_no_volume_factor(self):
        t=y.THETA_MAX
        for s in (F(0),F(1,7),F(1,2),F(9,10),F(1)):
            beta=y.plane_budget(t,s*t,s*t)
            self.assertEqual(beta,2*t*(1+s))
            self.assertGreaterEqual(y.gap_lower(beta),F(32,3))

    def test_short_winding_is_not_a_vacuum_mode(self):
        r=y.geometry_record((3,3,3))
        self.assertEqual(r['triangle_count'],27)
        self.assertEqual(r['centre_parities_of_triangles'],[(0,0,1),(0,1,0),(1,0,0)])
        self.assertTrue(r['elementary_faces_centre_even'])
        self.assertEqual(y.geometry_record((4,4,4))['triangle_count'],0)

    def test_resolvent_relative_bound_for_support_sizes(self):
        # Physical free energy simultaneously >=12 and >=3*|I|.
        for size,z in product(range(1,30),(F(-10),F(0),F(2),F(7),F(31,3))):
            floor=max(F(12),3*size)
            lhs=F(size)/(floor-z)
            rhs=F(1,3) if z<0 else F(1,3)/(1-z/12)
            self.assertLessEqual(lhs,rhs)

    def test_creator_and_locality_witnesses(self):
        checks=y.algebra_checks()
        self.assertTrue(all(checks.values()),[name for name,ok in checks.items() if not ok])

    def test_entangled_qutrit_clusters_not_qubit_assumption(self):
        # Each excited single-site space has dimension 2; the clusters
        # are not tensor products of one-site excitation vectors.
        c=sp.zeros(9,1); c[4]=sp.Rational(1,100); c[8]=sp.Rational(2,100)
        d=sp.zeros(9,1); d[4]=sp.Rational(1,200); d[5]=sp.Rational(-1,200); d[7]=sp.Rational(1,200)
        vac=sp.eye(9)[:,0]
        ca=sp.kronecker_product(c*vac.T,sp.eye(3))
        cb=sp.kronecker_product(sp.eye(3),d*vac.T)
        self.assertEqual(ca*cb,sp.zeros(27))
        self.assertEqual(cb*ca,sp.zeros(27))
        V=sp.kronecker_product(sp.eye(3),sp.ones(3)/3,sp.eye(3))
        vector=ca*V*cb*sp.eye(27)[:,0]
        self.assertNotEqual(vector,sp.zeros(27,1))
        patterns=set()
        for row,value in enumerate(vector):
            if value:
                a,b,cindex=row//9,(row//3)%3,row%3
                self.assertNotEqual(a,0)
                self.assertNotEqual(cindex,0)
                patterns.add((a!=0,b!=0,cindex!=0))
        self.assertLessEqual(len(patterns),2)  # only the X={middle} status can vary

    def test_fixed_point_and_excitation_similarity_on_exact_model(self):
        t=sp.Rational(1,160)
        x=sp.Matrix([[0,1],[1,0]])
        H0=sp.diag(0,3,3,6)
        V=-t*sp.kronecker_product(x,x)
        c=(sp.sqrt(9+t*t)-3)/t
        C=sp.zeros(4); C[3,0]=c
        omega=sp.eye(4)[:,0]
        W=(sp.eye(4)-C)*V*(sp.eye(4)+C)
        self.assertEqual(sp.simplify(c+W[3,0]/6),0)
        E=-t*c
        psi=(sp.eye(4)+C)*omega
        self.assertEqual((H0*psi+V*psi-E*psi).applyfunc(sp.simplify),sp.zeros(4,1))
        hat=((sp.eye(4)-C)*(H0+V-E*sp.eye(4))*(sp.eye(4)+C)).applyfunc(sp.simplify)
        self.assertEqual(hat[:,0],sp.zeros(4,1))
        self.assertTrue(3-t-E>sp.Rational(8,3))

    def test_iteration_tail_decreases_and_is_local(self):
        errors=[y.iteration_error(y.BETA_MAX,m) for m in range(8)]
        self.assertTrue(all(a>b>0 for a,b in zip(errors,errors[1:])))
        self.assertEqual(errors[0],F(3,1120))
        self.assertEqual(y.iteration_error(0,0),0)

    def test_sector_free_floor_and_first_return_normalization(self):
        self.assertEqual(4*F(1,2)*(F(1,2)+1),3)
        self.assertEqual(3*4*1*(1+1),24)
        self.assertEqual(4*3,12)
        # Each plaquette cluster has norm |t|/24 and weight four.
        self.assertEqual(4*F(1,24),F(1,6))
        self.assertEqual(F(1,4)/12,F(1,48))

    def test_scope_rejections(self):
        for beta in (F(-1,100),y.BETA_MAX+F(1,10**6)):
            with self.assertRaises(ValueError): y.gap_lower(beta)
        for theta in (F(-1,100),y.THETA_MAX+F(1,10**6)):
            with self.assertRaises(ValueError): y.isotropic_gap(theta)
        for m in (-1,True,F(1,2)):
            with self.assertRaises(ValueError): y.iteration_error(0,m)
        for shape in ((2,3,3),(3,3),(3,3,True)):
            with self.assertRaises(ValueError): y.geometry_record(shape)


if __name__=='__main__':
    unittest.main()
