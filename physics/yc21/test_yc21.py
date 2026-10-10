"""Independent spectral, source, incidence and removed-direction controls."""
from fractions import Fraction as Q
from itertools import product
import unittest
import sympy as sp
import yc21_two_cube_block as y


class YC21Tests(unittest.TestCase):
    def test_joint_haar_source_norm_with_signed_unequal_couplings(self):
        # Exact cubature to degree two on each bridge, at unit internal edges.
        axes = [tuple(s*int(i == j) for i in range(4))
                for j in range(4) for s in (-1,1)]
        couplings = (Q(1,2),Q(-1,3),Q(1,5),Q(1,7))
        mean, square = Q(0), Q(0)
        for bridges in product(axes, repeat=4):
            v = sum(t*sum(a*b for a,b in zip(bridges[i],bridges[(i+1)%4]))
                    for i,t in enumerate(couplings))
            mean += v
            square += v*v
        self.assertEqual(mean,0)
        self.assertEqual(square/len(axes)**4,sum(t*t for t in couplings)/4)

    def test_even_floor_must_not_be_applied_to_the_odd_carrier(self):
        self.assertEqual(sum(n*(n+2) for n in (1,0,0,0)),3)
        self.assertEqual(sum(n*(n+2) for n in (1,1,0,0)),6)
        self.assertEqual(sum(n*(n+2) for n in (2,0,0,0)),8)
        # Adding a one-bridge interaction would destroy the selected symmetry.
        original = {frozenset(s) for s in y.SIGNATURES}
        self.assertTrue(all(len(s)%2 == 0 for s in original))
        self.assertFalse(all(len(s)%2 == 0 for s in original|{frozenset((0,))}))

    def test_hidden_source_sign_changes_do_not_remove_return(self):
        D = sp.Matrix([[4,sp.Rational(1,3)],[sp.Rational(1,3),5]])
        R = sp.Matrix([[sp.Rational(1,4),0],[0,sp.Rational(1,5)]])
        full = R.H*D.inv()*R
        self.assertEqual((-R).H*D.inv()*(-R),full)
        self.assertNotEqual(full[0,1],0)
        self.assertNotEqual(D.inv()*R,sp.diag(sp.Rational(1,16),sp.Rational(1,25)))

    def test_gap_budget_covers_unequal_couplings_and_signs(self):
        for t in (Q(0),Q(1,8),Q(1,4),Q(1,2)):
            for ts in ((t,)*4,(t,-t,t/2,0),(t,0,0,0)):
                row = y.block_bounds(ts)
                self.assertTrue(row['positive_gap_certified'])
                self.assertGreaterEqual(row['spectral_reserve'],Q(3572617,443888640))
                self.assertGreaterEqual(row['ground_energy_lower'],-Q(1,16))

    def test_full_finite_spectrum_and_generalized_lift_not_just_schur_root(self):
        H,z,F,W = y.spectral_control()
        zz = sp.Rational(1,5)
        lifted = W.subs(z,zz)
        self.assertTrue(y.y20.equal(lifted.H*(H-zz*sp.eye(4))*lifted,F.subs(z,zz)))
        self.assertEqual(y.y20.inertia(H-zz*sp.eye(4)),(1,0,3))
        self.assertEqual(y.y20.inertia(F.subs(z,zz)),(1,0,1))
        self.assertNotEqual(lifted[3,0],0)

    def test_metric_is_the_complete_physical_norm_and_not_identity(self):
        H,z,F,W = y.spectral_control()
        row = W.subs(z,sp.Rational(1,2))
        M = row.H*row
        u = sp.Matrix([1,sp.I])
        self.assertTrue(y.y20.equal((row*u).H*(row*u),u.H*M*u))
        self.assertTrue(y.y20.positive_semidefinite(M-sp.eye(2)))
        self.assertTrue(y.y20.positive_semidefinite(sp.Rational(50,49)*sp.eye(2)-M))
        self.assertNotEqual(M,sp.eye(2))

    def test_metric_bound_applies_below_zero_as_well_as_at_excitation_energy(self):
        for z in (Q(-10),Q(-1,16),Q(0),Q(1,5),Q(1,2)):
            row = y.block_bounds((Q(1,2),)*4,z)
            self.assertGreaterEqual(row['metric_excess'],0)
            self.assertLessEqual(row['metric_excess'],Q(1,49))

    def test_actual_rectangular_lattice_partition_accounts_for_every_face(self):
        for shape in ((8,4,6),(16,4,4)):
            row = y.rectangular_tiling(shape)
            vertices = shape[0]*shape[1]*shape[2]
            self.assertEqual(sum(row['factor_links'].values()),3*vertices)
            self.assertEqual(len(row['internal'])+len(row['external']),3*vertices)
            blocks = [f for f in row['factor_links'] if f[0]=='block']
            self.assertEqual(len(blocks),vertices//16)
            self.assertEqual(len(row['internal']),16*len(blocks))
            self.assertEqual(sum(n for f,n in row['incidence'].items() if f[0]=='block'),48*len(blocks))
            self.assertIn((3,1),row['bridge_profiles'])

    def test_invalid_tilings_and_unproved_windows_are_rejected(self):
        for shape in ((4,4,4),(8,3,4),(10,4,4),(8,4)):
            with self.assertRaises(ValueError):
                y.rectangular_tiling(shape)
        with self.assertRaises(ValueError):
            y.block_bounds((Q(3,4),)*4)
        with self.assertRaises(ValueError):
            y.block_bounds((0,)*4,6)
        with self.assertRaises(ValueError):
            y.join_bounds(Q(1,10000))

    def test_join_keeps_single_factor_nonlinear_creators(self):
        row = y.join_bounds(y.EPS_MAX)
        self.assertEqual(row['minimum_nonlinear_support'],1)
        self.assertEqual(row['contraction'],Q(81,100))
        self.assertEqual(row['mapping_reserve'],Q(1,3200))
        self.assertEqual(row['gap_lower'],Q(21,125))
        self.assertEqual(y.join_bounds(0)['gap_lower'],Q(1,5))

    def test_curvature_identity_and_nonorthonormal_frame_covariance(self):
        x,z = sp.symbols('x z',real=True)
        W = sp.Matrix([[1,0],[0,1],[x,z]])
        original = y.graph_connection(W,(x,z))
        T = sp.Matrix([[2,x],[0,1]])
        changed = y.graph_connection(W*T,(x,z))
        self.assertTrue(y.y20.equal(changed['metric'],T.H*original['metric']*T))
        self.assertTrue(y.y20.equal(changed['curvature'][0,1],T.inv()*original['curvature'][0,1]*T))
        self.assertTrue(y.y20.equal(changed['metric']*changed['curvature'][0,1],changed['returned'][0,1]))

    def test_real_line_can_have_varying_norm_without_curvature(self):
        x,z = sp.symbols('x z',real=True)
        row = y.graph_connection(sp.Matrix([1,x+z]),(x,z))
        self.assertNotEqual(row['metric'].diff(x),sp.zeros(1))
        self.assertEqual(row['curvature'][0,1],sp.zeros(1))

    def test_distinct_bridge_signatures_force_zero_origin_curvature(self):
        # A diagonal kinetic resolvent on four orthogonal signature carriers.
        # Equal energy is not needed; only their invariance and orthogonality.
        D = sp.diag(6,9,17,26)
        sources = [sp.eye(4)[:,i]/2 for i in range(4)]
        for i,j in product(range(4),repeat=2):
            if i != j:
                self.assertEqual(sources[i].T*D.inv()**2*sources[j],sp.zeros(1))

    def test_complex_projection_retains_both_ordered_paths(self):
        x,z = sp.symbols('x z',real=True)
        W = sp.Matrix([1,x+sp.I*z])
        row = y.graph_connection(W,(x,z))
        self.assertEqual(sp.simplify(row['curvature'][0,1][0]-2*sp.I/(1+x*x+z*z)**2),0)
        self.assertEqual(sp.simplify(row['metric'][0]-(1+x*x+z*z)),0)
        # Ordinary transpose would not give the positive physical pairing.
        self.assertNotEqual(W.T*W,W.H*W)


if __name__ == '__main__':
    unittest.main()
