"""Full-operator, free-point and compensated-reference controls for YC16."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import yc16_absorbed_cube_return as y


class YC16Tests(unittest.TestCase):
    def test_exact_split_with_rotated_reference_keeps_operator_and_source(self):
        a = sp.Matrix([sp.Rational(3,5), sp.Rational(4,5)])
        b = sp.Matrix([sp.Rational(5,13), sp.Rational(12,13)])
        R = sp.Matrix([[4,1,0,1],[1,5,2,0],[0,2,6,1],[1,0,1,7]])/300
        parts = y.split_return(R,a,b)
        C = parts['connected']
        reconstructed = parts['alpha']*sp.eye(4)+sp.kronecker_product(parts['single_a'],sp.eye(2))+sp.kronecker_product(sp.eye(2),parts['single_b'])+C
        self.assertEqual(reconstructed,R)
        self.assertEqual(parts['va'].T*C*parts['va'],sp.zeros(2))
        self.assertEqual(parts['vb'].T*C*parts['vb'],sp.zeros(2))
        qa,qb=sp.eye(2)-a*a.T,sp.eye(2)-b*b.T
        self.assertEqual(sp.kronecker_product(qa,qb)*C*parts['omega'],C*parts['omega'])

    def test_general_positive_matrix_does_not_supply_heat_product_norm_bound(self):
        R=sp.diag(0,1,1,0)
        parts=y.split_return(R,sp.Matrix([1,0]),sp.Matrix([1,0]))
        self.assertEqual(parts['connected'],sp.diag(0,0,0,-2))
        self.assertGreater(2,1)  # Norm C exceeds norm R: (4) needs its heat proof.

    def test_connected_operator_can_be_nonzero_on_physical_free_faces(self):
        d={9:Q(1,4),17:Q(3,4),3:Q(-1)}
        self.assertEqual(y.square_integral(d,6)/4,Q(3931,599040))
        self.assertGreater(y.square_integral(d,6),0)
        self.assertEqual(y.exponential_integral(d,9)/4,Q(-19,1872))
        self.assertEqual(y.source_bounds(0,0)['connected_source_squared_upper'],0)

    def test_exact_connected_mixed_derivative_uses_the_full_heat_response(self):
        t=sp.symbols('t',nonnegative=True)
        f=sp.exp(-3*t)/84-sp.exp(-9*t)/48+sp.exp(-17*t)/112
        self.assertEqual(sp.integrate(sp.exp(-6*t)*f*f/4,(t,0,sp.oo)),Q(19,86261760))

    def test_charged_duhamel_integrals_have_the_stated_constants(self):
        t=sp.symbols('t',nonnegative=True)
        g=sp.Rational(1,4)
        J=(sp.exp(-g*t)-sp.exp(-3*t))/(3-g)
        self.assertEqual(sp.integrate(sp.exp(-(6+g)*t)*J/2,(t,0,sp.oo)),Q(4,481))
        self.assertEqual(sp.integrate(sp.exp(-6*t)*J*J,(t,0,sp.oo)),Q(4,1443))

    def test_either_free_cube_removes_the_connected_source_not_its_operator(self):
        for theta in (Q(0),Q(1,4),Q(1,2),Q(1)):
            self.assertEqual(y.source_bounds(0,theta)['connected_source_squared_upper'],0)
            self.assertEqual(y.source_bounds(theta,0)['connected_source_squared_upper'],0)
        row=y.source_bounds(1,1)
        self.assertEqual(row['connected_source_squared_upper'],Q(1,8658)**2)
        self.assertEqual(row['connected_operator_norm_upper'],Q(1,24))

    def test_reference_ground_shift_needs_the_offdiagonal_source(self):
        # Exact two-level control of the general complete-complement argument.
        h=sp.diag(0,sp.Rational(1,3))
        B=sp.Matrix([[0,sp.Rational(1,1000)],[sp.Rational(1,1000),sp.Rational(1,2000)]])
        eigenvalues=sorted((h-B).eigenvals(),key=lambda v:float(v))
        ground=eigenvalues[0]
        self.assertLess(ground,0)
        self.assertLessEqual(-ground,sp.Rational(1,1000)**2/(sp.Rational(1,3)-sp.Rational(3,2000)))
        self.assertGreaterEqual(eigenvalues[1]-ground,sp.Rational(1,3)-sp.Rational(3,2000))
        self.assertNotEqual(h-B,h)

    def test_compensation_preserves_the_original_operator_exactly(self):
        h=sp.diag(0,2)
        B=sp.Matrix([[0,sp.Rational(1,100)],[sp.Rational(1,100),sp.Rational(1,200)]])
        interaction=sp.Matrix([[0,sp.Rational(-1,10)],[sp.Rational(-1,10),0]])
        self.assertEqual((h-B)+(interaction+B),h+interaction)
        self.assertNotEqual((h-B)+interaction,h+interaction)

    def test_one_cube_second_order_coefficients_cancel_after_reference_shift(self):
        # Three-state carrier: vacuum, one-cube excitation, bridge source.
        z=sp.symbols('z')
        H0=sp.diag(0,2,6)
        V=sp.Matrix([[0,0,-1],[0,0,-sp.Rational(1,2)],[-1,-sp.Rational(1,2),0]])
        omega=sp.Matrix([1,0,0])
        inverse=sp.diag(0,sp.Rational(1,2),sp.Rational(1,6))
        first=-inverse*V*omega
        second=-inverse*V*first
        A=sp.Matrix([[0,sp.Rational(1,12),0],[sp.Rational(1,12),0,0],[0,0,0]])
        reference_second=inverse*A*omega
        self.assertEqual(first[1],0)
        self.assertEqual(second[1],sp.Rational(1,24))
        self.assertEqual(second[1],reference_second[1])
        self.assertEqual((H0-z*z*A)+(z*V+z*z*A),H0+z*V)

    def test_polynomial_inverse_and_counterterm_envelopes_cover_window(self):
        for i in range(13):
            eta=y.ETA_MAX*i/12
            row=y.reference_bounds(eta)
            eps=row['counterterm_norm_upper']
            s=eta**2/20
            gh=row['new_cube_gap_lower']
            delta=Q(10,27)*eta**2
            shift=eta**4/108
            self.assertGreater(gh,Q(27,100))
            self.assertLessEqual(s*s/gh,shift)
            self.assertLessEqual(2*s/gh,delta)
            self.assertLessEqual(delta/6+(eps+shift)/36,eta**2/13)
            self.assertLessEqual((s+eps*delta)/gh,Q(3,16)*eta**2)

    def test_whole_compensated_window_has_full_cluster_reserve(self):
        for i in range(13):
            eta=y.ETA_MAX*i/12
            row=y.join_bounds(eta)
            self.assertGreater(row['mapping_reserve'],Q(1,40000))
            self.assertLess(row['contraction'],1)
            self.assertLess(row['relative_return'],1)
            self.assertGreater(row['full_gap_lower'],Q(9,40))
            self.assertEqual(row['full_gap_lower'],y.GAP-Q(1024,5)*eta-Q(517,108)*eta**2)
            self.assertEqual(row['nonlinear_minimum_support'],1)
            self.assertLess(row['fixed_point_norm_upper'],y.RADIUS)

    def test_scalar_return_is_not_subtracted_twice(self):
        row=y.join_bounds(y.ETA_MAX)
        self.assertEqual(row['beta'],24*y.ETA_MAX+row['reference']['counterterm_norm_upper'])
        self.assertGreater(y.y15.join_bounds(y.ETA_MAX)['full_gap_lower'],row['full_gap_lower'])

    def test_connected_source_and_expectations_are_recentered_with_the_ground(self):
        ref=y.reference_bounds(y.ETA_MAX)
        self.assertGreater(ref['connected_source_after_reference_upper'],Q(1,8658))
        self.assertLess(ref['connected_source_after_reference_upper'],Q(1,8657))
        delta=ref['new_ground_vector_distance_upper']
        self.assertEqual(ref['recentered_connected_one_cube_cost_per_block_upper'],
                         24*y.ETA_MAX**2*(delta/12+delta**2/6))
        self.assertEqual(ref['recentered_connected_scalar_cost_per_face_upper'],
                         y.ETA_MAX**2*delta**2/6)

    def test_outside_windows_and_wrong_pairings_rejected(self):
        for eta in (-1,y.ETA_MAX+Q(1,10**9)):
            with self.assertRaises(ValueError): y.reference_bounds(eta)
        with self.assertRaises(ValueError): y.source_bounds(2,0)
        with self.assertRaises(ValueError): y.split_return(sp.eye(4),sp.Matrix([2,0]),sp.Matrix([1,0]))
        with self.assertRaises(ValueError): y.split_return(sp.zeros(3),sp.Matrix([1,0]),sp.Matrix([1,0]))

    def test_frozen_source_pins_match(self):
        for path in ('physics/yc12/YC12_RESULT.json','physics/yc15/YC15_RESULT.json'):
            y.y15.verify_predecessor(path)


if __name__=='__main__':
    unittest.main()
