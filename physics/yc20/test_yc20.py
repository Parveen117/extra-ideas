"""YC20: independent generalized-spectrum, residual and frame-loss controls."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import yc20_metric_retained_return as y


class YC20Tests(unittest.TestCase):
    def test_primitive_cut_survives_balanced_central_reading(self):
        # In old coordinates the balanced source is (1,0), not (1,1).
        G = sp.Matrix([[2,1],[1,1]])
        frame = y.metric_frame(3*G,G,1)
        cut = y.compatible_cut(frame)
        source = sp.Matrix([1,0])
        plus,minus = cut['plus']*source,cut['minus']*source
        self.assertEqual((plus.T*G*plus)[0],1)
        self.assertEqual((minus.T*G*minus)[0],1)
        self.assertEqual((source.T*G*cut['cut']*source)[0],0)
        self.assertNotEqual(cut['cut']*source,sp.zeros(2,1))
        self.assertTrue(y.equal(cut['cut']*cut['cut'],sp.eye(2)))
        self.assertTrue(y.equal(cut['cut'].T*G,G*cut['cut']))

    def test_exchange_symmetry_supplies_balance_but_opposite_signs_alone_do_not(self):
        Kcut = sp.diag(1,-1)
        E = sp.Matrix([[0,1],[1,0]])
        source = sp.ones(2,1)
        self.assertEqual(E*Kcut*E,-Kcut)
        self.assertEqual(E*source,source)
        generator = 3*sp.eye(2)+E
        self.assertEqual(generator*E,E*generator)
        self.assertEqual((source.T*(generator.T*Kcut+Kcut*generator)*source)[0],0)
        asymmetric = sp.Matrix([2,1])
        self.assertEqual((asymmetric.T*Kcut*asymmetric)[0],3)
        different_rates = sp.diag(1,2)
        self.assertNotEqual((source.T*(different_rates.T*Kcut+Kcut*different_rates)*source)[0],0)

    def test_cut_weights_are_covariant_under_a_general_observer_change(self):
        K,G = y.complex_fixture()
        cut = y.compatible_cut(y.metric_frame(K,G,1))
        Z = sp.Matrix([[1,sp.I,0],[0,2,sp.Rational(1,3)],[0,0,1]])
        v = sp.Matrix([2,1+sp.I,-1])
        Gnew,vnew = Z.H*G*Z,Z.inv()*v
        Knew = Z.inv()*cut['cut']*Z
        for sign in (1,-1):
            P = (sp.eye(3)+sign*cut['cut'])/2
            Pnew = (sp.eye(3)+sign*Knew)/2
            self.assertTrue(y.equal((P*v).H*G*(P*v),(Pnew*vnew).H*Gnew*(Pnew*vnew)))

    def test_hidden_metric_minimizer_and_complete_square(self):
        K, G = y.complex_fixture()
        frame = y.metric_frame(K, G, 1)
        u = sp.Matrix([2-sp.I])
        v = sp.Matrix([sp.Rational(2,3), 1+sp.I])
        vector = u.col_join(v)
        adjusted = v+frame['X']*u
        expected = u.H*frame['g']*u+adjusted.H*frame['Gh']*adjusted
        self.assertTrue(y.equal(vector.H*G*vector, expected))
        self.assertNotEqual(frame['g'], frame['Gp'])
        optimum = u.col_join(-frame['X']*u)
        self.assertTrue(y.equal(optimum.H*G*optimum, u.H*frame['g']*u))

    def test_reset_changes_the_source_and_preserves_full_pencil(self):
        K, G = y.complex_fixture()
        frame = y.metric_frame(K, G, 1)
        self.assertNotEqual(frame['b'], K[:1,1:])
        z = sp.Rational(1,3)
        T = frame['T']
        expected = (frame['k']-z*frame['g']).row_join(frame['b']).col_join(
            frame['b'].H.row_join(frame['Kh']-z*frame['Gh']))
        self.assertTrue(y.equal(T.H*(K-z*G)*T, expected))
        wrong = K[:1,:1]-z*G[:1,:1]-K[:1,1:]*(K[1:,1:]-z*G[1:,1:]).inv()*K[1:,:1]
        self.assertFalse(y.equal(wrong, y.retained_return(frame,z)['F']))
        A = G.inv()*K
        cut = y.compatible_cut(frame)
        Kcut,P,Qh = cut['cut'],cut['plus'],cut['minus']
        self.assertTrue(y.equal(Qh*A*P,Qh*(A*Kcut-Kcut*A)*P/2))
        self.assertTrue(y.equal((T.inv()*A*T)[1:,:1],frame['Gh'].inv()*frame['b'].H))
        # A scalar dynamics commutes with every cut and has no returned source.
        commuting = y.metric_frame(3*G,G,1)
        self.assertEqual(commuting['b'],sp.zeros(1,2))

    def test_energy_derivatives_match_independent_symbolic_differentiation(self):
        K = sp.Matrix([[7,2],[2,11]])
        G = sp.Matrix([[2,1],[1,3]])
        frame = y.metric_frame(K,G,1)
        z = sp.symbols('z',real=True)
        F, _ = y.schur(K-z*G,1)
        row = y.retained_return(frame,sp.Rational(1,4))
        self.assertTrue(y.equal(-F.diff(z).subs(z,sp.Rational(1,4)),row['slope']))
        self.assertTrue(y.equal(F.diff(z,2).subs(z,sp.Rational(1,4)),row['second']))
        self.assertTrue(y.positive_semidefinite(row['slope']-frame['g']))

    def test_two_energy_step_with_noncommuting_hidden_metric_and_energy(self):
        K,G = y.complex_fixture()
        frame = y.metric_frame(K,G,1)
        self.assertFalse(y.equal(frame['Kh']*frame['Gh'],frame['Gh']*frame['Kh']))
        x,z = sp.Rational(1,7),sp.Rational(2,5)
        left,right = y.retained_return(frame,x),y.retained_return(frame,z)
        self.assertTrue(y.equal(right['F']-left['F'],-(z-x)*right['W'].H*G*left['W']))

    def test_residual_completes_the_return_for_a_complex_inexact_solve(self):
        K,G = y.complex_fixture()
        frame = y.metric_frame(K,G,1)
        z = sp.Rational(1,3)
        Y = sp.Matrix([sp.Rational(1,11)+sp.I/13,-sp.Rational(2,17)])
        enclosure = y.residual_enclosure(frame,z,1,Y)
        D = frame['Kh']-z*frame['Gh']
        full = frame['b']*D.inv()*frame['b'].H
        self.assertTrue(y.equal(full,enclosure['trial']+enclosure['exact_error']))
        self.assertTrue(y.positive_semidefinite(enclosure['error']-enclosure['exact_error']))
        self.assertTrue(y.positive_semidefinite(enclosure['exact']-enclosure['lower']))
        self.assertTrue(y.positive_semidefinite(enclosure['upper']-enclosure['exact']))

    def test_exact_hidden_solve_collapses_both_bounds_to_the_full_return(self):
        K,G = y.complex_fixture()
        frame = y.metric_frame(K,G,1)
        z = sp.Rational(1,2)
        Y = (frame['Kh']-z*frame['Gh']).inv()*frame['b'].H
        enclosure = y.residual_enclosure(frame,z,1,Y)
        self.assertEqual(enclosure['residual'],sp.zeros(2,1))
        self.assertEqual(enclosure['lower'],enclosure['exact'])
        self.assertEqual(enclosure['upper'],enclosure['exact'])

    def test_generalized_inertia_counts_known_physical_levels_including_a_root(self):
        S = sp.Matrix([[1,sp.Rational(1,3)],[0,1]])
        H = sp.diag(2,9)
        K,G = S.T*H*S,S.T*S
        frame = y.metric_frame(K,G,1)
        for z in (0,1,2,3,5,8):
            expected_negative = int(2<z)+int(9<z)
            expected_zero = int(z==2)+int(z==9)
            full = y.inertia(K-z*G)
            reduced = y.inertia(y.retained_return(frame,z)['F'])
            self.assertEqual(full[0],expected_negative)
            self.assertEqual(full[1],expected_zero)
            self.assertEqual(reduced[0],expected_negative)
            self.assertEqual(reduced[1],expected_zero)

    def test_inertia_handles_zero_diagonals_and_complex_two_by_two_pivots(self):
        M = sp.Matrix([[0,sp.I,0],[-sp.I,0,0],[0,0,0]])
        self.assertEqual(y.inertia(M),(1,1,1))
        self.assertEqual(y.inertia(sp.diag(-3,0,4,7)),(1,1,2))
        with self.assertRaises(ValueError):y.inertia(sp.Matrix([[1,1],[0,1]]))

    def test_two_elimination_orders_match_determinants_and_harmonic_lifts(self):
        K,G = y.nested_fixture()
        z = sp.Rational(2,3)
        L = K-z*G
        first,W1 = y.schur(L,2)
        second,W2 = y.schur(first,1)
        total,W = y.schur(L,1)
        self.assertTrue(y.equal(second,total))
        self.assertTrue(y.equal(W1*W2,W))
        self.assertEqual(sp.simplify(L.det()-L[1:,1:].det()*total.det()),0)
        returned = W2.H*(W1.H*G*W1)*W2
        self.assertTrue(y.equal(returned,W.H*G*W))
        metric_only,_ = y.schur(G,1)
        self.assertFalse(y.equal(returned,metric_only))

    def test_metric_schur_is_associative_under_a_different_cut_order(self):
        _,G = y.nested_fixture()
        # Retain original coordinate 2; eliminate other coordinates in a new order.
        order = [2,0,3,1]
        G = G.extract(order,order)
        first,_ = y.schur(G,3)
        final,_ = y.schur(first,1)
        direct,_ = y.schur(G,1)
        self.assertTrue(y.equal(final,direct))

    def test_multichannel_projective_reserve_is_frame_invariant(self):
        K,G = y.nested_fixture()
        Bp = sp.Matrix([[2,sp.Rational(1,3)],[0,1]])
        Bq = sp.Matrix([[1,0],[sp.Rational(-1,5),3]])
        Z = sp.diag(Bp,Bq)
        f = y.metric_frame(K,G,2)
        h = y.metric_frame(Z.T*K*Z,Z.T*G*Z,2)
        # Similar to C C*, without irrational matrix square roots.
        correlation = f['Gp'].inv()*f['N']*f['Gh'].inv()*f['N'].H
        changed = h['Gp'].inv()*h['N']*h['Gh'].inv()*h['N'].H
        self.assertEqual(correlation.charpoly().all_coeffs(),changed.charpoly().all_coeffs())
        self.assertTrue(y.equal(h['g'],Bp.T*f['g']*Bp))

    def test_iteration_reserves_need_a_product_not_just_positive_steps(self):
        row = y.iteration_reserve([Q(1,16),Q(1,32),Q(1,64)])
        self.assertGreaterEqual(row['product'],Q(57,64))
        self.assertEqual(y.iteration_reserve([])['product'],1)
        small = y.iteration_reserve([Q(1,2)]*100)
        self.assertEqual(small['product'],Q(1,2**100))
        self.assertEqual(small['additive_lower'],0)
        for bad in (-1,1,2):
            with self.assertRaises(ValueError):y.iteration_reserve([bad])

    def test_ratio_and_curvature_do_not_supply_spectral_scale(self):
        theta,phi,N,P,F = y.projector_control()
        self.assertTrue(y.equal(N*N,sp.eye(2)))
        self.assertTrue(y.equal(P*P,P))
        for epsilon in (sp.Rational(1,2),sp.Rational(1,1000)):
            H = 2*sp.eye(2)+epsilon*N
            self.assertTrue(y.equal(H*P,(2+epsilon)*P))
            self.assertEqual(H.trace(),4)
            self.assertEqual(sp.trigsimp(H.det()),4-epsilon**2)
        self.assertEqual(sp.trigsimp(F-sp.sin(theta)/2),0)

    def test_metric_ill_conditioning_alone_does_not_close_physical_gap(self):
        for t in (1,2,10):
            S = sp.Matrix([[1,t],[0,1]])
            H = sp.diag(1,4)
            A,G = S.inv()*H*S,S.T*S
            short,_ = y.schur(G,1)
            self.assertEqual(short[0,0]/G[0,0],sp.Rational(1,1+t*t))
            self.assertEqual(A.eigenvals(),{sp.Integer(1):1,sp.Integer(4):1})
            self.assertTrue(y.equal(A.T*G,G*A))
        self.assertEqual(y.inertia((A+A.T)/2)[0],1)

    def test_certificate_rejects_invalid_metric_hidden_floor_and_domains(self):
        K,G = y.complex_fixture()
        frame = y.metric_frame(K,G,1)
        with self.assertRaises(ValueError):y.metric_frame(K,sp.diag(1,1,-1),1)
        with self.assertRaises(ValueError):y.metric_frame(K,G,0)
        with self.assertRaises(ValueError):y.retained_return(frame,sp.I)
        with self.assertRaises(ValueError):y.residual_enclosure(frame,1,1,sp.zeros(2,1))
        with self.assertRaises(ValueError):y.residual_enclosure(frame,0,100,sp.zeros(2,1))
        with self.assertRaises(ValueError):y.residual_enclosure(frame,0,1,sp.zeros(1,1))
        with self.assertRaises(ValueError):y.schur(sp.diag(1,0),1)
        with self.assertRaises(ValueError):y.metric_frame(sp.eye(2),sp.diag(1,1.0),1)

    def test_frozen_predecessor_remains_pinned(self):
        y.y19.y15.verify_predecessor('physics/yc19/YC19_RESULT.json')


if __name__ == '__main__':
    unittest.main()
