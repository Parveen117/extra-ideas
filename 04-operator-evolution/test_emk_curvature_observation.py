"""Meaningful exact R8 checks; general proofs are separate from finite PASS."""

import importlib.util
import itertools
import os
from pathlib import Path
import unittest

import emk_curvature_observation as c

Q, I, K, R, RK = c.Q, c.identity(2), c.K, c.R, c.RK
Z = c.scale(I, 0)


class CurvatureObservationTests(unittest.TestCase):
    def test_native_nonzero_curvature_and_trace_blindness(self):
        f = c.commutator(K, R)
        self.assertEqual(f, c.matrix([[2, 0], [0, -2]]))
        self.assertEqual(c.trace(f), 0)
        self.assertFalse(c.is_zero(f))
        for a, b in itertools.product((I, K, R, RK), repeat=2):
            self.assertEqual(c.trace(c.commutator(a, b)), 0)

    def test_compression_identity_rational_grid_and_oblique_projectors(self):
        p0 = c.matrix([[1, 0], [0, 0]])
        p_plus, p_minus = c.parity_projectors(K)
        oblique = c.matrix([[1, Q(2, 3)], [0, 0]])
        operators = (I, K, R, RK, c.add(K, c.scale(R, Q(3, 5))))
        for p, a, b in itertools.product((p0, p_plus, p_minus, oblique), operators, operators):
            with self.subTest(p=p, a=a, b=b):
                self.assertEqual(c.compression_report(a, b, p)['identity_residual'], Z)

    def test_full_curved_reduced_flat_with_excursion_cancellation(self):
        p = c.matrix([[1, 0], [0, 0]])
        report = c.compression_report(K, R, p)
        self.assertFalse(report['full_flat'])
        self.assertTrue(report['reduced_flat'])
        self.assertEqual(report['compressed_full'], c.scale(p, 2))
        self.assertEqual(report['excursions'], c.scale(p, 2))

    def test_full_flat_reduced_curved_three_mode_witness(self):
        report = c.exact_witnesses()['full_flat_reduced_curved']
        expected = c.scale(c.matrix([[0, -1, 1], [1, 0, -1], [-1, 1, 0]]), Q(1, 9))
        self.assertTrue(report['full_flat'])
        self.assertFalse(report['reduced_flat'])
        self.assertEqual(report['reduced'], expected)
        self.assertEqual(report['excursions'], c.scale(expected, -1))
        self.assertTrue(c.is_zero(report['identity_residual']))

    def test_compression_covariance_with_observer_transport(self):
        p = c.matrix([[1, Q(2, 3)], [0, 0]])
        s = c.matrix([[1, 2], [0, 1]])
        original = c.compression_report(K, R, p)
        changed = c.compression_report(c.conjugate(K, s), c.conjugate(R, s), c.conjugate(p, s))
        for key in ('full', 'compressed_full', 'reduced', 'excursions'):
            self.assertEqual(changed[key], c.conjugate(original[key], s))

    def test_parity_same_sector_flatness_hides_cross_sector_coupling(self):
        p, q = c.parity_projectors(K)
        f = c.commutator(K, R)
        self.assertEqual(c.mul(c.mul(p, f), p), Z)
        self.assertEqual(c.mul(c.mul(q, f), q), Z)
        self.assertNotEqual(c.mul(c.mul(p, f), q), Z)
        self.assertNotEqual(c.mul(c.mul(q, f), p), Z)

    def test_two_readouts_recover_every_declared_odd_target(self):
        for x, y in itertools.product((Q(-2), Q(-1, 3), Q(0), Q(3, 5), Q(2)), repeat=2):
            f = c.add(c.scale(R, x), c.scale(RK, y))
            self.assertEqual(c.odd_readouts(f), (x-y, -x-y))
            self.assertEqual(c.recover_odd(*c.odd_readouts(f)), f)

    def test_one_channel_has_nonzero_blind_target(self):
        f = c.add(R, RK)
        self.assertEqual(c.odd_readouts(f)[0], 0)
        self.assertNotEqual(c.odd_readouts(f)[1], 0)
        self.assertFalse(c.is_zero(f))

    def test_active_cut_sign_and_passive_reference_covariance(self):
        f = c.add(c.scale(R, Q(2, 3)), c.scale(RK, Q(-4, 5)))
        original = c.odd_readouts(f)
        self.assertEqual(c.odd_readouts(c.mul(c.mul(K, f), K)), tuple(-x for x in original))
        s = c.matrix([[2, 1], [1, 1]])
        changed = c.conjugate(f, s)
        for ell, w, expected in ((c.ELL_PLUS, c.W_MINUS, original[0]),
                                 (c.ELL_MINUS, c.W_PLUS, original[1])):
            left, right = c.transport_probe(ell, w, s)
            self.assertEqual(c.pair(left, changed, right), expected)

    def test_curvature_does_not_identify_generator_factors(self):
        def native(b, cc, d):
            return c.commutator(c.scale(K, b), c.add(c.scale(R, cc), c.scale(RK, d)))
        self.assertEqual(native(1, 3, 2), native(2, Q(3, 2), 1))

    def test_static_two_channels_not_future_complete(self):
        # Same current odd readout, distinguishable after lawful native left action.
        self.assertEqual(c.pair(c.ELL_PLUS, I, c.W_MINUS), 0)
        self.assertEqual(c.pair(c.ELL_MINUS, I, c.W_PLUS), 0)
        self.assertNotEqual(c.pair(c.ELL_PLUS, c.mul(R, I), c.W_MINUS), 0)

    def test_dimension_counts_and_three_direction_couplings(self):
        for d in range(1, 6):
            jet = c.ConnectionJet.constant(tuple(f'a{i}' for i in range(d)), (K,)*d)
            self.assertEqual(len(jet.components()), d*(d-1)//2)
        jet = c.ConnectionJet.constant(('r', 't1', 't2'), (K, R, RK))
        self.assertEqual(jet.curvature(0, 1), c.scale(RK, -2))
        self.assertEqual(jet.curvature(0, 2), c.scale(R, -2))
        self.assertEqual(jet.curvature(1, 2), c.scale(K, -2))
        self.assertEqual(jet.curvature(2, 0), c.scale(R, 2))

    def test_pullback_wedge_identity_with_derivative_terms(self):
        jet = c.ConnectionJet(('a', 'b', 'c'), (K, R, RK), ((I, K, R), (RK, Z, I), (K, R, Z)))
        j = ((1, 2), (Q(1, 3), -1), (2, Q(3, 5)))
        changed = jet.pullback(j, ('u', 'v'))
        expected = Z
        for i in range(3):
            for k in range(i+1, 3):
                expected = c.add(expected, c.scale(jet.curvature(i, k), j[i][0]*j[k][1]-j[k][0]*j[i][1]))
        self.assertEqual(changed.curvature(0, 1), expected)

    def test_one_direction_local_flatness_retains_lost_curvature(self):
        jet = c.ConnectionJet.constant(('a', 'b'), (K, R))
        self.assertNotEqual(jet.curvature(0, 1), Z)
        path = jet.pullback(((1,), (2,)), ('path',))
        self.assertEqual(path.components(), {})
        self.assertEqual(path.curvature(0, 0), Z)

    def test_invertible_base_change_preserves_full_local_flatness(self):
        jet = c.ConnectionJet.constant(('a', 'b'), (K, R))
        j = ((1, 2), (3, 5))
        changed = jet.pullback(j, ('u', 'v'))
        restored = changed.pullback(c.inverse(c.matrix(j)), ('a', 'b'))
        self.assertEqual(restored, jet)
        self.assertEqual(changed.curvature(0, 1), c.scale(jet.curvature(0, 1), -1))

    def test_constant_reference_change_curvature_covariance(self):
        jet = c.ConnectionJet(('a', 'b'), (K, R), ((I, RK), (K, Z)))
        s = c.matrix([[1, 2], [0, 1]])
        self.assertEqual(jet.change_reference(s).curvature(0, 1), c.conjugate(jet.curvature(0, 1), s))

    def test_constant_projector_full_connection_curvature_law(self):
        jet = c.ConnectionJet(('a', 'b'), (K, R), ((I, RK), (K, Z)))
        p = c.matrix([[1, 0], [0, 0]])
        expected = c.sub(c.mul(c.mul(p, jet.curvature(0, 1)), p),
                         c.compression_report(K, R, p)['excursions'])
        self.assertEqual(jet.compress(p).curvature(0, 1), expected)

    def test_noncommuting_coefficients_can_have_zero_connection_curvature(self):
        for t in (Q(-2), Q(-1, 3), Q(0), Q(2, 3), Q(3)):
            jet = c.derivative_cancellation_jet(t)
            self.assertNotEqual(c.commutator(*jet.operators), Z)
            self.assertEqual(jet.curvature(0, 1), Z)

    def test_constant_onsager_equivalence_and_orientation(self):
        for a, b, cc, d in itertools.product((-1, 0, 2), repeat=4):
            jet = c.linear_response_jet(((a, b), (cc, d)), (2, 3))
            self.assertEqual(jet.curvature(0, 1), c.matrix([[cc-b]]))
            self.assertEqual(c.is_zero(jet.curvature(0, 1)), b == cc)
            self.assertEqual(c.commutator(*jet.operators), c.matrix([[0]]))

    def test_symmetric_variable_response_is_not_sufficient(self):
        witness = c.exact_witnesses()['state_dependent_symmetric_response']
        l = witness['response']
        self.assertEqual(l[0][1], l[1][0])
        self.assertEqual(witness['curvature'], c.matrix([[-2]]))

    def test_existing_sheet_memory_survives_local_readout_flatness(self):
        # Reuse R7's independently declared integer sector; do not consume the
        # source EMKG3 compensator convention flagged by the master review.
        import emk_tensor_calculus as tc
        loop = tc.Move('p', 'p', (('V', I),), 2)
        audit = tc.return_report(loop, ((tc.Slot('V', 1),),))
        self.assertTrue(audit['full_carrier_return'])
        self.assertEqual(audit['sheet_residue'], 2)
        self.assertFalse(audit['full_unledgered_return'])
        self.assertTrue(tc.return_report(loop, ((tc.Slot('V', 1),),), sheet_ledger=2)['full_return_after_sheet_ledger'])

    def test_exact_source_curvature_agreement(self):
        source = Path(os.environ['RKF_SOURCE_ROOT'])/'proof_lab/emk_algebra_cut_graded_curvature.py'
        spec = importlib.util.spec_from_file_location('r8_pinned_source', source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for b, cc, d in itertools.product((Q(-1), Q(0), Q(2, 3)), repeat=3):
            a = c.scale(K, b)
            odd = c.add(c.scale(R, cc), c.scale(RK, d))
            source_result = module.commutator(a, odd)
            formula = c.scale(c.add(c.scale(R, d), c.scale(RK, cc)), -2*b)
            self.assertEqual(c.commutator(a, odd), source_result)
            self.assertEqual(source_result, formula)

    def test_contract_errors_are_rejected(self):
        with self.assertRaises(ValueError):
            c.compression_report(K, R, c.scale(I, 2))
        with self.assertRaises(ValueError):
            c.compression_report(K, R, c.identity(3))
        with self.assertRaises(ValueError):
            c.odd_readouts(K)
        with self.assertRaises(TypeError):
            c.recover_odd(0.5, 1)
        with self.assertRaises(ValueError):
            c.ConnectionJet.constant(('x', 'x'), (K, R))
        with self.assertRaises(ValueError):
            c.ConnectionJet.constant(('x',), (K,)).pullback(((1, 2), (3, 4)), ('u', 'v'))


if __name__ == '__main__':
    unittest.main()
