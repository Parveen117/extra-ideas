"""Exact R9 identities and source integration; run through verify_r9.py."""

import itertools
import unittest
import emk_curvature_balance as b

Q, K, R, RK = b.Q, b.K, b.R, b.RK
I, Z = b.identity(2), b.scale(b.identity(2), 0)
SOURCES = None


class CurvatureBalanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if SOURCES is None:
            raise RuntimeError('Run verify_r9.py to pin-check upstream sources before execution')
        cls.qth, cls.master = SOURCES['qth'], SOURCES['master']

    def test_cut_swap_is_the_existing_derived_master_operator(self):
        self.assertEqual(b.matrix(self.master.cut_swap()), K)
        self.assertEqual(b.mul(K, K), I)

    def test_exchange_invariance_selects_unique_equal_weight(self):
        for theta in (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)):
            observed = b.recognize(R, theta)
            self.assertEqual(observed == b.cut_action(observed), theta == Q(1, 2))
            self.assertEqual(b.recognize(b.cut_action(R), theta) == observed, theta == Q(1, 2))

    def test_balanced_observer_is_idempotent_and_preserves_even_sector(self):
        for a in (I, K, R, RK, b.add(K, b.scale(R, Q(2, 3)))):
            even, odd = b.grade(a)
            self.assertEqual(b.recognize(a), even)
            self.assertEqual(b.recognize(even), even)
            self.assertEqual(b.recognize(odd), Z)
            self.assertEqual(b.cut_action(even), even)

    def test_clock_free_iteration_closed_form(self):
        initial = b.add(K, b.add(b.scale(R, Q(2, 3)), RK))
        weights = (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))
        for sequence in itertools.product(weights, repeat=2):
            report = b.iterate(initial, sequence)
            self.assertEqual(report['closed_form_residual'], Z)
            self.assertEqual(report['odd_multiplier'], (1-2*sequence[0])*(1-2*sequence[1]))

    def test_finite_balance_requires_balancing_event_when_odd_content_present(self):
        for sequence in itertools.product((Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)), repeat=2):
            report = b.iterate(b.axis_state(Q(3, 5)), sequence)
            self.assertEqual(report['cut_balanced'], Q(1, 2) in sequence)
        self.assertTrue(b.iterate(b.axis_state(0), (Q(1, 4),))['cut_balanced'])

    def test_sign_flip_and_exact_contraction_without_clock(self):
        initial = b.axis_state(Q(3, 5))
        self.assertEqual(b.iterate(initial, (Q(3, 4),))['odd_multiplier'], Q(-1, 2))
        sequence = b.iterate(initial, (Q(1, 4),)*6)
        self.assertEqual(sequence['odd_multiplier'], Q(1, 64))
        self.assertFalse(sequence['cut_balanced'])

    def test_conditional_expectation_commutator_defect_grid(self):
        operators = (I, K, R, RK, b.add(K, R))
        for a, bb in itertools.product(operators, repeat=2):
            self.assertEqual(b.curvature_selection(a, bb)['identity_residual'], Z)

    def test_radial_tangential_curvature_is_cut_odd_and_cancels(self):
        for tangent, expected in ((R, b.scale(RK, -2)), (RK, b.scale(R, -2))):
            report = b.curvature_selection(K, tangent)
            self.assertEqual(report['raw'], expected)
            self.assertNotEqual(expected, Z)
            self.assertEqual(report['recognized_raw'], Z)

    def test_tangential_tangential_curvature_survives_balance(self):
        report = b.curvature_selection(R, RK)
        self.assertEqual(report['raw'], b.scale(K, -2))
        self.assertEqual(report['recognized_raw'], b.scale(K, -2))
        self.assertEqual(report['recomputed_from_recognized_generators'], Z)
        self.assertEqual(report['odd_pair_memory'], b.scale(K, -2))

    def test_cut_balanced_information_state_can_read_even_curvature(self):
        rho = b.scale(b.add(I, b.scale(K, Q(1, 2))), Q(1, 2))
        self.assertEqual(b.recognize(rho), rho)
        self.assertEqual(b.trace(b.mul(rho, b.commutator(K, R))), 0)
        self.assertEqual(b.trace(b.mul(rho, b.commutator(R, RK))), -1)

    def test_cut_covariance_in_generic_algebra_chart(self):
        from aghora_return import conjugate
        change = b.matrix([[1, 2], [0, 1]])
        kk = conjugate(K, change)
        for a in (K, R, RK):
            self.assertEqual(b.recognize(conjugate(a, change), cut=kk), conjugate(b.recognize(a), change))

    def test_signed_mean_variance_budget_is_exact(self):
        h = b.scale(b.commutator(K, R), Q(1, 2))
        for m, theta in itertools.product((Q(-3, 5), Q(-1, 3), Q(0), Q(1, 3), Q(3, 5)),
                                          (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))):
            before = b.curvature_budget(b.axis_state(m), h)
            after = b.curvature_budget(b.recognize(b.axis_state(m), theta), h)
            self.assertEqual(after['mean'], (1-2*theta)*m)
            self.assertEqual(after['second_moment'], 1)
            self.assertEqual(after['variance'], 1-after['mean']**2)
            self.assertEqual(after['variance']-before['variance'], 4*theta*(1-theta)*m*m)

    def test_trace_duality_of_state_and_observable_recognition(self):
        for a, f, theta in itertools.product((I, K, R, RK), (K, R, RK), (Q(1, 4), Q(1, 2), Q(3, 4))):
            self.assertEqual(b.trace(b.mul(b.recognize(a, theta), f)),
                             b.trace(b.mul(a, b.recognize(f, theta))))

    def test_source_master_mixed_channel_curvature_is_retained_in_audit(self):
        channels = [(K, Z), (Z, R)]
        raw = b.matrix(self.master.master_curvature(channels))
        self.assertEqual(raw, b.scale(RK, -2))
        self.assertEqual(tuple(b.matrix(x) for x in self.master.channel_curvatures(channels)), (Z, Z))
        self.assertEqual(b.matrix(self.master.mixed_curvature(channels)), raw)
        self.assertEqual(b.recognize(raw), Z)
        even_channels = [(R, Z), (Z, RK)]
        self.assertEqual(b.recognize(b.matrix(self.master.master_curvature(even_channels))), b.scale(K, -2))

    def test_source_full_family_pushforward_matches_analytic_tensor(self):
        for m, theta in itertools.product((Q(-3, 5), Q(-1, 3), Q(0), Q(1, 3), Q(3, 5)),
                                          (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1))):
            report = b.source_information_probe(self.qth, m, theta)
            self.assertEqual(report['source_information_after'], report['information_after'])
            self.assertEqual(report['source_curvature_after'], report['recomputed_tensor_curvature_after'])
            self.assertEqual(report['source_flagged_branch_information'], I)
            self.assertEqual(report['source_SLD_commutator_nonzero'], report['post_channel_SLD_pair_noncommutes'])

    def test_fixed_probe_and_recomputed_tensor_have_different_scaling(self):
        report = b.information_pushforward(Q(3, 5), Q(1, 4))
        self.assertEqual(report['fixed_probe_curvature_mean_after'], Q(3, 10))
        self.assertEqual(report['recomputed_tensor_curvature_after'], Q(3, 20))
        self.assertNotEqual(report['fixed_probe_curvature_mean_after'], report['recomputed_tensor_curvature_after'])

    def test_balancing_transforms_derivatives_and_records_information_cost(self):
        report = b.source_information_probe(self.qth, Q(3, 5), Q(1, 2))
        self.assertEqual(report['source_information_after'], b.matrix([[1, 0], [0, 0]]))
        self.assertEqual(report['discarded_information'], b.matrix([[0, 0], [0, 1]]))
        self.assertEqual(report['source_curvature_after'], 0)

    def test_recorded_cut_branch_allows_exact_state_recovery(self):
        rho = b.axis_state(Q(3, 5))
        for branch in (0, 1):
            out = rho if branch == 0 else b.cut_action(rho)
            recovered = out if branch == 0 else b.cut_action(out)
            self.assertEqual(recovered, rho)

    def test_averaging_without_branch_record_merges_distinct_states(self):
        a, bb = b.axis_state(Q(3, 5)), b.axis_state(Q(-3, 5))
        self.assertNotEqual(a, bb)
        self.assertEqual(b.recognize(a), b.recognize(bb))
        self.assertEqual(b.recognize(a), b.scale(I, Q(1, 2)))

    def test_finite_measurement_information_bound(self):
        directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (Q(3, 5), Q(4, 5), 0), (Q(3, 5), 0, Q(4, 5)))
        for m, n in itertools.product((Q(-3, 5), Q(-1, 3), Q(0), Q(1, 3), Q(3, 5)), directions):
            report = b.measurement_information(m, b.projective_effects(n))
            self.assertLessEqual(report['classical_trace'], 1)
            self.assertGreaterEqual(report['joint_trace_gap'], 1)
            self.assertEqual(sum(report['probabilities']), 1)

    def test_measurement_matrix_agrees_with_source_directional_information(self):
        directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (Q(3, 5), Q(4, 5), 0), (Q(3, 5), 0, Q(4, 5)))
        for m, n in itertools.product((Q(-3, 5), Q(-1, 3), Q(0), Q(1, 3), Q(3, 5)), directions):
            c = b.measurement_information(m, b.projective_effects(n))['classical_information']
            for dx, dy in ((1, 0), (0, 1), (1, 1)):
                predicted = dx*dx*c[0][0]+2*dx*dy*c[0][1]+dy*dy*c[1][1]
                self.assertEqual(predicted, self.qth.cfi_projective((Q(0), Q(0), m), (Q(dx), Q(dy), Q(0)), n))

    def test_balanced_mean_curvature_does_not_imply_joint_information_saturation(self):
        static = b.source_information_probe(self.qth, Q(0), Q(0))
        self.assertEqual(static['source_curvature_after'], 0)
        self.assertTrue(static['source_SLD_commutator_nonzero'])
        self.assertEqual(static['source_information_after'], I)
        report = b.measurement_information(0, b.projective_effects((1, 0, 0)))
        self.assertEqual(report['classical_information'], b.matrix([[1, 0], [0, 0]]))
        self.assertEqual(report['joint_trace_gap'], 1)

    def test_equal_axis_and_unsharp_measurement_ledgers(self):
        equal = b.measurement_information(0, ((Q(1, 4), (1, 0, 0)), (Q(1, 4), (-1, 0, 0)),
                                             (Q(1, 4), (0, 1, 0)), (Q(1, 4), (0, -1, 0))))
        self.assertEqual(equal['classical_information'], b.scale(I, Q(1, 2)))
        self.assertEqual(equal['discarded_information'], b.scale(I, Q(1, 2)))
        unsharp = b.measurement_information(Q(1, 3), ((Q(1, 2), (Q(1, 2), 0, 0)),
                                                     (Q(1, 2), (Q(-1, 2), 0, 0))))
        self.assertEqual(unsharp['classical_trace'], Q(1, 4))

    def test_invalid_balance_and_measurement_contracts_are_rejected(self):
        with self.assertRaises(TypeError):
            b.recognize(K, 0.5)
        with self.assertRaises(ValueError):
            b.recognize(K, Q(5, 4))
        with self.assertRaises(ValueError):
            b.recognize(K, cut=b.scale(K, 2))
        with self.assertRaises(ValueError):
            b.information_pushforward(1, Q(1, 2))
        with self.assertRaises(ValueError):
            b.curvature_budget(b.axis_state(0), R)
        with self.assertRaises(ValueError):
            b.measurement_information(0, ((Q(1), (1, 0, 0)),))
        with self.assertRaises(ValueError):
            b.measurement_information(0, ((Q(1, 2), (2, 0, 0)), (Q(1, 2), (-2, 0, 0))))
