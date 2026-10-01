"""Exact R13 complement, memory, event law and information checks."""

import itertools
import unittest
import complement_memory_noise as m

n, Q = m.n, m.Q
SOURCES = None
PARAMETERS = (Q(-2), Q(-3, 4), Q(-2, 3), Q(-1, 2), Q(-1, 3), Q(-1, 5),
              Q(1, 5), Q(1, 3), Q(1, 2), Q(2, 3), Q(3, 4), Q(2))


class ComplementMemoryNoiseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if SOURCES is None:
            raise RuntimeError('Use verify_r13.py to check the unchanged public sources')
        cls.carrier, cls.quantum = SOURCES['carrier'], SOURCES['quantum']

    def test_rotation_is_the_existing_native_even_plus_odd_carrier(self):
        for parameter in PARAMETERS:
            data = m.cayley_rotation(parameter)
            c, s, t = data['cosine'], data['sine'], data['transport']
            self.assertEqual(t, n.matrix(self.carrier.block(c, 0, s, 0)))
            self.assertEqual(self.carrier.coeffs([list(row) for row in t]), (c, 0, s, 0))
            self.assertEqual(n.mul(n.transpose(t), t), n.identity(2))

    def test_observer_cut_is_derived_from_the_existing_RK_basis(self):
        j = n.scale(n.RK, -1)
        p = n.scale(n.add(n.identity(2), j), Q(1, 2))
        q = n.sub(n.identity(2), p)
        self.assertEqual(p, n.matrix(((1, 0), (0, 0))))
        self.assertEqual(n.mul(j, j), n.identity(2))
        self.assertEqual(n.mul(n.mul(j, n.R), j), n.scale(n.R, -1))
        self.assertTrue(n.is_zero(n.mul(n.mul(p, n.R), p)))
        self.assertEqual(n.mul(n.mul(n.mul(n.mul(p, n.R), q), n.R), p), n.scale(p, -1))

    def test_cut_noncommutation_and_complement_strength_have_the_same_native_amplitude(self):
        for parameter in PARAMETERS:
            data = m.cayley_rotation(parameter)
            c, s = data['cosine'], data['sine']
            result = m.cut_incompatibility(c, s)
            self.assertEqual(result['cut_transport_commutator'], n.scale(n.K, -s))
            self.assertEqual(result['commutator_square'], n.scale(n.identity(2), s*s))
            self.assertEqual(result['complementary_norm'], result['retained_norm_loss'])

    def test_complement_memory_is_distinct_from_curvature_of_commuting_transports(self):
        t1 = m.rotation(Q(3, 5), Q(4, 5))
        t2 = m.rotation(Q(5, 13), Q(12, 13))
        self.assertTrue(n.is_zero(n.commutator(t1, t2)))
        self.assertFalse(n.is_zero(m.cut_incompatibility(Q(3, 5), Q(4, 5))['cut_transport_commutator']))
        self.assertNotEqual(m.memory_kernels(t1, 1, 1), (((Q(0),),),))

    def test_hidden_elimination_reproduces_full_rotation_histories(self):
        for parameter, x, z in itertools.product(PARAMETERS[:4], (Q(-1), Q(0), Q(2)),
                                                  (Q(-2), Q(0), Q(3))):
            t = m.cayley_rotation(parameter)['transport']
            history = m.reduced_history(t, 1, ((x,),), ((z,),), 7)
            self.assertTrue(all(n.is_zero(r) for r in history['residuals']))
            c, s = t[0][0], t[1][0]
            self.assertEqual(history['complement_forces'], tuple(((-s*(c**k)*z,),) for k in range(7)))
            self.assertEqual(history['kernels'], tuple(((-s*s*(c**k),),) for k in range(7)))

    def test_memory_term_is_load_bearing_even_when_the_initial_hidden_state_is_zero(self):
        t = m.rotation(Q(3, 5), Q(4, 5))
        result = m.reduced_history(t, 1, ((1,),), ((0,),), 3)
        self.assertEqual(result['visible_states'][2], ((Q(-7, 25),),))
        self.assertEqual(result['memory_terms'][1], ((Q(-16, 25),),))
        self.assertTrue(all(n.is_zero(r) for r in result['residuals']))
        self.assertTrue(all(n.is_zero(r) for r in result['complement_forces']))

    def test_two_visible_time_records_close_the_second_order_native_recurrence(self):
        for parameter in PARAMETERS:
            t = m.cayley_rotation(parameter)['transport']
            c = t[0][0]
            result = m.reduced_history(t, 1, ((Q(5, 13),),), ((Q(12, 13),),), 6)
            xs = [x[0][0] for x in result['visible_states']]
            for k in range(5):
                self.assertEqual(xs[k+2], 2*c*xs[k+1]-xs[k])

    def test_one_visible_value_does_not_determine_the_complementary_force(self):
        a = m.reduced_history(m.rotation(Q(3, 5), Q(4, 5)), 1,
                              ((Q(5, 13),),), ((Q(12, 13),),), 1)
        b = m.reduced_history(m.rotation(Q(3, 5), Q(4, 5)), 1,
                              ((Q(5, 13),),), ((Q(-12, 13),),), 1)
        self.assertEqual(a['visible_states'][0], b['visible_states'][0])
        self.assertNotEqual(a['visible_states'][1], b['visible_states'][1])

    def test_a_future_visible_probe_recovers_the_hidden_initial_component(self):
        for parameter, x, z in itertools.product(PARAMETERS, (Q(-1), Q(0), Q(2)),
                                                  (Q(-2), Q(0), Q(3))):
            t = m.cayley_rotation(parameter)['transport']
            c, s = t[0][0], t[1][0]
            x1 = c*x-s*z
            self.assertEqual(m.recover_complement(c, s, x, x1), z)

    def test_complement_recovery_error_bound_is_sharp_and_gain_dependent(self):
        c, s, e0, e1 = Q(3, 5), Q(4, 5), Q(1, 20), Q(1, 30)
        radius = m.complement_error_radius(c, s, e0, e1)
        errors = [abs(m.recover_complement(c, s, dx, dy))
                  for dx, dy in itertools.product((-e0, e0), (-e1, e1))]
        self.assertEqual(radius, max(errors))
        tiny = m.cayley_rotation(Q(1, 100))
        self.assertGreater(m.complement_error_radius(tiny['cosine'], tiny['sine'], e0, e1), radius)

    def test_balanced_unit_preparation_derives_the_hidden_variance(self):
        for x, z in ((Q(3, 5), Q(4, 5)), (Q(5, 13), Q(12, 13)), (Q(1), Q(0))):
            p = m.balanced_unit_pair(x, z)
            self.assertEqual(p['hidden_mean'], 0)
            self.assertEqual(p['hidden_variance'], z*z)
            self.assertEqual(sum(w*y*y for w, y in zip(p['weights'], p['hidden_branches'])), 1-x*x)

    def test_noise_covariance_is_derived_from_the_same_complementary_leg(self):
        result = m.endogenous_noise(Q(3, 5), Q(4, 5), Q(5, 13), Q(12, 13), 5)
        expected = Q(16, 25)*Q(144, 169)
        self.assertEqual(result['covariance'][0][0][0][0], expected)
        for k, l in itertools.product(range(5), repeat=2):
            self.assertEqual(result['covariance'][k][l][0][0], expected*(Q(3, 5)**(k+l)))
            self.assertEqual(result['covariance_by_branches'][k][l], result['covariance'][k][l][0][0])
        self.assertEqual(result['noise_mean'], (Q(0),)*5)

    def test_memory_and_noise_have_an_exact_shared_kernel_relation(self):
        for parameter in PARAMETERS:
            t = m.cayley_rotation(parameter)['transport']
            result = m.endogenous_noise(t[0][0], t[1][0], Q(5, 13), Q(12, 13), 5)
            self.assertEqual(result['memory_noise_residual'], (Q(0),)*5)

    def test_noise_is_colored_and_temporally_rank_one(self):
        result = m.endogenous_noise(Q(3, 5), Q(4, 5), Q(5, 13), Q(12, 13), 5)
        self.assertNotEqual(result['covariance'][0][1][0][0], 0)
        self.assertEqual(result['temporal_rank'], 1)
        self.assertTrue(m.o.positive_semidefinite(result['covariance_by_branches']))

    def test_balanced_complement_does_not_generate_a_Gaussian_likelihood(self):
        p = m.balanced_unit_pair(Q(5, 13), Q(12, 13))
        self.assertFalse(p['gaussian_hidden_law'])
        self.assertEqual(p['hidden_fourth_moment'], p['hidden_variance']**2)
        self.assertNotEqual(p['hidden_fourth_moment'], 3*p['hidden_variance']**2)

    def test_general_block_elimination_keeps_multiple_hidden_modes_and_feedback(self):
        t = n.matrix([[Q(3, 5), Q(-4, 5), 0], [Q(4, 5), Q(3, 5), 1], [0, 0, 1]])
        result = m.reduced_history(t, 1, ((2,),), ((-1,), (3,)), 8)
        self.assertTrue(all(n.is_zero(r) for r in result['residuals']))
        self.assertNotEqual(result['complement_forces'][1], ((Q(12, 25),),))
        cov = m.complementary_covariance(t, 1, ((2, 1), (1, 2)), 5)
        scalar_cov = tuple(tuple(item[0][0] for item in row) for row in cov)
        self.assertTrue(m.o.positive_semidefinite(scalar_cov))
        self.assertEqual(len(m.o.independent_rows(scalar_cov)), 2)

    def test_a_nonreturning_hidden_record_generates_no_visible_noise_or_memory(self):
        t = n.matrix([[1, 0], [1, 1]])
        result = m.reduced_history(t, 1, ((2,),), ((3,),), 4)
        self.assertTrue(all(n.is_zero(k) for k in result['kernels']))
        self.assertTrue(all(n.is_zero(k) for k in result['complement_forces']))
        self.assertEqual(result['visible_states'], (((Q(2),),),)*5)
        self.assertNotEqual(result['full_states'][0], result['full_states'][-1])

    def test_first_exit_weights_are_norms_of_actual_native_branch_amplitudes(self):
        p, q = n.matrix(((1, 0), (0, 0))), n.matrix(((0, 0), (0, 1)))
        for parameter in PARAMETERS:
            data = m.cayley_rotation(parameter)
            c, s, t = data['cosine'], data['sine'], data['transport']
            result = m.first_exit(c, s, 5)
            survival = n.mul(n.mul(p, t), p)
            for k, expected in enumerate(result['probabilities']):
                branch = n.mul(n.mul(n.mul(q, t), m.power(survival, k)), p)
                branch_norm = n.mul(n.transpose(branch), branch)
                self.assertEqual(branch_norm, n.scale(p, expected))

    def test_the_finite_event_ledger_retains_its_exact_unseen_tail(self):
        for c, s, length in itertools.product((Q(3, 5),), (Q(4, 5),), range(7)):
            result = m.first_exit(c, s, length)
            self.assertEqual(result['total_mass_with_tail'], 1)
            self.assertEqual(result['unresolved_tail'], Q(9, 25)**(length+1))
            self.assertGreater(result['unresolved_tail'], 0)

    def test_event_degeneracy_distinguishes_immediate_exit_and_no_exit(self):
        immediate = m.first_exit(0, 1, 3)
        self.assertEqual(immediate['probabilities'], (Q(1), Q(0), Q(0), Q(0)))
        self.assertEqual(immediate['unresolved_tail'], 0)
        closed = m.first_exit(1, 0, 3)
        self.assertFalse(closed['terminates_almost_surely'])
        self.assertEqual(closed['unresolved_tail'], 1)
        self.assertIsNone(closed['mean_continuations'])
        for c, s in ((0, 1), (1, 0)):
            with self.assertRaises(ValueError):
                m.event_information(c, s)

    def test_varying_native_moves_derive_a_nonconstant_event_hazard(self):
        angles = ((Q(3, 5), Q(4, 5)), (Q(5, 13), Q(12, 13)), (Q(4, 5), Q(3, 5)))
        result = m.varying_first_exit(angles)
        self.assertEqual(result['probabilities'][0], Q(16, 25))
        self.assertEqual(result['probabilities'][1], Q(9, 25)*Q(144, 169))
        self.assertEqual(result['total_mass_with_tail'], 1)
        reversed_result = m.varying_first_exit(tuple(reversed(angles)))
        self.assertEqual(result['unresolved_tail'], reversed_result['unresolved_tail'])
        self.assertNotEqual(result['probabilities'], reversed_result['probabilities'])

    def test_coherent_retention_and_monitored_reset_have_different_return_weights(self):
        result = m.coherent_vs_monitored(Q(3, 5), Q(4, 5), 2)
        self.assertEqual(result['coherent_visible_amplitude'], Q(-7, 25))
        self.assertEqual(result['monitored_survival_amplitude'], Q(9, 25))
        self.assertEqual(result['coherent_final_visible_weight'], Q(49, 625))
        self.assertEqual(result['monitored_survival_weight'], Q(81, 625))
        self.assertEqual(result['two_step_compression_defect'], result['two_step_excursion'])
        self.assertEqual(result['two_step_excursion_coefficient'], Q(-16, 25))

    def test_binary_and_completed_event_information_match_the_unchanged_QTH_readout(self):
        for parameter in PARAMETERS:
            data = m.cayley_rotation(parameter)
            c, s = data['cosine'], data['sine']
            info = m.event_information(c, s)
            bloch = (2*c*s, Q(0), c*c-s*s)
            angular_derivative = (2*(c*c-s*s), Q(0), -4*c*s)
            source = self.quantum.cfi_projective(bloch, angular_derivative, (Q(0), Q(0), Q(1)))
            self.assertEqual(info['binary_fisher_angular'], source)
            self.assertEqual(source, 4)
            self.assertEqual(info['waiting_fisher_angular'], 4/(s*s))
            self.assertEqual(info['information_per_expected_trial'], 4)

    def test_completed_waiting_moments_are_checked_against_finite_tail_composition(self):
        r, q = Q(9, 25), Q(16, 25)
        mean, variance = r/q, r/(q*q)
        second = variance+mean*mean
        for length in range(5):
            event = m.first_exit(Q(3, 5), Q(4, 5), length)
            tail = event['unresolved_tail']
            mu = sum(k*p for k, p in enumerate(event['probabilities']))+tail*(length+1+mean)
            moment = sum(k*k*p for k, p in enumerate(event['probabilities']))+tail*(
                (length+1)**2+2*(length+1)*mean+second)
            self.assertEqual(mu, mean)
            self.assertEqual(moment-mu*mu, variance)
            info = m.event_information(Q(3, 5), Q(4, 5))
            self.assertEqual(variance/(r*r), info['waiting_fisher_continuation'])

    def test_probability_fluctuations_and_hidden_force_covariance_are_distinct_derived_objects(self):
        weights = m.branch_weights(Q(3, 5), Q(4, 5))
        self.assertEqual(weights['binary_event_variance'], Q(144, 625))
        noise = m.endogenous_noise(Q(3, 5), Q(4, 5), Q(5, 13), Q(12, 13), 1)
        self.assertNotEqual(weights['binary_event_variance'], noise['covariance'][0][0][0][0])
        self.assertEqual(weights['native_split_norm'], 1)

    def test_r12_event_bridge_eliminates_the_independent_stopping_and_gain_parameters(self):
        report = m.r12_event_bridge(Q(3, 5), Q(4, 5), Q(1, 3))
        self.assertEqual(report['derived_continuation'], Q(9, 25))
        self.assertEqual(report['derived_record_gain'], Q(4, 5))
        self.assertEqual(report['derived_kappa'], Q(9, 25))
        self.assertEqual(report['tangent_metric'], n.matrix(((Q(26, 25), 0), (0, 1))))
        self.assertTrue(report['Gaussian_calibration_retained_from_R12'])
        self.assertFalse(report['Gaussian_calibration_derived_from_complement'])

    def test_native_balanced_binary_meter_generates_its_own_probability_and_noise(self):
        for amplitude in (Q(-4), Q(-1), Q(0), Q(1, 3), Q(2)):
            meter = m.native_binary_response(amplitude)
            p0, p1 = meter['probabilities']
            self.assertEqual(p0+p1, 1)
            self.assertGreaterEqual(min(p0, p1), 0)
            self.assertEqual(meter['contrast_variance'], 4*p0*p1)
            rho0 = self.quantum.bloch_state((Q(1), Q(0), Q(0)))
            transport = m.cayley_rotation(amplitude/4)['transport']
            matrix0 = n.matrix([[pair[0] for pair in row] for row in rho0])
            rho = m.o.product(transport, matrix0, n.transpose(transport))
            self.assertEqual((rho[0][0], rho[1][1]), meter['probabilities'])
        origin = m.native_binary_response(0)
        self.assertEqual(origin['probabilities'], (Q(1, 2), Q(1, 2)))
        self.assertEqual(origin['contrast_variance'], 1)

    def test_native_binary_score_metric_matches_independent_QTH_score_evaluations(self):
        rows = ((Q(1), Q(2)), (Q(-1), Q(3)), (Q(1, 2), Q(-1, 3)))
        result = m.native_binary_score_metric(rows)
        self.assertFalse(result['Gaussian_assumption'])
        self.assertEqual(result['record_covariance_at_origin'], n.identity(3))
        for h in itertools.product((Q(-1), Q(0), Q(2)), repeat=2):
            value = sum(self.quantum.cfi_projective((Q(1), Q(0), Q(0)),
                (Q(0), Q(0), -sum(a*b for a, b in zip(row, h))), (Q(0), Q(0), Q(1)))
                for row in rows)
            column = tuple((x,) for x in h)
            self.assertEqual(value, m.o.product(n.transpose(column), result['metric'], column)[0][0])

    def test_a_native_binary_likelihood_realizes_the_seam_metric_without_Gaussian_noise(self):
        for parameter in PARAMETERS:
            data = m.cayley_rotation(parameter)
            c, s, v = data['cosine'], data['sine'], Q(1, 3)
            result = m.native_binary_seam(c, s, v)
            kappa = c*c
            self.assertFalse(result['Gaussian_assumption'])
            self.assertFalse(result['external_record_noise_precision'])
            self.assertEqual(result['derived_kappa'], kappa)
            self.assertEqual(result['tangent_metric'], n.matrix(((1+kappa*v*v, 0), (0, 1))))
            self.assertEqual(result['future_observer'], n.projection(3, 4))
            self.assertEqual(result['gate_information']['waiting_fisher_angular'], 4/(s*s))
            self.assertEqual(result['gate_response_cross_score_at_origin'], (Q(0),)*4)
            self.assertEqual(result['untagged_response_metric_at_origin'], result['seed_information'])
            self.assertTrue(m.o.positive_semidefinite(result['discarded_tag_information']))
            self.assertEqual(result['tangent_metric_first'][1], n.matrix(((2*kappa*v, 0), (0, 0))))
            self.assertEqual(result['tangent_metric_second'][1][1], n.matrix(((2*kappa, 0), (0, 0))))

    def test_event_norms_and_noise_covariance_do_not_recover_orientation(self):
        positive = m.first_exit(Q(3, 5), Q(4, 5), 3)
        negative = m.first_exit(Q(3, 5), Q(-4, 5), 3)
        self.assertEqual(positive, negative)
        a = m.endogenous_noise(Q(3, 5), Q(4, 5), Q(5, 13), Q(12, 13), 3)
        b = m.endogenous_noise(Q(3, 5), Q(-4, 5), Q(5, 13), Q(12, 13), 3)
        self.assertEqual(a['covariance'], b['covariance'])
        self.assertNotEqual(a['noise_paths'], b['noise_paths'])

    def test_invalid_preparations_covariances_and_noncoupled_recovery_are_rejected(self):
        with self.assertRaises(ValueError):
            m.rotation(Q(1, 2), Q(1, 2))
        with self.assertRaises(ValueError):
            m.balanced_unit_pair(1, 1)
        with self.assertRaises(ValueError):
            m.recover_complement(1, 0, 1, 1)
        with self.assertRaises(ValueError):
            m.complementary_covariance(m.rotation(Q(3, 5), Q(4, 5)), 1, ((-1,),), 2)
        with self.assertRaises(TypeError):
            m.cayley_rotation(0.5)
        with self.assertRaises(ValueError):
            m.first_exit(0, 1, -1)
