"""Exact tests of R12 selection, forgetting, covariance and derived geometry."""

from dataclasses import replace
import itertools
import unittest
import observer_metric_foundation as o
from emk_curvature_balance import recognize

n, Q = o.n, o.Q
SOURCES = None


def diagonal(*values):
    return n.matrix([[value if i == j else 0 for j in range(len(values))]
                     for i, value in enumerate(values)])


def diagonal_experiment():
    return o.EventExperiment(((1, 1),), ((1,),), (diagonal(1, Q(1, 2)),),
                             (Q(1),), Q(1, 2), n.identity(2), Q(3, 4))


class ObserverMetricFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if SOURCES is None:
            raise RuntimeError('Use verify_r12.py to check upstream source bytes before execution')
        cls.geometry = SOURCES['geometry']

    def test_dynamic_metric_is_selected_by_the_entire_event_law(self):
        report = o.derive_observer_metric(diagonal_experiment())
        self.assertEqual(report['metric'], n.matrix([[1, Q(2, 3)], [Q(2, 3), Q(4, 7)]]))
        self.assertTrue(n.is_zero(report['stein_residual']))

    def test_future_observer_is_derived_beyond_the_seed_probe(self):
        report = o.derive_observer_metric(diagonal_experiment())
        self.assertEqual(len(o.independent_rows(diagonal_experiment().readout)), 1)
        self.assertEqual(report['rank'], 2)
        self.assertEqual(report['observer'], n.identity(2))
        self.assertEqual(report['kernel'], ())

    def test_feedback_markers_select_four_state_modes_from_two_seed_modes(self):
        e = o.feedback_experiment()
        result = o.derive_observer_metric(e)
        self.assertEqual(result['rank'], 4)
        self.assertEqual(result['metric'], diagonal(Q(64, 63), Q(64, 63), Q(8, 63), Q(8, 63)))
        for action, t in zip(result['descended_transports'], e.transports):
            self.assertEqual(o.product(action, result['observer']), o.product(result['observer'], t))

    def test_positive_probability_and_noise_changes_keep_the_future_kernel(self):
        base = diagonal_experiment()
        for r, gain in itertools.product((Q(1, 4), Q(1, 2), Q(2, 3)), (Q(1, 3), Q(1), Q(7))):
            e = replace(base, continuation=r, precision=((gain,),), contraction=(1+r)/2)
            report = o.derive_observer_metric(e)
            self.assertEqual(report['observer'], n.identity(2))
            self.assertTrue(o.positive_definite(report['metric']))

    def test_zero_readout_selects_the_zero_observer(self):
        e = replace(diagonal_experiment(), readout=((0, 0),))
        result = o.derive_observer_metric(e)
        self.assertEqual(result['rank'], 0)
        self.assertEqual(result['observer'], ())
        self.assertEqual(len(result['kernel']), 2)
        self.assertTrue(n.is_zero(result['metric']))

    def test_quotient_factorization_is_independent_of_the_section(self):
        e = o.seam_experiment(2, Q(1, 3), Q(1, 2))
        result = o.derive_observer_metric(e)
        c, b = result['observer'], result['section']
        alternate = tuple(tuple(x+(j+1 if i == 3 else 0) for j, x in enumerate(row))
                          for i, row in enumerate(b))
        self.assertEqual(o.product(c, alternate), n.identity(3))
        self.assertEqual(o.product(n.transpose(alternate), result['metric'], alternate),
                         result['quotient_metric'])
        self.assertEqual(o.product(n.transpose(c), result['quotient_metric'], c), result['metric'])

    def test_state_reference_change_carries_the_metric_covariantly(self):
        e = o.seam_experiment(2, Q(1, 3), Q(1, 2))
        s = n.matrix([[1, 1, 0, 1], [0, 1, 1, 0], [0, 0, 2, 0], [0, 0, 0, 1]])
        si = n.inverse(s)
        transformed = replace(e, readout=o.product(e.readout, si),
                              transports=tuple(o.product(s, t, si) for t in e.transports),
                              stability_form=o.product(n.transpose(si), e.stability_form, si))
        a, b = o.derive_observer_metric(e), o.derive_observer_metric(transformed)
        self.assertEqual(b['metric'], o.product(n.transpose(si), a['metric'], si))
        self.assertEqual(b['rank'], a['rank'])
        self.assertTrue(n.is_zero(o.product(b['observer'], o.product(s, ((0,), (0,), (0,), (1,))))))

    def test_record_coordinate_change_is_a_calibration_change(self):
        e = o.seam_experiment(2, Q(1, 3), Q(1, 2))
        u = n.matrix([[1, 1, 0], [0, 2, 0], [1, 0, 1]])
        ui = n.inverse(u)
        transformed = replace(e, readout=o.product(u, e.readout),
                              precision=o.product(n.transpose(ui), e.precision, ui))
        self.assertEqual(o.tagged_metric(transformed), o.tagged_metric(e))

    def test_finite_event_apertures_have_a_positive_exact_remainder(self):
        e = diagonal_experiment()
        g = o.tagged_metric(e)
        for length in range(5):
            finite = o.finite_tagged_metric(e, length)
            remainder = n.sub(g, finite)
            predicted = g
            for _ in range(length+1):
                predicted = o.event_operator(e, predicted)
            self.assertEqual(remainder, predicted)
            self.assertTrue(o.positive_semidefinite(remainder))
            self.assertFalse(n.is_zero(remainder))

    def test_stability_is_load_bearing_and_an_algebraic_solve_is_not_enough(self):
        with self.assertRaisesRegex(ValueError, 'stability inequality'):
            o.EventExperiment(((1,),), ((1,),), (((2,),),), (Q(1),),
                              Q(3, 4), ((1,),), Q(7, 8))
        with self.assertRaises(ValueError):
            replace(diagonal_experiment(), contraction=Q(1))

    def test_precision_must_be_positive_and_transport_must_be_lawful(self):
        for precision in (((0,),), ((-1,),)):
            with self.assertRaises(ValueError):
                replace(diagonal_experiment(), precision=precision)
        with self.assertRaises(ValueError):
            replace(diagonal_experiment(), transports=(diagonal(1, 0),))
        with self.assertRaises(ValueError):
            replace(diagonal_experiment(), continuation=Q(1))

    def test_catalogue_support_and_rational_input_are_explicit(self):
        with self.assertRaises(ValueError):
            replace(diagonal_experiment(), probabilities=(Q(0),))
        with self.assertRaises(TypeError):
            o.seam_experiment(0.5, 1, Q(1, 2))
        with self.assertRaises(ValueError):
            o.rect(((1, 0), (1,)))

    def test_singular_positive_forms_and_zero_pivots_are_handled_exactly(self):
        self.assertTrue(o.positive_semidefinite(n.matrix([[0, 0], [0, 2]])))
        self.assertTrue(o.positive_semidefinite(n.matrix([[1, 1], [1, 1]])))
        self.assertFalse(o.positive_definite(n.matrix([[1, 1], [1, 1]])))
        self.assertFalse(o.positive_semidefinite(n.matrix([[0, 1], [1, 0]])))
        self.assertFalse(o.positive_semidefinite(n.matrix([[1, 2], [2, 1]])))

    def test_balanced_sheet_metric_is_derived_instead_of_declared(self):
        for c, v, r in itertools.product((Q(0), Q(1), Q(-2)),
                                         (Q(0), Q(1, 3), Q(-3, 4)),
                                         (Q(1, 4), Q(1, 2), Q(2, 3))):
            e = o.seam_experiment(c, v, r)
            kappa = c*c*r/(1-r)
            report = o.derive_observer_metric(e)
            self.assertEqual(report['metric'], diagonal(1+kappa*v*v, 1, 1, 0))
            self.assertEqual(report['observer'], n.projection(3, 4))

    def test_derived_metric_two_jet_contains_the_hidden_second_moment(self):
        result = o.seam_foundation(2, Q(1, 3), Q(1, 2))
        g = result['metric_jet']
        z = n.scale(n.identity(2), 0)
        self.assertEqual(g.value, diagonal(Q(13, 9), 1))
        self.assertEqual(g.first, (z, diagonal(Q(8, 3), 0)))
        self.assertEqual(g.second, ((z, z), (z, diagonal(8, 0))))

    def test_two_jet_differentiates_every_cross_and_second_transport_term(self):
        e = o.EventExperiment(((1,),), ((1,),), (((Q(1, 2),),),),
                              (Q(1),), Q(1, 4), ((1,),), Q(1, 2))
        first = ((((2,),),), (((3,),),))
        second = (((((0,),),), (((4,),),)), ((((4,),),), (((0,),),)))
        jet = o.response_two_jet(e, first, second)
        self.assertEqual(jet['value'], ((Q(4, 5),),))
        self.assertEqual(jet['first'], (((Q(32, 75),),), ((Q(16, 25),),)))
        self.assertEqual(jet['second'][0][1], ((Q(512, 125),),))
        self.assertEqual(jet['second'][0][0], ((Q(2432, 1125),),))
        self.assertEqual(jet['second'][1][1], ((Q(608, 125),),))

    def test_seed_derivatives_enter_the_same_selection_equation(self):
        e = o.EventExperiment(((1,),), ((1,),), (((1,),),),
                              (Q(1),), Q(1, 2), ((1,),), Q(3, 4))
        jet = o.response_two_jet(e, ((((0,),),),), (((((0,),),),),),
                                 first_seed=(((2,),),), second_seed=((((2,),),),))
        self.assertEqual(jet['first'], (((Q(2),),),))
        self.assertEqual(jet['second'], ((((Q(2),),),),))

    def test_inconsistent_jets_and_invisible_tangent_modes_are_rejected(self):
        e = o.seam_experiment(1, 1, Q(1, 2))
        z = n.scale(n.identity(4), 0)
        with self.assertRaises(ValueError):
            o.response_two_jet(e, ((z, z),), ())
        with self.assertRaises(ValueError):
            o.tangent_metric({'value': o.tagged_metric(e), 'first': (z,), 'second': ((z,),)},
                              ((0,), (0,), (0,), (1,)), ('u',))

    def test_derived_christoffel_and_curvature_match_the_existing_metric_source(self):
        for c, v, r in itertools.product((Q(0), Q(1), Q(2)),
                                         (Q(0), Q(1, 3), Q(-3, 4)),
                                         (Q(1, 3), Q(1, 2))):
            result = o.seam_foundation(c, v, r)
            kappa, g, lc = result['kappa'], result['metric_jet'], result['levi_civita']
            a = self.geometry.christoffel_u_uv(v, kappa)
            b = self.geometry.christoffel_v_uu(v, kappa)
            self.assertEqual(lc.operators, (n.matrix([[0, a], [b, 0]]), diagonal(a, 0)))
            lowered = o.product(g.value, lc.curvature(0, 1))
            self.assertEqual(lowered[0][1], self.geometry.riemann_uvuv(v, kappa))
            self.assertEqual(lowered[0][1]/g.value[0][0], self.geometry.gaussian_curvature_closed(v, kappa))

    def test_balanced_native_curvature_can_vanish_while_derived_metric_is_curved(self):
        s = o.seam_foundation(2, 0, Q(1, 2))
        self.assertFalse(n.is_zero(s['raw_native_curvature']))
        self.assertTrue(n.is_zero(recognize(s['raw_native_curvature'], cut=o.spectral.cut(2))))
        self.assertFalse(n.is_zero(s['levi_civita'].curvature(0, 1)))
        self.assertEqual(s['kappa'], 4)

    def test_event_tag_forgetting_has_a_positive_information_ledger_at_the_origin(self):
        s = o.seam_foundation(2, Q(1, 3), Q(1, 2))
        report = s['untagged_at_origin']
        self.assertEqual(report['mean_transport'], n.identity(4))
        self.assertEqual(report['metric'], diagonal(1, 1, 1, 0))
        self.assertEqual(report['discard'], diagonal(Q(4, 9), 0, 0, 0))
        self.assertEqual(report['state_scope'], 'zero state only')

    def test_forgetting_the_history_tag_can_destroy_future_complete_observation(self):
        e = diagonal_experiment()
        result = o.origin_untagged_metric(e)
        self.assertEqual(result['mean_response'], ((Q(1), Q(2, 3)),))
        self.assertEqual(result['metric'], n.matrix([[1, Q(2, 3)], [Q(2, 3), Q(4, 9)]]))
        self.assertEqual(len(o.independent_rows(result['metric'])), 1)
        self.assertEqual(o.derive_observer_metric(e)['rank'], 2)
        self.assertTrue(o.positive_semidefinite(result['discard']))

    def test_marginal_gaussian_score_is_the_average_score_only_at_a_common_mean(self):
        # At the origin every conditional density is the same. Two tagged
        # responses +/-1 have Fisher 1; their marginal first derivative is 0.
        responses = (Q(1), Q(-1))
        tagged = sum(x*x for x in responses)/2
        marginal_score = sum(responses)/2
        self.assertEqual(tagged, 1)
        self.assertEqual(marginal_score*marginal_score, 0)
        # This control forbids substituting the average of squared scores.
        self.assertNotEqual(tagged, marginal_score*marginal_score)

    def test_native_connection_need_not_be_the_derived_metrics_levi_civita_connection(self):
        s = o.seam_foundation(2, Q(1, 3), Q(1, 2))
        d = s['native_visible_decomposition']
        self.assertFalse(s['native_connection_is_levi_civita'])
        self.assertTrue(n.is_zero(d['visible_full_curvature']))
        self.assertFalse(n.is_zero(d['Riemann']))
        self.assertEqual(d['distortion'], n.scale(d['Riemann'], -1))
        self.assertTrue(n.is_zero(d['excursion']))
        self.assertTrue(n.is_zero(d['reconstruction_residual']))

    def test_the_derived_metric_can_supply_an_exact_riemann_special_sector(self):
        s = o.seam_foundation()
        lc, g = s['levi_civita'], s['metric_jet']
        native = n.extend_connection(lc, n.ConnectionJet.constant(('u', 'v'), (n.K, n.R)))
        self.assertTrue(n.riemann_reduction(native, g, n.projection(2, 4))['Riemann_sector_certified'])
        self.assertFalse(n.is_zero(native.curvature(0, 1)))

    def test_the_event_law_retains_a_physical_selection_freedom(self):
        first = o.seam_foundation(1, 1, Q(1, 2))
        second = o.seam_foundation(1, 1, Q(1, 3))
        self.assertEqual(first['observer_metric']['observer'], second['observer_metric']['observer'])
        self.assertEqual(first['metric_jet'].value, diagonal(2, 1))
        self.assertEqual(second['metric_jet'].value, diagonal(Q(3, 2), 1))
        # The change is anisotropic; one overall units conversion cannot fix it.
        self.assertNotEqual(first['metric_jet'].value[0][0], second['metric_jet'].value[0][0])

    def test_positive_record_cost_derives_only_the_nonnegative_quadratic_warp_branch(self):
        for c in (Q(-2), Q(0), Q(1, 3)):
            s = o.seam_foundation(c, 0, Q(2, 3))
            self.assertGreaterEqual(s['kappa'], 0)
            self.assertTrue(o.positive_definite(s['metric_jet'].value))
        with self.assertRaises(ValueError):
            o.seam_experiment(1, 0, Q(-1, 2))

    def test_response_rank_changes_are_reported_instead_of_assuming_a_global_quotient(self):
        results = []
        for v in (Q(0), Q(1)):
            t = n.matrix([[1, v], [0, 1]])
            h = diagonal(1, 1+2*v*v)
            e = o.EventExperiment(((1, 0),), ((1,),), (t, n.inverse(t)),
                                  (Q(1, 2),)*2, Q(1, 2), h, Q(3, 4))
            results.append(o.derive_observer_metric(e)['rank'])
        self.assertEqual(results, [1, 2])
