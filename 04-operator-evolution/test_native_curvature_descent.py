"""R10 exact native conservation, descent, geometry and memory checks."""

import itertools
import unittest
import native_curvature_descent as n

Q, K, R, RK = n.Q, n.K, n.R, n.RK
I, Z = n.identity(2), n.scale(n.identity(2), 0)
SOURCES = None


def e3(i, j):
    return n.matrix([[int(a == i and b == j) for b in range(3)] for a in range(3)])


class NativeCurvatureDescentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if SOURCES is None:
            raise RuntimeError('Run verify_r10.py to pin-check source bytes before execution')
        cls.geometry = SOURCES['geometry']

    def test_visible_bianchi_has_exact_hidden_source(self):
        cut = n.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
        report = n.sourced_bianchi((e3(0, 1), e3(1, 2), e3(2, 0)), cut)
        current = n.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]])
        self.assertEqual(report['visible_bianchi'], current)
        self.assertEqual(report['hidden_source'], n.scale(current, -1))
        self.assertTrue(n.is_zero(report['full_bianchi']))
        self.assertTrue(n.is_zero(report['sourced_residual']))

    def test_no_odd_sector_has_no_hidden_bianchi_source(self):
        report = n.sourced_bianchi((e3(0, 1), e3(1, 2), e3(2, 0)), n.identity(3))
        self.assertTrue(n.is_zero(report['visible_bianchi']))
        self.assertTrue(n.is_zero(report['hidden_source']))

    def test_native_sourced_bianchi_on_distinct_cut_ranks(self):
        cuts = (n.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),
                n.matrix([[1, 0, 0], [0, -1, 0], [0, 0, -1]]))
        triples = ((e3(0, 1), e3(1, 2), e3(2, 0)),
                   (n.add(e3(0, 1), e3(2, 1)), e3(1, 2), n.add(e3(2, 0), e3(0, 0))),
                   (e3(0, 0), e3(1, 1), e3(2, 2)))
        for cut, operators in itertools.product(cuts, triples):
            report = n.sourced_bianchi(operators, cut)
            self.assertTrue(n.is_zero(report['full_bianchi']))
            self.assertTrue(n.is_zero(report['sourced_residual']))

    def test_smooth_sourced_bianchi_includes_curvature_derivatives(self):
        operators = (e3(0, 1), e3(1, 2), e3(2, 0))
        partials = tuple(tuple(n.scale(e3((i+j)%3, (i+2*j+1)%3), i+j+1)
                               for j in range(3)) for i in range(3))
        second = tuple(tuple(tuple(n.scale(e3((i+j+k)%3, (2*i+j+k+1)%3), (i+1)*(j+k+1))
                                    for k in range(3)) for j in range(3)) for i in range(3))
        jet = n.ConnectionJet(('x', 'y', 'z'), operators, partials)
        cut = n.matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
        report = n.sourced_bianchi_jet(jet, second, cut)
        self.assertTrue(any(not n.is_zero(x) for x in report['cyclic_curvature_derivatives']))
        self.assertFalse(n.is_zero(report['visible_bianchi']))
        self.assertTrue(n.is_zero(report['full_bianchi']))
        self.assertEqual(report['visible_bianchi'], n.scale(report['hidden_source'], -1))
        self.assertTrue(n.is_zero(report['sourced_residual']))

    def test_levi_civita_values_match_unchanged_geometry_source(self):
        for kappa, v in itertools.product((Q(3), Q(1, 5), Q(0), Q(-2), Q(7, 4)),
                                          (Q(0), Q(1, 2), Q(2), Q(-3, 4), Q(1, 3))):
            if 1+kappa*v*v <= 0:
                continue
            lc = n.warp_metric(kappa, v).levi_civita()
            a = self.geometry.christoffel_u_uv(v, kappa)
            b = self.geometry.christoffel_v_uu(v, kappa)
            self.assertEqual(lc.operators, (n.matrix([[0, a], [b, 0]]), n.matrix([[a, 0], [0, 0]])))

    def test_riemann_two_jet_matches_source_curvature_by_three_routes(self):
        for kappa, v in itertools.product((Q(3), Q(1, 5), Q(0), Q(-2), Q(7, 4)),
                                          (Q(0), Q(1, 2), Q(2), Q(-3, 4), Q(1, 3))):
            w = 1+kappa*v*v
            if w <= 0:
                continue
            metric = n.warp_metric(kappa, v)
            curvature = metric.levi_civita().curvature(0, 1)
            lowered = n.mul(metric.value, curvature)[0][1]
            self.assertEqual(lowered, self.geometry.riemann_uvuv(v, kappa))
            self.assertEqual(lowered/w, self.geometry.gaussian_curvature_closed(v, kappa))
            self.assertEqual(lowered/w, self.geometry.gaussian_curvature(v, kappa))

    def test_nondiagonal_metric_jet_is_metric_compatible_and_torsion_free(self):
        a, b, c = n.matrix([[1, 0], [0, 2]]), n.matrix([[0, 1], [1, 1]]), n.matrix([[1, 2], [2, 0]])
        metric = n.MetricJet(('x', 'y'), n.matrix([[2, 1], [1, 3]]), (a, b), ((a, c), (c, b)))
        lc = metric.levi_civita()
        report = n.tangent_compatibility(lc, metric)
        self.assertTrue(all(n.is_zero(x) for x in report['nonmetricity']))
        self.assertTrue(all(x == 0 for row in report['torsion'] for vector in row for x in vector))
        lowered = n.mul(metric.value, lc.curvature(0, 1))
        self.assertEqual(n.add(lowered, n.transpose(lowered)), Z)

    def test_flat_classical_adapters_in_one_two_and_three_dimensions(self):
        for dimension in (1, 2, 3):
            metric = n.flat_metric(dimension)
            self.assertTrue(all(n.is_zero(x) for x in metric.levi_civita().components().values()))
            self.assertEqual(len(metric.levi_civita().components()), dimension*(dimension-1)//2)

    def test_nondegenerate_indefinite_metric_is_an_admissible_tangent_adapter(self):
        metric = n.MetricJet(('x', 'y'), n.matrix([[-1, 0], [0, 1]]), (Z, Z), ((Z, Z), (Z, Z)))
        self.assertEqual(metric.levi_civita().curvature(0, 1), Z)

    def test_riemann_curvature_is_an_exact_descended_sector(self):
        metric = n.warp_metric(Q(3), Q(1, 3))
        hidden = n.ConnectionJet.constant(metric.directions, (K, R))
        native = n.extend_connection(metric.levi_civita(), hidden, lower=(K, R))
        report = n.riemann_reduction(native, metric, n.projection(2, 4))
        self.assertTrue(report['Riemann_sector_certified'])
        self.assertTrue(all(n.is_zero(x) for x in report['curvature_intertwining']))

    def test_same_classical_geometry_has_different_native_curvatures(self):
        metric = n.warp_metric(Q(3), Q(1, 3))
        curved = n.extend_connection(metric.levi_civita(), n.ConnectionJet.constant(metric.directions, (K, R)))
        quiet = n.extend_connection(metric.levi_civita(), n.ConnectionJet.constant(metric.directions, (Z, Z)))
        c = n.projection(2, 4)
        self.assertNotEqual(curved.curvature(0, 1), quiet.curvature(0, 1))
        self.assertEqual(n.rproduct(c, curved.curvature(0, 1)), n.rproduct(c, quiet.curvature(0, 1)))
        self.assertTrue(n.riemann_reduction(curved, metric, c)['Riemann_sector_certified'])
        self.assertTrue(n.riemann_reduction(quiet, metric, c)['Riemann_sector_certified'])

    def test_flat_riemann_sector_retains_nonzero_native_curvature(self):
        metric = n.flat_metric()
        native = n.extend_connection(metric.levi_civita(), n.ConnectionJet.constant(metric.directions, (K, R)))
        expected = n.block(Z, Z, Z, n.scale(RK, -2))
        self.assertEqual(native.curvature(0, 1), expected)
        self.assertTrue(n.riemann_reduction(native, metric, n.projection(2, 4))['Riemann_sector_certified'])

    def test_visible_to_hidden_memory_does_not_break_exact_descent(self):
        metric = n.flat_metric()
        native = n.extend_connection(metric.levi_civita(), n.ConnectionJet.constant(metric.directions, (K, R)), lower=(K, R))
        lower = n.matrix([row[:2] for row in native.curvature(0, 1)[2:]])
        self.assertEqual(lower, n.scale(RK, -2))
        self.assertTrue(n.riemann_reduction(native, metric, n.projection(2, 4))['Riemann_sector_certified'])

    def test_hidden_feedback_is_an_exact_excursion_curvature(self):
        metric = n.flat_metric()
        native = n.extend_connection(metric.levi_civita(), n.ConnectionJet.constant(metric.directions, (K, R)),
                                     upper=(I, Z), lower=(Z, K))
        report = n.visible_decomposition(native, metric)
        self.assertEqual(report['Riemann'], Z)
        self.assertEqual(report['distortion'], Z)
        self.assertEqual(report['excursion'], K)
        self.assertEqual(report['visible_full_curvature'], K)
        self.assertEqual(report['reconstruction_residual'], Z)
        self.assertFalse(n.riemann_reduction(native, metric, n.projection(2, 4))['Riemann_sector_certified'])

    def test_nonmetric_native_connection_has_a_distortion_curvature(self):
        metric = n.flat_metric()
        native = n.ConnectionJet.constant(metric.directions, (K, R))
        report = n.visible_decomposition(native, metric)
        self.assertEqual(report['Riemann'], Z)
        self.assertEqual(report['distortion'], n.scale(RK, -2))
        self.assertEqual(report['reconstruction_residual'], Z)
        self.assertFalse(n.tangent_compatibility(native, metric)['Levi_Civita_connection_jet'])

    def test_pointwise_metric_and_torsion_conditions_do_not_certify_curvature(self):
        metric = n.flat_metric()
        native = n.ConnectionJet(metric.directions, (Z, Z), ((Z, K), (Z, Z)))
        report = n.tangent_compatibility(native, metric)
        self.assertTrue(all(n.is_zero(x) for x in report['nonmetricity']))
        self.assertTrue(all(x == 0 for row in report['torsion'] for vector in row for x in vector))
        self.assertEqual(native.curvature(0, 1), n.scale(K, -1))
        self.assertFalse(report['Levi_Civita_connection_jet'])
        self.assertFalse(n.riemann_reduction(native, metric, I)['Riemann_sector_certified'])

    def test_accidental_curvature_equality_is_not_connection_identification(self):
        from emk_curvature_observation import derivative_cancellation_jet
        native = derivative_cancellation_jet(Q(2, 3))
        metric = n.flat_metric(directions=native.directions)
        report = n.riemann_reduction(native, metric, I)
        self.assertTrue(report['curvature_descends'])
        self.assertFalse(report['connection_jet_descends'])
        self.assertFalse(report['Riemann_sector_certified'])

    def test_nonsurjective_tangent_observer_is_rejected(self):
        metric = n.flat_metric()
        native = n.extend_connection(metric.levi_civita(), n.ConnectionJet.constant(metric.directions, (K, R)))
        with self.assertRaises(ValueError):
            n.riemann_reduction(native, metric, n.rzero(2, 4))

    def test_descent_covariance_under_constant_native_and_visible_reframes(self):
        metric = n.warp_metric(Q(3), Q(1, 3))
        lc = metric.levi_civita()
        native = n.extend_connection(lc, n.ConnectionJet.constant(metric.directions, (K, R)), lower=(K, R))
        s = n.matrix([[1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        t = n.matrix([[1, 2], [0, 1]])
        observer = n.rproduct(n.rproduct(t, n.projection(2, 4)), n.inverse(s))
        report = n.descent_report(native.change_reference(s), lc.change_reference(t), observer)
        self.assertTrue(report['connection_jet_descends'])
        self.assertTrue(report['curvature_descends'])

    def test_differentiated_observer_intertwining_is_load_bearing(self):
        metric = n.flat_metric()
        z4 = n.scale(n.identity(4), 0)
        derivative = n.block(Z, I, Z, Z)
        native = n.ConnectionJet(metric.directions, (z4, z4), ((z4, derivative), (z4, z4)))
        report = n.riemann_reduction(native, metric, n.projection(2, 4))
        self.assertTrue(all(n.is_zero(x) for x in report['value_intertwining']))
        self.assertFalse(report['connection_jet_descends'])
        self.assertFalse(report['curvature_descends'])

    def test_finite_visible_return_retains_hidden_order_holonomy(self):
        loop = n.lifted_order_loop()
        self.assertEqual(dict(loop.operators)['V'], n.block(I, Z, Z, n.scale(I, -1)))
        self.assertEqual(n.rproduct(n.projection(2, 4), dict(loop.operators)['V']), n.projection(2, 4))

    def test_independent_sheet_memory_survives_identity_transport(self):
        from emk_tensor_calculus import Slot, return_report
        loop = n.Move('p', 'p', (('V', n.identity(4)),), 2)
        signature = ((Slot('V', 1, 4),),)
        report = return_report(loop, signature)
        self.assertTrue(report['full_carrier_return'])
        self.assertEqual(report['sheet_residue'], 2)
        self.assertFalse(report['full_unledgered_return'])
        self.assertTrue(return_report(loop, signature, sheet_ledger=2)['full_return_after_sheet_ledger'])

    def test_finite_path_composition_descends_with_retained_hidden_memory(self):
        v1, v2 = n.add(I, n.scale(K, Q(1, 3))), n.add(I, n.scale(R, Q(1, 4)))
        t1, t2 = n.block(v1, Z, K, R), n.block(v2, Z, R, K)
        observer = n.projection(2, 4)
        self.assertEqual(n.rproduct(observer, n.mul(t2, t1)), n.rproduct(n.mul(v2, v1), observer))

    def test_all_three_visible_curvature_terms_reconstruct_full_result(self):
        metric = n.warp_metric(Q(3), Q(1, 3))
        lc = metric.levi_civita()
        visible = n.ConnectionJet(metric.directions, (n.add(lc.operators[0], K), n.add(lc.operators[1], R)),
                                  ((lc.partials[0][0], n.add(lc.partials[0][1], RK)), lc.partials[1]))
        native = n.extend_connection(visible, n.ConnectionJet.constant(metric.directions, (K, R)),
                                     upper=(I, K), lower=(R, RK))
        report = n.visible_decomposition(native, metric)
        self.assertTrue(all(not n.is_zero(report[key]) for key in ('Riemann', 'distortion', 'excursion')))
        self.assertEqual(report['reconstruction_residual'], Z)

    def test_invalid_geometry_and_transport_contracts_are_rejected(self):
        with self.assertRaises(ValueError):
            n.warp_metric(Q(-4), Q(1, 2))
        with self.assertRaises(TypeError):
            n.warp_metric(0.5, Q(1))
        with self.assertRaises(ValueError):
            n.MetricJet(('x', 'y'), n.matrix([[1, 2], [0, 1]]), (Z, Z), ((Z, Z), (Z, Z)))
        with self.assertRaises(ValueError):
            n.MetricJet(('x', 'y'), I, (Z, Z), ((Z, K), (Z, Z)))
        with self.assertRaises(ValueError):
            n.sourced_bianchi((K, R), K)
        with self.assertRaises(ValueError):
            n.descent_report(n.ConnectionJet.constant(('a', 'b'), (K, R)), n.flat_metric().levi_civita(), I)
