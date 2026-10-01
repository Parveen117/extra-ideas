"""Focused R11 tests: full curvature, blind spectra, active markers and limits."""

import itertools
from math import comb
import unittest
import spectral_curvature_observer as s

n, Q = s.n, s.Q
SOURCE_PINS_VERIFIED = False
LOWER = ((Q(2), Q(-3, 4)), (Q(5, 7), Q(1)))
SAMPLES = (Q(-2), Q(-1, 3), Q(0), Q(1, 2), Q(3))


def lower_grid():
    for entries in itertools.product((-1, 0, 1), repeat=4):
        yield (entries[:2], entries[2:])


class SpectralCurvatureObserverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not SOURCE_PINS_VERIFIED:
            raise RuntimeError('Run verify_r11.py to check pinned spectral PDFs and native source bytes')

    def test_curvature_is_nonzero_nilpotent_and_has_a_flat_riemann_quotient(self):
        for h in (1, 2, 3):
            lower = tuple((Q(i+1), Q(-i-2)) for i in range(h))
            actual = s.native_connection(lower).curvature(0, 1)
            self.assertEqual(actual, s.nilpotent_carrier(lower))
            self.assertFalse(n.is_zero(actual))
            self.assertTrue(n.is_zero(n.mul(actual, actual)))
            self.assertTrue(s.classical_report(lower)['Riemann_sector_certified'])
            self.assertTrue(all(n.is_zero(x) for x in s.classical_report(lower)['tangent_curvature'].values()))

    def test_all_curvature_characteristic_coefficients_forget_the_whole_family(self):
        for lower in lower_grid():
            self.assertEqual(s.spectral_coefficients(s.nilpotent_carrier(lower)), (1, 0, 0, 0, 0))

    def test_independent_determinants_replay_unmarked_blindness_at_all_samples(self):
        for lower, q in itertools.product(lower_grid(), SAMPLES):
            self.assertEqual(s.spectral_value(s.nilpotent_carrier(lower), q), 1)

    def test_all_diagonal_entries_and_trace_powers_can_vanish_for_curved_transport(self):
        operator = s.nilpotent_carrier(LOWER)
        self.assertTrue(all(operator[i][i] == 0 for i in range(4)))
        power = n.identity(4)
        for _ in range(1, 9):
            power = n.mul(power, operator)
            self.assertEqual(s.trace(power), 0)
        self.assertFalse(n.is_zero(s.native_connection(LOWER).curvature(0, 1)))

    def test_two_ordered_factors_have_equal_spectral_shadows_and_nonzero_residue(self):
        x, y = s.order_factors(LOWER, Q(1, 3), Q(2, 5))
        xy, yx = n.mul(x, y), n.mul(y, x)
        self.assertNotEqual(xy, yx)
        self.assertEqual(n.sub(xy, yx), n.scale(s.nilpotent_carrier(LOWER), Q(2, 15)))
        self.assertEqual(s.spectral_coefficients(xy), s.spectral_coefficients(yx))
        for q in SAMPLES:
            self.assertEqual(s.spectral_value(xy, q), s.spectral_value(yx, q))

    def test_exact_finite_loop_retains_curvature_including_both_orientation_signs(self):
        for a, b in itertools.product((Q(-2), Q(-1, 2), Q(0), Q(1, 3), Q(1)),
                                      (Q(-3, 4), Q(0), Q(1, 2))):
            loop = s.order_loop(LOWER, a, b)
            self.assertEqual(dict(loop.operators)['V'],
                             n.add(n.identity(4), n.scale(s.nilpotent_carrier(LOWER), a*b)))
            self.assertEqual(loop.sheet, 0)

    def test_nonidentity_order_loop_is_spectrally_identical_to_identity(self):
        returned = dict(s.order_loop(LOWER).operators)['V']
        self.assertNotEqual(returned, n.identity(4))
        self.assertEqual(s.spectral_coefficients(returned), tuple((-1)**k*comb(4, k) for k in range(5)))
        for q in SAMPLES:
            self.assertEqual(s.spectral_value(returned, q), (1-q)**4)
        self.assertEqual(n.rproduct(n.projection(2, 4), returned), n.projection(2, 4))

    def test_source_insertion_markers_recover_every_lower_coupling_entry(self):
        for lower in lower_grid():
            operator = s.nilpotent_carrier(lower)
            self.assertEqual(tuple(s.marked_trace(marker, operator) for marker in s.marker_bank(2)),
                             tuple(Q(x) for row in lower for x in row))

    def test_source_marked_determinants_have_the_exact_linear_coefficient(self):
        for lower, q in itertools.product(lower_grid(), SAMPLES):
            operator = s.nilpotent_carrier(lower)
            for marker, value in zip(s.marker_bank(2), (x for row in lower for x in row)):
                self.assertEqual(s.marked_determinant(marker, operator, q), 1-q*value)

    def test_complete_marked_loop_bank_recovers_the_original_curvature(self):
        a, b, q = Q(-1, 2), Q(3, 5), Q(7, 3)
        returned = dict(s.order_loop(LOWER, a, b).operators)['V']
        values = tuple(s.marked_determinant(marker, returned, q) for marker in s.marker_bank(2))
        self.assertEqual(s.decode_determinants(values, 2, q, a*b), s.lower_data(LOWER))

    def test_hermitian_markers_are_complete_on_this_real_curvature_family(self):
        operator = s.nilpotent_carrier(LOWER)
        for marker, expected in zip(s.symmetric_marker_bank(2), (x for row in LOWER for x in row)):
            self.assertEqual(n.transpose(marker), marker)
            self.assertEqual(s.marked_trace(marker, operator), expected)
            for q in SAMPLES:
                self.assertEqual(s.marked_determinant(marker, operator, q), 1-q*expected)

    def test_minimum_is_two_scalar_probes_per_hidden_mode(self):
        for h in (1, 2, 3):
            columns = 2*h
            zero_rows = ((Q(0),)*columns,)
            report = s.observer_report(zero_rows, n.identity(columns), columns)
            self.assertEqual(report['minimum_arbitrary_scalar_supplements'], columns)
            self.assertEqual(report['blind_dimension'], columns)
            self.assertEqual(s.analysis_rows(s.marker_bank(h), h), n.identity(columns))

    def test_every_subbank_has_its_exact_remaining_curvature_blindness(self):
        markers = s.marker_bank(2)
        for bits in itertools.product((0, 1), repeat=4):
            chosen = tuple(marker for marker, bit in zip(markers, bits) if bit)
            report = s.observer_report(s.analysis_rows(chosen, 2), n.identity(4), 4)
            self.assertEqual(report['target_blindness'], 4-sum(bits))

    def test_each_missing_cross_sector_marker_has_a_nonzero_blind_witness(self):
        markers = s.marker_bank(2)
        for missing in range(4):
            lower = tuple(tuple(Q(int(2*r+c == missing)) for c in range(2)) for r in range(2))
            operator = s.nilpotent_carrier(lower)
            self.assertFalse(n.is_zero(operator))
            self.assertTrue(all(s.marked_trace(marker, operator) == 0 for i, marker in enumerate(markers) if i != missing))
            self.assertEqual(s.marked_trace(markers[missing], operator), 1)

    def test_nonfeedback_catalogue_cannot_repair_curvature_even_with_spectral_sampling(self):
        operator = s.nilpotent_carrier(LOWER)
        catalogue = s.nonfeedback_catalogue(2)
        self.assertEqual(len(catalogue), 12)
        report = s.observer_report(s.analysis_rows(catalogue, 2), n.identity(4), 4)
        self.assertEqual(report['target_blindness'], 4)
        self.assertIsNone(s.catalogue_minimum(catalogue, 2, n.identity(4))['minimum_catalogue_channels'])
        for marker, q in itertools.product(catalogue, SAMPLES):
            self.assertEqual(s.marked_determinant(marker, operator, q), 1)

    def test_marker_and_transport_covariance_retains_the_same_responses(self):
        operator = s.nilpotent_carrier(LOWER)
        change = n.matrix([[1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 1, 2], [0, 0, 0, 1]])
        inverse = n.inverse(change)
        transformed = n.mul(n.mul(change, operator), inverse)
        observer = n.rproduct(n.projection(2, 4), inverse)
        native = s.native_connection(LOWER).change_reference(change)
        self.assertTrue(n.riemann_reduction(native, n.flat_metric(directions=('u', 'v')), observer)['Riemann_sector_certified'])
        for marker in s.marker_bank(2):
            moved_marker = n.mul(n.mul(change, marker), inverse)
            self.assertEqual(s.marked_trace(moved_marker, transformed), s.marked_trace(marker, operator))
            for q in SAMPLES:
                self.assertEqual(s.marked_determinant(moved_marker, transformed, q),
                                 s.marked_determinant(marker, operator, q))

    def test_active_cut_reverses_signed_marker_response_but_preserves_spectrum(self):
        operator, j = s.nilpotent_carrier(LOWER), s.cut(2)
        changed = n.mul(n.mul(j, operator), j)
        self.assertEqual(changed, n.scale(operator, -1))
        self.assertEqual(s.spectral_coefficients(changed), s.spectral_coefficients(operator))
        for marker in s.marker_bank(2):
            self.assertEqual(s.marked_trace(marker, changed), -s.marked_trace(marker, operator))

    def test_balance_erases_this_four_dimensional_readout_without_flattening_raw_transport(self):
        operator = s.native_connection(LOWER).curvature(0, 1)
        recognized = s.recognize(operator, Q(1, 2), s.cut(2))
        self.assertTrue(n.is_zero(recognized))
        self.assertFalse(n.is_zero(operator))
        self.assertEqual(s.balance_response_report(2, Q(1, 2))['observer']['target_blindness'], 4)

    def test_recorded_cut_branch_recovery_is_exact(self):
        operator, j = s.nilpotent_carrier(LOWER), s.cut(2)
        for branch in (0, 1):
            realized = n.mul(n.mul(j, operator), j) if branch else operator
            restored = n.mul(n.mul(j, realized), j) if branch else realized
            self.assertEqual(restored, operator)

    def test_response_metric_and_target_blindness_track_conditional_balance(self):
        for theta in (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)):
            report = s.balance_response_report(2, theta)
            multiplier = 1-2*theta
            self.assertEqual(report['response_Gram'], n.scale(n.identity(4), multiplier**2))
            self.assertEqual(report['observer']['target_blindness'], 4 if multiplier == 0 else 0)
            self.assertEqual(report['odd_multiplier'], multiplier)

    def test_cubic_excitation_replays_the_papers_connected_order_formula(self):
        operator, projector = s.nilpotent_carrier(LOWER), s.hidden_projector(2)
        for marker, q in itertools.product(s.marker_bank(2), (Q(-2), Q(0), Q(1, 3), Q(2))):
            computed = s.cubic_connected_difference(projector, operator, marker, q)
            expected = q/(1-q)**2*s.marked_trace(marker, operator)
            self.assertEqual(computed, expected)

    def test_cubic_channel_still_requires_a_marker_that_closes_the_hidden_path(self):
        for marker in s.nonfeedback_catalogue(2):
            self.assertEqual(s.cubic_connected_difference(s.hidden_projector(2), s.nilpotent_carrier(LOWER), marker, Q(2)), 0)

    def test_marker_closure_is_the_same_upper_lower_excursion_in_r10(self):
        for marker, expected in zip(s.marker_bank(2), (x for row in LOWER for x in row)):
            report = s.feedback_report(LOWER, marker)
            self.assertTrue(n.is_zero(report['Riemann']))
            self.assertTrue(n.is_zero(report['distortion']))
            self.assertTrue(n.is_zero(report['reconstruction_residual']))
            self.assertEqual(report['trace_visible_excursion'], expected)
            self.assertEqual(report['marked_curvature_response'], expected)
            self.assertEqual(report['visible_full_curvature'], report['excursion'])

    def test_deterministic_noise_bound_is_sharp_for_the_complete_bank(self):
        q, factor, error = Q(2), Q(-1, 3), Q(1, 100)
        operator = n.scale(s.nilpotent_carrier(LOWER), factor)
        true = tuple(s.marked_determinant(marker, operator, q) for marker in s.marker_bank(2))
        measured = tuple(value+error for value in true)
        report = s.classify_determinants(measured, 2, q, error, factor)
        actual = sum(((x-y)**2 for xr, yr in zip(report['lower_estimate'], LOWER) for x, y in zip(xr, yr)), Q(0))
        self.assertEqual(actual, report['squared_Frobenius_error_bound'])

    def test_noisy_zero_responses_do_not_certify_flatness(self):
        report = s.classify_determinants((1, 1, 1, 1), 2, Q(2), Q(1, 100))
        self.assertEqual(report['status'], 'UNRESOLVED_CURVATURE_AT_DECLARED_ERROR')
        self.assertTrue(all(lo < 0 < hi for lo, hi in report['entry_intervals']))

    def test_exact_flatness_and_robust_nonzero_curvature_have_distinct_certificates(self):
        self.assertEqual(s.classify_determinants((1, 1, 1, 1), 2, 1)['status'],
                         'CERTIFIED_ZERO_CURVATURE_IN_DECLARED_FAMILY')
        self.assertEqual(s.classify_determinants((Q(3, 4), 1, 1, 1), 2, 1, Q(1, 100))['status'],
                         'CERTIFIED_NONZERO_CURVATURE_IN_DECLARED_FAMILY')

    def test_the_marker_bank_is_model_relative_and_refuses_other_hidden_curvature(self):
        operator = n.block(n.scale(n.identity(2), 0), n.rzero(2, 2), n.rzero(2, 2), n.RK)
        self.assertFalse(n.is_zero(operator))
        self.assertTrue(all(s.marked_trace(marker, operator) == 0 for marker in s.marker_bank(2)))
        with self.assertRaises(ValueError):
            s.extract_lower(operator)

    def test_endpoint_markers_cannot_reconstruct_erased_histories_or_sheet_registers(self):
        operator, unit = s.nilpotent_carrier(LOWER), n.identity(4)
        endpoint = n.mul(n.add(unit, operator), n.sub(unit, operator))
        self.assertEqual(endpoint, unit)
        for marker, q in itertools.product(s.marker_bank(2), SAMPLES):
            self.assertEqual(s.marked_determinant(marker, endpoint, q), s.marked_determinant(marker, unit, q))
        path = n.Move('p', 'p', (('V', endpoint),), 2)
        from emk_tensor_calculus import Slot, return_report
        audit = return_report(path, ((Slot('V', 1, 4),),))
        self.assertFalse(audit['full_unledgered_return'])
        self.assertEqual(audit['sheet_residue'], 2)

    def test_target_specific_catalogue_overhead_is_not_full_curvature_probe_count(self):
        target = ((Q(1), Q(0), Q(0), Q(1)),)
        zero_rows = ((Q(0),)*4,)
        self.assertEqual(s.observer_report(zero_rows, target, 4)['minimum_arbitrary_scalar_supplements'], 1)
        self.assertEqual(s.catalogue_minimum(s.marker_bank(2), 2, target)['minimum_catalogue_channels'], 2)

    def test_repeated_erasure_is_counted_once_in_the_target_ledger(self):
        zero_rows = ((Q(0),)*4,)
        ledger = s.blindness_ledger((n.identity(4), zero_rows, zero_rows, zero_rows), n.identity(4), 4)
        self.assertEqual(ledger['target_blindness_by_stage'], (0, 4, 4, 4))
        self.assertEqual(ledger['increments'], (4, 0, 0))
        with self.assertRaises(ValueError):
            s.blindness_ledger((zero_rows, n.identity(4)), n.identity(4), 4)

    def test_degenerate_gain_and_invalid_declared_contracts_are_refused(self):
        with self.assertRaises(ValueError):
            s.order_loop(LOWER, -1, 1)
        with self.assertRaises(ValueError):
            s.decode_determinants((1, 1, 1, 1), 2, 0)
        with self.assertRaises(ValueError):
            s.decode_determinants((1, 1, 1, 1), 2, 1, 0)
        with self.assertRaises(ValueError):
            s.cubic_connected_difference(s.hidden_projector(2), s.nilpotent_carrier(LOWER), s.marker_bank(2)[0], 1)
        with self.assertRaises(ValueError):
            s.classify_determinants((1, 1, 1, 1), 2, 1, -1)
        with self.assertRaises(TypeError):
            s.nilpotent_carrier(((0.5, 1),))
        with self.assertRaises(ValueError):
            s.nilpotent_carrier(((1, 2, 3),))
