"""R6 error propagation and native prediction-corner checks."""

from fractions import Fraction as Q
from itertools import product
import os
import unittest

from interval_profile_recovery import (
    Interval as I, predict_interval, prediction_corners, recover_intervals,
)
from profile_recovery import finite_jet, native_queries, recover_prefix


def observations(last=216, error=Q(1, 100)):
    return tuple(I(c-error, c+error) for c in (3, -18, last))


def requests():
    queries = []
    for m in (2, 3):
        cells = recover_intervals(observations()[:m])['cell_intervals']
        for probe in (Q(1, 10), Q(1, 2)):
            for name, corner in prediction_corners(cells, probe, 3).items():
                queries.append({'id': f'corner_{m}_{probe}_{name}', 'kind': 'finite',
                                'cells': list(map(str, corner['cells'])),
                                'z': str(probe), 'tail': str(corner['tail'])})
    for name, prefix, period in (('base', [], [3, 2]), ('changed', [3, 2], [2, 1])):
        queries.append({'id': name, 'kind': 'periodic_prefix', 'prefix': prefix,
                        'pattern': period, 'z': '1/10', 'tail_cells': 48})
    return queries


class IntervalProfileRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = native_queries(os.environ['RKF_SOURCE_ROOT'], requests())
        cls.native = {item['id']: item['result'] for item in cls.packet['results']}

    def test_point_intervals_reproduce_exact_r5_recovery(self):
        for cells in ((3, 2, 3, 2), (Q(3, 2), Q(5, 3), Q(7, 4), Q(9, 5))):
            coefficients = finite_jet(cells, len(cells))
            result = recover_intervals([(c, c) for c in coefficients])
            self.assertEqual(result['status'], 'CERTIFIED_PREFIX')
            self.assertEqual(result['cell_intervals'], tuple(I(c, c) for c in cells))

    def test_noise_box_encloses_exact_recovery_at_all_coefficient_corners(self):
        box = observations()
        intervals = recover_intervals(box)['cell_intervals']
        for corner in product(*[(b.lower, b.upper) for b in box]):
            recovered = recover_prefix(corner)
            self.assertTrue(all(b.contains(a) for b, a in zip(intervals, recovered)))
        self.assertTrue(all(b.contains(a) for b, a in zip(intervals, (3, 2, 3))))

    def test_first_two_coupling_bounds_match_the_exact_monotone_formula(self):
        c0, c1 = observations()[:2]
        result = recover_intervals((c0, c1))
        self.assertEqual(result['cell_intervals'], (
            c0, I(-c1.upper/c0.upper**2, -c1.lower/c0.lower**2)))

    def test_third_depth_stays_distinguishable_with_supplied_error_bounds(self):
        base = recover_intervals(observations(216))['cell_intervals'][2]
        changed = recover_intervals(observations(180))['cell_intervals'][2]
        self.assertLess(changed.upper, base.lower)
        self.assertLessEqual(Q(29, 10), base.lower)
        self.assertGreaterEqual(Q(31, 10), base.upper)
        self.assertLessEqual(Q(19, 10), changed.lower)
        self.assertGreaterEqual(Q(21, 10), changed.upper)

    def test_crossing_zero_returns_a_partial_prefix_without_guessing(self):
        result = recover_intervals([(3, 3), (-18, -18), (107, 109)])
        self.assertEqual(result['status'], 'UNRESOLVED_POSITIVITY')
        self.assertEqual(result['stopped_depth'], 2)
        self.assertEqual(result['cell_intervals'], (I(3, 3), I(2, 2)))
        self.assertEqual(result['attempted_cell_interval'], I(Q(-1, 36), Q(1, 36)))
        self.assertTrue(result['unresolved_tail'])

    def test_an_entirely_nonpositive_depth_certifies_incompatibility(self):
        result = recover_intervals([(3, 3), (-18, -18), (100, 107)])
        self.assertEqual(result['status'], 'INCOMPATIBLE_POSITIVE_PROFILE')
        self.assertEqual(result['stopped_depth'], 2)
        self.assertLess(result['attempted_cell_interval'].upper, 0)
        zero = recover_intervals([(3, 3), (-18, -18), (108, 108)])
        self.assertEqual(zero['status'], 'INCOMPATIBLE_POSITIVE_PROFILE')

    def test_tighter_data_intervals_narrow_each_certified_coupling_range(self):
        wide = recover_intervals(observations())['cell_intervals']
        narrow = recover_intervals(observations(error=Q(1, 1000)))['cell_intervals']
        for w, n in zip(wide, narrow):
            self.assertLessEqual(w.lower, n.lower)
            self.assertGreaterEqual(w.upper, n.upper)
            self.assertLess(n.width, w.width)

    def test_even_and_odd_prefix_prediction_bounds_equal_native_corner_values(self):
        for m in (2, 3):
            cells = recover_intervals(observations()[:m])['cell_intervals']
            for probe in (Q(1, 10), Q(1, 2)):
                prediction = predict_interval(cells, probe, 3)
                self.assertEqual(prediction.lower, Q(self.native[f'corner_{m}_{probe}_lower']['value']))
                self.assertEqual(prediction.upper, Q(self.native[f'corner_{m}_{probe}_upper']['value']))

    def test_predictions_retain_completed_native_responses_and_hidden_tail(self):
        for name, last in (('base', 216), ('changed', 180)):
            cells = recover_intervals(observations(last))['cell_intervals']
            prediction = predict_interval(cells, '1/10', 3)
            native = self.native[name]['interval']
            self.assertTrue(prediction.contains(native['lower']))
            self.assertTrue(prediction.contains(native['upper']))
            exact_prediction = predict_interval([(a, a) for a in recover_prefix([3, -18, last])], '1/10', 3)
            self.assertLessEqual(prediction.lower, exact_prediction.lower)
            self.assertGreaterEqual(prediction.upper, exact_prediction.upper)
            self.assertGreater(exact_prediction.width, 0)

    def test_invalid_or_unresolved_arithmetic_is_rejected(self):
        with self.assertRaises(ValueError):
            I(2, 1)
        with self.assertRaises(TypeError):
            I(1.0, 2)
        with self.assertRaises(ZeroDivisionError):
            I(-1, 1).reciprocal()
        with self.assertRaises(ValueError):
            recover_intervals([])
        for cells, z, q in (([(0, 1)], 1, 3), ([(1, 2)], 0, 3), ([(1, 2)], 1, 0)):
            with self.assertRaises(ValueError):
                predict_interval(cells, z, q)


if __name__ == '__main__':
    unittest.main()
