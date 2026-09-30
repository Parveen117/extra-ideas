"""Focused R5 checks, including the unchanged native forward solver."""

from fractions import Fraction as Q
import os
from pathlib import Path
import tempfile
import unittest

from profile_recovery import (
    finite_fraction, finite_jet, first_difference, native_queries,
    recover_prefix, weak_tail_bound,
)


PROFILES = ((3, 2, 3, 2), (3, 2, 2, 1),
            (Q(3, 2), Q(5, 3), Q(7, 4), Q(9, 5), Q(11, 6)))
PROBES = (Q(1, 10), Q(1, 2), Q(7, 4))
BOUND_CASES = (((3, 2), Q(1, 10)), ((3, 2, 3), Q(1, 10)),
               ((3, 2), Q(1, 2)))


def evaluate(polynomial, argument):
    value = Q(0)
    for coefficient in reversed(polynomial):
        value = coefficient + argument * value
    return value


def requests():
    queries = []
    for i, cells in enumerate(PROFILES):
        for j, probe in enumerate(PROBES):
            queries.append({'id': f'poly_{i}_{j}', 'kind': 'finite',
                            'cells': list(map(str, cells)), 'z': str(probe)})
    for m in range(1, 7):
        prefix = [3 if j % 2 == 0 else 2 for j in range(m)]
        target = '1' if m % 2 == 0 else '2/3'
        queries.append({'id': f'closed_prefix_{m}', 'kind': 'finite',
                        'cells': prefix, 'tail': target})
    for parity, target in (('even', '1'), ('odd', '2/3')):
        queries.append({'id': f'design_{parity}', 'kind': 'synthesize',
                        'desired': target, 'second': '1'})
    for name, prefix, period in (('base', [], [3, 2]), ('changed', [3, 2], [2, 1])):
        for label, probe in (('closed', '1'), ('heldout', '1/2')):
            queries.append({'id': f'{name}_{label}', 'kind': 'periodic_prefix',
                            'prefix': prefix, 'pattern': period, 'z': probe,
                            'tail_cells': 48})
    for i, (prefix, probe) in enumerate(BOUND_CASES):
        queries.extend([
            {'id': f'bound_{i}', 'kind': 'prefix_interval', 'prefix': prefix,
             'z': str(probe), 'lower': '0', 'upper': str(3*probe)},
            {'id': f'transfer_{i}', 'kind': 'finite', 'cells': prefix,
             'z': str(probe)},
        ])
    gain, dial, probe = Q(3, 2), Q(2, 3), Q(1, 4)
    cells = PROFILES[2]
    gauged = [dial * (gain if j % 2 == 0 else 1/gain) * cell
              for j, cell in enumerate(cells)]
    queries.extend([
        {'id': 'calibration_base', 'kind': 'finite', 'cells': list(map(str, cells)),
         'z': str(dial*probe)},
        {'id': 'calibration_gauged', 'kind': 'finite', 'cells': list(map(str, gauged)),
         'z': str(probe)},
        {'id': 'native_schur', 'kind': 'schur', 'cell': '3', 'z': '1/2',
         'tail': '2/3'},
    ])
    return queries


class ProfileRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = native_queries(os.environ['RKF_SOURCE_ROOT'], requests())
        cls.native = {item['id']: item['result'] for item in cls.packet['results']}

    def test_known_coefficients_recover_the_declared_depth(self):
        self.assertEqual(recover_prefix([3, -18, 216]), (3, 2, 3))
        self.assertEqual(finite_jet([3, 2, 3, 2], 4), (3, -18, 216, -3240))
        self.assertEqual(finite_jet([1]*6, 6), (1, -1, 2, -5, 14, -42))
        # These are coefficients; the third derivative itself would be -108.
        self.assertEqual(recover_prefix([3, -18]), (3, 2))

    def test_polynomial_companion_agrees_with_native_finite_inverses(self):
        for i, cells in enumerate(PROFILES):
            numerator, denominator = finite_fraction(cells)
            for j, probe in enumerate(PROBES):
                with self.subTest(cells=cells, probe=probe):
                    rational = probe * evaluate(numerator, probe**2) / evaluate(denominator, probe**2)
                    self.assertEqual(rational, Q(self.native[f'poly_{i}_{j}']['value']))

    def test_heterogeneous_prefixes_recover_without_a_period_assumption(self):
        cells = PROFILES[2] + (Q(13, 7), Q(17, 8))
        for m in range(1, len(cells)+1):
            with self.subTest(depth=m):
                self.assertEqual(recover_prefix(finite_jet(cells, m)), cells[:m])
                self.assertEqual(finite_jet(cells[:m], m), finite_jet(cells, m))

    def test_first_different_depth_has_the_proved_first_visible_order(self):
        for m in range(7):
            prefix = tuple(3 if j % 2 == 0 else 2 for j in range(m))
            left, right = Q(7, 3), Q(11, 5)
            predicted = first_difference(prefix, left, right)
            a = finite_jet(prefix + (left,), m+1)
            b = finite_jet(prefix + (right,), m+1)
            self.assertEqual(a[:m], b[:m])
            self.assertEqual(predicted, {'power': 2*m+1, 'coefficient': a[m]-b[m]})

    def test_alternating_signs_alone_do_not_certify_a_positive_profile(self):
        for coefficients in ([], [0], [1, 1], [1, -1, '1/2'], [1, -1, 1]):
            with self.subTest(coefficients=coefficients):
                with self.assertRaises(ValueError):
                    recover_prefix(coefficients)

    def test_exact_input_contract_rejects_floats_and_nonpositive_cells(self):
        for coefficients in ([3.0, -18], [True, -1]):
            with self.assertRaises(TypeError):
                recover_prefix(coefficients)
        for cells in ([], [1, 0], [1, -2]):
            with self.assertRaises(ValueError):
                finite_fraction(cells)
        with self.assertRaises(ValueError):
            finite_jet([1], 0)
        with self.assertRaises(ValueError):
            first_difference([1], 2, 2)

    def test_native_tail_design_keeps_exact_closure_at_arbitrary_prefix_depth(self):
        self.assertEqual(self.native['design_even']['period'], ['2', '1'])
        self.assertEqual(self.native['design_odd']['period'], ['10/9', '1'])
        for m in range(1, 7):
            self.assertEqual(Q(self.native[f'closed_prefix_{m}']['value']), 1)
        # Closure is proved by source synthesis plus the prefix, not by a
        # finite interval containing one. These finite intervals remain open.
        for name in ('base', 'changed'):
            interval = self.native[f'{name}_closed']['interval']
            self.assertLessEqual(Q(interval['lower']), 1)
            self.assertGreaterEqual(Q(interval['upper']), 1)
            self.assertGreater(Q(interval['width']), 0)

    def test_same_closed_bond_and_two_coefficients_hide_a_changed_tail(self):
        base = finite_jet([3, 2, 3], 3)
        changed = finite_jet([3, 2, 2], 3)
        self.assertEqual(base, (3, -18, 216))
        self.assertEqual(changed, (3, -18, 180))
        self.assertEqual(base[:2], changed[:2])
        self.assertEqual(base[2]-changed[2], 36)

    def test_heldout_native_probe_separates_the_completed_profiles(self):
        base = self.native['base_heldout']['interval']
        changed = self.native['changed_heldout']['interval']
        blo, bhi = Q(base['lower']), Q(base['upper'])
        clo, chi = Q(changed['lower']), Q(changed['upper'])
        self.assertLess(chi, blo)
        self.assertLess(bhi-blo, Q(1, 10**12))
        self.assertLess(chi-clo, Q(1, 10**12))
        # Independent algebra: base=(-1+sqrt(7))/2; changed=6-3*sqrt(3).
        self.assertLess(2*blo*blo+2*blo-3, 0)
        self.assertGreater(2*bhi*bhi+2*bhi-3, 0)
        # The second polynomial decreases throughout the certified interval.
        self.assertLess(chi, 1)
        self.assertGreater(clo*clo-12*clo+9, 0)
        self.assertLess(chi*chi-12*chi+9, 0)

    def test_tail_prior_yields_the_exact_transfer_width_and_weak_bound(self):
        for i, (prefix, probe) in enumerate(BOUND_CASES):
            interval = self.native[f'bound_{i}']
            transfer = self.native[f'transfer_{i}']['transfer']
            c, d = Q(transfer['C']), Q(transfer['D'])
            width = 3*probe / (d*(c*3*probe+d))
            self.assertEqual(width, Q(interval['width']))
            self.assertLessEqual(width, weak_tail_bound(prefix, probe, 3))
        with self.assertRaises(ValueError):
            weak_tail_bound([3, 2], '1/2', 0)

    def test_unknown_probe_and_readout_gains_change_recovered_products(self):
        gain, dial = Q(3, 2), Q(2, 3)
        base = Q(self.native['calibration_base']['value'])
        gauged = Q(self.native['calibration_gauged']['value'])
        self.assertEqual(gauged, gain*base)
        cells = PROFILES[2]
        coefficients = finite_jet(cells, len(cells))
        transformed = [gain*dial**(2*j+1)*value for j, value in enumerate(coefficients)]
        expected = tuple(dial*(gain if j % 2 == 0 else 1/gain)*a
                         for j, a in enumerate(cells))
        self.assertEqual(recover_prefix(transformed), expected)
        self.assertNotEqual(expected, cells)

    def test_aperture_shift_is_observable_in_calibrated_coefficients(self):
        self.assertEqual(finite_jet([2, 3, 2], 3), (2, -12, 144))
        self.assertEqual(recover_prefix([2, -12, 144]), (2, 3, 2))
        self.assertNotEqual(recover_prefix([2, -12]), recover_prefix([3, -18]))

    def test_native_schur_proofs_replay_for_the_scaled_cell(self):
        result = self.native['native_schur']
        self.assertEqual(result['response']['u'], '1')
        self.assertEqual(result['response']['v'], '3/4')
        self.assertEqual(set(result['proofs']), {'schur', 'leftInverse', 'rightInverse', 'response'})
        self.assertTrue(all(proof['replayed_by_source'] for proof in result['proofs'].values()))

    def test_changed_native_runtime_is_refused_before_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'operator_foundation/core/native_operator.cjs'
            path.parent.mkdir(parents=True)
            path.write_text('// changed source must never execute\n')
            with self.assertRaisesRegex(ValueError, 'Upstream source hash mismatch'):
                native_queries(directory, [])


if __name__ == '__main__':
    unittest.main()
