"""R4 source-integration checks; set RKF_SOURCE_ROOT to the pinned checkout."""

from fractions import Fraction as Q
import os
from pathlib import Path
import tempfile
import unittest

from aghora_return import (
    add, cayley_flow, conjugate, identity, inverse, matrix, mul,
    return_operator, scale, sub, transpose,
)
from native_bond_bridge import (
    NATIVE_K, NATIVE_L, NATIVE_R, check_runtime_pins, closure_budget,
    load_native_probe, period_cut_condition, response_bridge,
)
from seam_bond import seam_geometry


I = identity(2)
ZERO = scale(I, 0)


class NativeBondBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source_root = Path(os.environ['RKF_SOURCE_ROOT'])
        cls.packet = load_native_probe(cls.source_root)

    def test_native_proof_replays_and_actual_schur_response(self):
        self.assertEqual(len(self.packet['checked_runtime_sha256']), 4)
        self.assertEqual(len(self.packet['algebra_checks']), 4)
        for case in self.packet['algebra_checks']:
            self.assertEqual(len(case['proofs']), 4)
            for proof in case['proofs'].values():
                self.assertTrue(proof['replayed'])
                self.assertEqual(proof['normal_terms'], 0)
        response = self.packet['source_schur_response']
        self.assertEqual((Q(response['u']), Q(response['v'])), (Q(1), Q(1)))
        self.assertEqual(set(self.packet['source_schur_proof_names']),
                         {'schur', 'leftInverse', 'rightInverse', 'response'})

    def test_source_supplies_closed_bond_and_opposite_seam(self):
        for design in self.packet['designs'][:2]:
            a, b = map(Q, design['period'])
            self.assertTrue(period_cut_condition(a, b)['exact_completed_cut'])
            result = response_bridge(Q(design['desiredReturn']))
            bond = result['bond_candidate']
            self.assertEqual(bond, matrix(((1, 0), (0, 0))))
            self.assertEqual(mul(mul(NATIVE_K, bond), NATIVE_K), sub(I, bond))
            self.assertEqual(result['commutator'], NATIVE_R)
            self.assertEqual(result['bond_defect'], ZERO)
            self.assertEqual(result['complex_defect'], ZERO)

    def test_full_r3_chart_and_metric_are_intertwined(self):
        native_bond = response_bridge(1)['bond_candidate']
        for gamma in (Q(2), Q(-3, 5), Q(1, 7)):
            with self.subTest(gamma=gamma):
                r3 = seam_geometry(gamma)
                change = r3['basis']
                self.assertEqual(conjugate(NATIVE_K, change), r3['aghora'])
                self.assertEqual(conjugate(native_bond, change), r3['bond'])
                self.assertEqual(conjugate(NATIVE_R, change), r3['complex_structure'])
                pulled_metric = mul(transpose(inverse(change)), inverse(change))
                self.assertEqual(scale(pulled_metric, 2), r3['metric'])

    def test_raw_defects_are_linked_for_source_responses(self):
        values = [Q(d['desiredReturn']) for d in self.packet['designs']]
        values += [Q(p['actual_zero_tail']) for p in self.packet['prefixes'].values()]
        for x in values:
            result = response_bridge(x)
            self.assertEqual(result['commutator'], scale(NATIVE_R, x))
            self.assertEqual(result['complex_defect'], scale(result['bond_defect'], -4))
            self.assertEqual(result['bond_defect'], scale(I, (x*x-1)/4))

    def test_normalized_iota_and_exchange_can_hide_nonclosure(self):
        design = self.packet['designs'][2]
        self.assertEqual(tuple(map(Q, design['period'])), (Q(30), Q(20)))
        x = Q(design['desiredReturn'])
        self.assertEqual(x, Q(6, 5))
        result = response_bridge(x)
        bond = result['bond_candidate']
        self.assertEqual(mul(mul(NATIVE_K, bond), NATIVE_K), sub(I, bond))
        normalized = scale(result['commutator'], 1/x)
        self.assertEqual(mul(normalized, normalized), scale(I, -1))
        self.assertEqual(result['bond_defect'], scale(I, Q(11, 100)))
        self.assertFalse(result['exact_cut'])
        # Sharpening also manufactures a different idempotent from every x>0.
        sharpened = scale(add(I, scale(sub(scale(bond, 2), I), 1/x)), Q(1, 2))
        self.assertEqual(mul(sharpened, sharpened), sharpened)
        self.assertNotEqual(sharpened, bond)

    def test_finite_exact_cut_does_not_certify_completed_cut(self):
        prefixes = self.packet['prefixes']
        first, later = prefixes['accidental_finite_cut'], prefixes['reopened_same_profile']
        self.assertEqual(Q(first['actual_zero_tail']), 1)
        self.assertTrue(response_bridge(Q(first['actual_zero_tail']))['exact_cut'])
        lo, hi = Q(later['interval']['lower']), Q(later['interval']['upper'])
        self.assertEqual((lo, hi), (Q(1, 2), Q(2, 3)))
        self.assertTrue(closure_budget(lo, hi)['excludes_exact_cut'])
        self.assertFalse(period_cut_condition(1, 1)['exact_completed_cut'])

    def test_native_enclosures_control_both_defects(self):
        for name, target in (('cut_32_short', Q(1)), ('cut_32_long', Q(1)),
                             ('noncut_30_20', Q(6, 5))):
            packet = self.packet['prefixes'][name]
            lo, hi = Q(packet['interval']['lower']), Q(packet['interval']['upper'])
            budget = closure_budget(lo, hi)
            self.assertLessEqual(lo, target)
            self.assertLessEqual(target, hi)
            expected = (target*target-1)/4
            dl, du = budget['bond_defect_interval']
            self.assertLessEqual(dl, expected)
            self.assertLessEqual(expected, du)
            cl, cu = budget['complex_defect_interval']
            self.assertEqual((cl, cu), (-4*du, -4*dl))
            self.assertLessEqual(abs(target-1)/2, budget['distance_to_cut_bound'])

    def test_source_refinement_shrinks_the_cut_error_budget(self):
        packets = self.packet['prefixes']
        short, long = [packets[name]['interval'] for name in ('cut_32_short', 'cut_32_long')]
        sl, su, ll, lu = map(Q, (short['lower'], short['upper'], long['lower'], long['upper']))
        self.assertLessEqual(sl, ll)
        self.assertLessEqual(lu, su)
        self.assertLess(closure_budget(ll, lu)['distance_to_cut_bound'],
                        closure_budget(sl, su)['distance_to_cut_bound'])

    def test_cell_mismatch_has_the_same_sign_as_bond_defect(self):
        for design in self.packet['designs']:
            a, b = map(Q, design['period'])
            x = Q(design['desiredReturn'])
            delta = period_cut_condition(a, b)['mismatch']
            self.assertEqual(delta, (x-1)*(b*(x+1)+1))
            self.assertEqual((x*x-1)/4, delta*(x+1)/(4*(b*(x+1)+1)))
        # The original source family includes a non-cut phase below one too.
        a, b, x = Q(2), Q(3), Q(2, 3)
        self.assertEqual(b*x*x+x, a)
        self.assertLess(period_cut_condition(a, b)['mismatch'], 0)
        self.assertLess((x*x-1)/4, 0)

    def test_selecting_bond_does_not_select_return_phase(self):
        bond = response_bridge(1)['bond_candidate']
        for sign in (-1, 1):
            flow = cayley_flow(NATIVE_R, sign)
            returned = return_operator(NATIVE_K, NATIVE_R, flow)
            self.assertEqual(mul(returned, bond), scale(bond, sign))

    def test_zero_response_and_invalid_source_contracts(self):
        zero = response_bridge(0)
        self.assertEqual(zero['commutator'], ZERO)
        self.assertFalse(zero['exact_cut'])
        with self.assertRaises(ValueError):
            response_bridge(-1)
        with self.assertRaises(ValueError):
            closure_budget(2, 1)
        with self.assertRaises(ValueError):
            period_cut_condition(0, 1)
        with self.assertRaises(TypeError):
            response_bridge(0.5)

    def test_changed_upstream_bytes_are_refused_before_execution(self):
        pins = self.packet['checked_runtime_sha256']
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for name in pins:
                path = root/name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes((self.source_root/name).read_bytes())
            target = root/next(iter(pins))
            target.write_bytes(target.read_bytes()+b'\n// changed\n')
            with self.assertRaisesRegex(ValueError, 'source hash mismatch'):
                check_runtime_pins(root)


if __name__ == '__main__':
    unittest.main()
