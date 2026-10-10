"""YC4 independent arithmetic, completeness and omission controls."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import copy
import json
import unittest

import yc4_harmonic_return as m


def character_singlets(ns):
    """Weight-0 minus weight-1 multiplicity, independent of triangle counting."""
    weights = {0: 1}
    for n in ns:
        updated = {}
        for old, count in weights.items():
            for weight in range(-n, n+1):
                updated[old+weight] = updated.get(old+weight, 0)+count*(n-abs(weight)+1)
        weights = updated
    return weights.get(0, 0)-weights.get(1, 0)


def mm(a, b):
    return [[sum(a[i][l]*b[l][j] for l in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return list(map(list, zip(*a)))


def rotation(p, i, j):
    out = {}
    for e, c in p.items():
        for a, b, sign in ((i, j, 1), (j, i, -1)):
            if e[b]:
                f = list(e)
                f[a] += 1
                f[b] -= 1
                f = tuple(f)
                out[f] = out.get(f, F(0))+c*sign*e[b]
    return {e: c for e, c in out.items() if c}


class YC4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()
        cls.records = json.loads((m.HERE/'YC4_ENDPOINTS.json').read_text())['endpoints']

    def test_packet_and_source_pins(self):
        self.assertEqual(self.result, json.loads((m.HERE/'YC4_RESULT.json').read_text()))
        self.assertTrue(all(self.result['checks'].values()))

    def test_independent_character_rank_formula(self):
        for ns in product(range(5), repeat=3):
            self.assertEqual(m.spin_count(ns), character_singlets(ns))
        self.assertEqual([sum(m.spin_count(ns) for ns in m.shells(k)) for k in range(4)], [13, 21, 32, 23])

    def test_no_parity_allowed_shell_below_floor_is_omitted(self):
        for k in range(4):
            expected = {ns for ns in product(range(8), repeat=3)
                        if all(n % 2 == int(i < k) for i, n in enumerate(ns))
                        and sum(n*(n+2) for n in ns) < 3*k+24}
            self.assertEqual(set(m.shells(k)), expected)
            # Any degree >=8 alone contributes at least 80, beyond every cut.
            self.assertGreater(8*10, 3*k+24)

    def test_fast_haar_pairing_against_unpacked_engine(self):
        samples = [m.y.lowest_basis(3)[-1], m.y.potential(), m.y.mono(0, 0, 4),
                   m.shell((2, 2, 0))[-1][0], m.shell((3, 1, 0))[-1][0]]
        for i, p in enumerate(samples):
            for q in samples[:i+1]:
                self.assertEqual(m.inner(p, q), m.y.inner(p, q))
        self.assertEqual(m.inner({}, samples[0]), 0)
        with self.assertRaises(ValueError):
            m.inner(m.y.mono(*([0]*16)), samples[0])

    def test_link_generator_against_rotation_casimir(self):
        p = m.shell((3, 1, 2))[-1][0]
        for link in range(3):
            out = {}
            for i in range(4*link, 4*link+4):
                for j in range(4*link, i):
                    out = m.y.padd(out, rotation(rotation(p, i, j), i, j), -1)
            diff = m.y.padd(m.hi(p, link), out, -1)
            self.assertEqual(m.inner(diff, diff), 0)

    def test_full_source_gram_by_explicit_projected_polynomial(self):
        for k in range(4):
            block = max(m.matrices(k).values(), key=lambda b: len(b['norms']))
            n = len(block['norms'])
            coeff = [F((i % 3)-1, i+1) for i in range(n)]
            p = {}
            for c, vector in zip(coeff, block['vectors']):
                p = m.y.padd(p, vector[0], c)
            r = m.y.pmul(m.y.potential(), p)
            for i, vector in enumerate(block['vectors']):
                projection = sum(block['V'][i][j]*coeff[j] for j in range(n))/block['norms'][i]
                r = m.y.padd(r, vector[0], -projection)
            for vector in block['vectors']:
                self.assertEqual(m.inner(r, vector[0]), 0)
            expected = sum(coeff[i]*block['B'][i][j]*coeff[j] for i in range(n) for j in range(n))
            self.assertEqual(m.inner(r, r), expected)
            self.assertGreater(expected, 0)

    def test_first_residual_is_retained_exactly_only_in_claimed_sectors(self):
        for k in range(4):
            signature = tuple(int(i < k) for i in range(3))
            block = m.matrices(k)[signature]
            if k < 2:
                self.assertEqual(block['B'][0][0], 0)
            else:
                self.assertGreater(block['B'][0][0], 0)

    def test_omitting_a_shell_invalidates_the_hidden_floor(self):
        # The old vacuum cut leaves a degree-two class harmonic at energy 8.
        p, norm, _ = m.shell((2, 0, 0))[0]
        self.assertEqual(m.inner(p, {m.y.ZERO: F(1)}), 0)
        self.assertEqual(m.inner(p, m.y.h0(p))/norm, 8)
        self.assertLess(8, 24)

    def test_exact_inertia_zero_diagonal_and_congruence(self):
        self.assertEqual(m.inertia([[0, 1], [1, 0]]), (1, 0, 1))
        self.assertEqual(m.inertia([[0, 0], [0, 0]]), (0, 2, 0))
        diagonal = [[F(-3), 0, 0], [0, 0, 0], [0, 0, F(5)]]
        change = [[1, 2, 1], [0, 1, -1], [0, 0, 2]]
        transformed = mm(transpose(change), mm(diagonal, change))
        self.assertEqual(m.inertia(transformed), (1, 1, 1))
        with self.assertRaises(ValueError):
            m.inertia([[1, 2], [0, 1]])

    def test_dropped_return_would_falsely_promote_ritz_to_lower(self):
        # H0=diag(0,12), V=ones is positive; H=[[1,1],[1,13]].
        z = F(19, 20)
        self.assertEqual(m.inertia([[1-z, 1], [1, 13-z]])[0], 1)
        self.assertGreater(1-z, 0)  # wrong lower if the lost source is dropped
        self.assertLess(1-z-1/(12-z), 0)  # full-source bound correctly refuses
        lower = F(9, 10)
        self.assertGreater(1-lower-1/(12-lower), 0)  # a true certified lower

    def test_metric_cannot_be_replaced_by_identity(self):
        block = m.matrices(0)[(0, 0, 0)]
        z = F(4)
        actual = m.pencil(block, 0, 2, z, True)
        self.assertEqual(m.inertia(actual)[0], 0)
        wrong = copy.deepcopy(actual)
        for i, norm in enumerate(block['norms']):
            wrong[i][i] += z*(norm-1)
        self.assertGreater(m.inertia(wrong)[0], 0)

    def test_bad_endpoint_and_outside_window_are_rejected(self):
        record = copy.deepcopy(self.records[16])  # theta=2
        record['ground'][1][0] = record['ground'][1][1]
        with self.assertRaises(AssertionError):
            m.validate_endpoint(record)
        record = copy.deepcopy(self.records[-1])
        record['theta'] = '9'
        with self.assertRaises(ValueError):
            m.validate_endpoint(record)
        record = copy.deepcopy(self.records[0])
        record['second_lower'].pop()
        with self.assertRaises(ValueError):
            m.validate_endpoint(record)

    def test_floor_and_parameter_gates(self):
        for k in range(4):
            block = next(iter(m.matrices(k).values()))
            with self.assertRaises(ValueError):
                m.pencil(block, k, 1, 3*k+24, True)
            with self.assertRaises(ValueError):
                m.pencil(block, k, -1, 0, True)
        smaller = next(iter(m.matrices(0, 8).values()))
        self.assertEqual(smaller['floor'], 8)
        with self.assertRaises(ValueError):
            m.pencil(smaller, 0, 1, 8, True)
        with self.assertRaises(ValueError):
            m.pencil(smaller, 1, 1, 0, True)
        for value in (0, -1, F(1, 2), True):
            with self.assertRaises(ValueError):
                m.shells(0, value)

    def test_whole_interval_and_multiplicity_margins(self):
        self.assertEqual(len(self.result['coupling_cells']), 64)
        self.assertEqual([F(r['theta']) for r in self.records], [F(n, 8) for n in range(65)])
        for c in self.result['coupling_cells']:
            self.assertGreaterEqual(F(c['gap_lower']), F(1, 5))
            self.assertGreater(F(c['other_centre_margin']), 0)
            self.assertGreater(F(c['second_level_margin']), 0)
        self.assertEqual(self.result['first_excitation_multiplicity'], 3)

    def test_improves_yc3_and_crosses_old_window(self):
        new = list(map(F, self.result['point_gap_enclosures']['2']))
        old = m.y.splitting(1, 2)
        self.assertGreater(new[0], old[0])
        self.assertLess(new[1], old[1])
        self.assertLess(new[1]-new[0], F(12, 1000))
        self.assertGreater(F(self.result['point_gap_enclosures']['8'][0]), F(3, 10))

    def test_no_continuum_or_weak_limit_upgrade(self):
        boundary = self.result['claim_boundary']
        self.assertTrue(boundary['one_site_only'])
        self.assertTrue(boundary['all_gauge_and_centre_sectors_covered'])
        self.assertTrue(boundary['uncomputed_complement_floor_proved'])
        for key in ('floating_point_used_to_certify', 'weak_limit_absolute_splitting_proved',
                    'volume_uniform_or_4d_gap', 'formal_proof_assistant_verification'):
            self.assertFalse(boundary[key])


if __name__ == '__main__':
    unittest.main()
