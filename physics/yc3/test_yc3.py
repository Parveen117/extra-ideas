"""Independent structural and exact-arithmetic checks for YC3."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
import unittest

import yc3_sector_splitting as m


def rational_rank(rows):
    rows = [list(map(F, row)) for row in rows]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [x/scale for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            scale = rows[i][col]
            if scale:
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def invariant_dimension(k):
    """Kernel on all 4^k degree-one tensors, without assuming an invariant basis."""
    words = list(product(range(4), repeat=k))
    indices = {word: i for i, word in enumerate(words)}
    matrices = []
    for a, b in combinations((1, 2, 3), 2):
        matrix = [[0]*len(words) for _ in words]
        for col, word in enumerate(words):
            for link, coord in enumerate(word):
                if coord in (a, b):
                    changed = list(word)
                    changed[link] = a if coord == b else b
                    matrix[indices[tuple(changed)]][col] += 1 if coord == b else -1
        matrices.extend(matrix)
    return len(words)-rational_rank(matrices)


def plane_rotation(p, i, j):
    """x_i partial_j - x_j partial_i, independent of the spherical generator."""
    out = {}
    for e, c in p.items():
        for a, b, sign in ((i, j, 1), (j, i, -1)):
            if e[b]:
                f = list(e)
                f[a] += 1
                f[b] -= 1
                f = tuple(f)
                out[f] = out.get(f, F(0))+sign*e[b]*c
    return {e: c for e, c in out.items() if c}


def casimir(p):
    out = {}
    for link in range(3):
        for i, j in combinations(range(4*link, 4*link+4), 2):
            out = m.padd(out, plane_rotation(plane_rotation(p, i, j), i, j), -1)
    return out


class YC3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_packet_and_source_pins(self):
        stored = json.loads((Path(__file__).parent/'YC3_RESULT.json').read_text())
        self.assertEqual(stored, self.result)
        self.assertTrue(all(self.result['checks'].values()))

    def test_independent_gauge_kernel_completeness(self):
        self.assertEqual([invariant_dimension(k) for k in range(4)], [1, 1, 2, 5])
        for k in range(4):
            self.assertTrue(all(m.sector_data(k)['gram'][i][i] > 0
                                for i in range(len(m.lowest_basis(k)))))

    def test_sphere_constraints_and_normalization(self):
        self.assertEqual(m.mean({m.ZERO: F(1)}), 1)
        for link in range(3):
            constraint = {m.ZERO: F(-1)}
            for i in range(4*link, 4*link+4):
                constraint = m.padd(constraint, m.mono(i, i))
                self.assertEqual(m.mean(m.mono(i, i)), F(1, 4))
                self.assertEqual(m.mean(m.mono(i, i, i, i)), F(1, 8))
            self.assertEqual(m.inner(constraint, constraint), 0)

    def test_independent_radial_angular_scalar_moments(self):
        # Under scalar-product trial squared, odd links have Beta(3/2,3/2)
        # vector radius squared; even links have Beta(3/2,1/2).
        # Angular cross-square means are 2/3, squares 8/15, shared pairs 4/9.
        for k in range(4):
            first = [F(1, 2) if i < k else F(3, 4) for i in range(3)]
            second = [F(5, 16) if i < k else F(5, 8) for i in range(3)]
            mean = F(4, 3)*sum(first[i]*first[j] for i, j in combinations(range(3), 2))
            raw_second = 4*(F(8, 15)*sum(second[i]*second[j] for i, j in combinations(range(3), 2))
                           + F(8, 9)*sum(second[i]*first[(i+1) % 3]*first[(i+2) % 3] for i in range(3)))
            data = m.sector_data(k)
            self.assertEqual(mean, data['a'])
            self.assertEqual(raw_second-mean**2, data['residual_squares'][0])

    def test_spherical_generator_against_rotation_casimir(self):
        for k in range(4):
            for p in (m.lowest_basis(k)[0], m.sector_data(k)['residuals'][0]):
                diff = m.padd(m.h0(p), casimir(p), -1)
                # Ambient polynomials can differ by sphere ideal elements.
                self.assertEqual(m.inner(diff, diff), 0)

    def test_h0_self_adjoint_on_residuals(self):
        for k in range(4):
            for f in m.sector_data(k)['residuals']:
                p = m.h0(f)
                self.assertEqual(m.inner(f, m.h0(p)), m.inner(p, p))

    def test_rank_one_cut_would_miss_degenerate_hidden_state(self):
        f, lost = m.lowest_basis(2)
        self.assertEqual(m.inner(f, lost), 0)
        self.assertEqual(m.inner(lost, m.h0(lost))/m.inner(lost, lost), 6)
        self.assertLess(6, m.sector_data(2)['e']+m.sector_data(2)['delta'])
        self.assertNotEqual(m.reflected_signature(f), m.reflected_signature(lost))

    def test_so3_triple_product_is_not_discarded(self):
        t = m.triple()
        self.assertGreater(m.inner(t, t), 0)
        self.assertEqual(m.inner(t, m.h0(t))/m.inner(t, t), 9)
        for a, b in combinations((1, 2, 3), 2):
            self.assertFalse(m.gauge_generator(t, a, b))

    def test_complete_residual_weights_and_kinetic_filter(self):
        expected = [{8: F(3, 4), 16: F(7, 16)},
                    {8: F(25, 72), 12: F(1, 4), 16: F(7, 48), 20: F(11, 72)},
                    {8: F(1, 9), 12: F(25, 72), 20: F(11, 72), 24: F(1, 24)},
                    {12: F(1, 3), 24: F(1, 8)}]
        for k, weights in enumerate(expected):
            self.assertEqual({s: w for s, w in m.source_components(k)['weights'].items() if w}, weights)
            d = m.sector_data(k)
            p = r = d['residuals'][0]
            norm = d['gram'][0][0]
            # Moment identities use the independent rotation Casimir, not
            # the polynomial projectors which produced the packet weights.
            for power in range(3):
                self.assertEqual(m.inner(r, p)/norm, sum(w*s**power for s, w in weights.items()))
                p = m.padd(casimir(p), p, -d['e'])
        self.assertEqual(m.source_components(0)['coefficient'], F(31, 256))
        self.assertNotEqual(m.source_components(0)['coefficient'], m.sector_data(0)['b2'])

    def test_exact_endpoint_gap_and_curvature(self):
        self.assertEqual(m.splitting(1, 1), (F(707, 300), F(249, 92)))
        self.assertEqual(m.splitting(1, 2), (F(65, 54), F(47, 14)))
        self.assertEqual(m.source_components(0)['coefficient']-m.source_components(1)['coefficient'], F(77, 1920))

    def test_finite_window_bounds_and_ordering(self):
        # Sampling is a regression check; the note proves the full interval
        # by monotonicity and exact endpoint margins.
        for n in range(41):
            theta = F(n, 20)
            lower_gap, upper_gap = m.splitting(1, theta)
            self.assertGreaterEqual(lower_gap, F(65, 54))
            self.assertLessEqual(lower_gap, upper_gap)
            for k in range(4):
                l, u = m.lower_upper(k, theta)
                refined_l, refined_u = m.second_order_window(k, theta)
                self.assertLessEqual(l, refined_l)
                self.assertLessEqual(refined_l, refined_u)
                self.assertLessEqual(refined_u, u)
            for k in (2, 3):
                self.assertGreater(m.lower_upper(k, theta)[0], m.lower_upper(1, theta)[1])

    def test_zero_coupling_and_cubic_order(self):
        for k in range(4):
            self.assertEqual(m.lower_upper(k, 0), (3*k, 3*k))
            self.assertEqual(m.second_order_window(k, 0), (3*k, 3*k))
            d = m.sector_data(k)
            t = F(1, 100)
            ratio = m.cubic_remainder(k, t)/t**3
            self.assertEqual(ratio, 6*d['residual_squares'][0]/(d['delta']*(d['delta']-t*d['a'])))

    def test_domains_reject_unsupported_continuations(self):
        for k in (-1, 4, True, F(1, 2)):
            with self.assertRaises(ValueError):
                m.lowest_basis(k)
        for k in range(4):
            with self.assertRaises(ValueError):
                m.lower_upper(k, -1)
            with self.assertRaises(ValueError):
                m.lower_upper(k, F(m.sector_data(k)['delta'])/m.sector_data(k)['a'])
            with self.assertRaises(ValueError):
                m.cubic_remainder(k, F(201, 100))
        with self.assertRaises(ValueError):
            m.moment((0,)*11)
        with self.assertRaises(ValueError):
            m.splitting(0, 1)

    def test_claim_boundary(self):
        boundary = self.result['claim_boundary']
        self.assertTrue(boundary['complete_hidden_space_floor_used'])
        for key in ('static_Gibbs_compass_substituted_for_quantum_resolvent',
                    'weak_limit_centre_splitting_rate', 'larger_lattice_or_continuum_mass_gap',
                    'DR2_column_count_identified_with_volume',
                    'twisted_trace_is_unique_possible_sector_observer'):
            self.assertFalse(boundary[key])


if __name__ == '__main__':
    unittest.main()
