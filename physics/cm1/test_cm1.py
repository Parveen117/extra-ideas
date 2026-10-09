"""CM1: physical-entry checks, metric controls and fail-closed boundaries."""
from fractions import Fraction as F
import unittest

import cm1_core_sector_memory as m


def derivative(p, j):
    out = {}
    for e, v in p.items():
        if e[j]:
            f = list(e)
            f[j] -= 1
            out[tuple(f)] = v*e[j]
    return out


def entry_mean(p, w):
    return sum(m.tc.free_mean({e: v})/w**(sum(e)//2) for e, v in p.items())


class CoreSectorMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_complete_exact_packet(self):
        self.assertGreaterEqual(len(self.result['checks']), 50)
        for name, passed in self.result['checks'].items():
            self.assertTrue(passed, name)

    def test_turning_energy_from_entry_gradient_form(self):
        # A second route to an actual matrix element: gradient energy rather
        # than the reduced generator, with Gaussian moments in C entries.
        for d in (3, 4):
            w = m.cr.RATE[d]
            p = m.spin_to_entries({(1, 1, 0, 0): F(1)}, d)
            kinetic = {}
            for j in range(3*d):
                e = [0]*(3*d)
                e[j] = 1
                dg = m.tc.padd(derivative(p, j), m.tc.pmul({tuple(e): F(1)}, p), -w)
                kinetic = m.tc.padd(kinetic, m.tc.pmul(dg, dg), F(1, 2))
            x = m.cr.to_entries(m.E2, m.tc.entry_polys(d))
            full = m.tc.padd(kinetic, m.tc.pmul(x, m.tc.pmul(p, p)))
            keys, s, h = m.spin_ladder(d, w, 3)
            i = keys.index((1, 1, 0, 0))
            self.assertEqual(entry_mean(full, w), h[i][i])
            self.assertEqual(entry_mean(m.tc.pmul(p, p), w), s[i][i])

    def test_column_swap_changes_turning_sign(self):
        for d in (3, 4):
            p = m.spin_to_entries({(1, 1, 0, 0): F(1), (2, 0, 0, 0): F(2)}, d)
            swapped = {}
            for e, v in p.items():
                f = list(e)
                for row in range(3):
                    f[row*d], f[row*d+1] = f[row*d+1], f[row*d]
                swapped[tuple(f)] = v
            self.assertEqual(swapped, {e: -v for e, v in p.items()})
            self.assertEqual(entry_mean(p, F(1)), 0)

    def test_source_gram_is_congruent_under_nonorthogonal_relabeling(self):
        packet = m.source_packet(3, m.cr.RATE[3], 2)
        n = len(packet['s'])
        r = [[F(i == j) for j in range(n)] for i in range(n)]
        r[0][1] = 2
        r[1][2] = -1
        def relabel(a):
            return m.matmul(m.transpose(r), m.matmul(a, r))
        s, h, t = [relabel(packet[k]) for k in ('s', 'h', 'second')]
        correct = m.subtract(t, m.matmul(h, m.solve_spd(s, h)))
        self.assertEqual(correct, relabel(packet['gram']))
        wrong = m.subtract(t, m.matmul(h, h))
        self.assertNotEqual(wrong, correct)

    def test_hidden_floor_is_a_required_hypothesis(self):
        packet = m.source_packet(3, m.cr.RATE[3], 2)
        for floor in (F(414, 100), F(6)):
            with self.assertRaises(ValueError):
                m.floor_count(packet, F(6), floor)
        low, high, count = m.required_floor(packet, F(6))
        self.assertGreater(m.floor_count(packet, F(6), low), 1)
        self.assertLessEqual(count, 1)
        self.assertLessEqual(high-low, F(1, 1000))

    def test_no_promotion_from_ritz_or_conditional_floor(self):
        boundary = self.result['claim_boundary']
        self.assertEqual(boundary['excited_energy_lower'], 'OPEN')
        self.assertEqual(boundary['quantitative_core_gap_lower'], 'OPEN')
        self.assertEqual(boundary['continuum_mass_gap'], 'OPEN')
        self.assertIn('NOT_AN_ACTUAL_HIDDEN_FLOOR', boundary['conditional_schur_floor_demand'])
        self.assertIn('NOT_A_FINITE_CERTIFICATE', boundary['qualitative_core_gap'])

    def test_gram_solver_rejects_indefinite_or_nonsymmetric_metric(self):
        for s in ([[F(1), F(0)], [F(0), F(-1)]],
                  [[F(1), F(1)], [F(0), F(1)]]):
            with self.assertRaises(ValueError):
                m.solve_spd(s, [[F(1)], [F(1)]])


if __name__ == '__main__':
    unittest.main()
