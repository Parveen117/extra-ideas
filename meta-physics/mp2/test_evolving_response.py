"""Exact controls and failure cases, supplementary to the seven written proofs."""
from fractions import Fraction as F
import itertools
import unittest
import evolving_response as e


class EvolvingResponseTests(unittest.TestCase):
    def setUp(self):
        self.u = e.u
        self.V = self.u.a.diag((4, 1))
        self.O = [[F(3, 5), F(-4, 5)], [F(4, 5), F(3, 5)]]

    def test_cofactor_geometry_grid(self):
        # 624 nonzero integer matrices, including rank-one and negative area.
        u = self.u
        count = 0
        for values in itertools.product(range(-2, 3), repeat=4):
            if not any(values):
                continue
            V = [list(values[:2]), list(values[2:])]
            s = e.state(V)
            self.assertEqual(e.cofactor(e.cofactor(V)), V)
            self.assertEqual(e.inner(e.cofactor(V), e.cofactor(V)), s['B'])
            self.assertEqual(e.inner(s['plus'], s['minus']), 0)
            self.assertEqual(s['a_plus']-s['a_minus'], 2*s['d'])
            v = e.velocity(V, F(2, 3), F(-3, 7))
            self.assertEqual(e.inner(e.cofactor(V), v), 0)
            self.assertEqual(2*e.inner(V, v), -F(4, 3)*(s['B']**2-4*s['d']**2)/s['B'])
            self.assertEqual(e.cofactor(u.mm(self.O, V)), u.mm(self.O, e.cofactor(V)))
            count += 1
        self.assertEqual(count, 624)

    def test_projection_is_nearest_tangent_velocity(self):
        u = self.u
        V = [[2, 3], [-1, 4]]
        s = e.state(V)
        Z = e.velocity(V, 3)
        residual = u.a.add(Z, u.a.scale(V, 3))
        for W in [u.I, u.K, u.L, u.R]:
            tangent = u.a.add(W, u.a.scale(e.cofactor(V), -e.inner(e.cofactor(V), W)/s['B']))
            self.assertEqual(e.inner(residual, tangent), 0)
            moved = u.a.add(residual, tangent)
            self.assertEqual(e.inner(moved, moved), e.inner(residual, residual)+e.inner(tangent, tangent))

    def test_hessian_integrability_and_native_marker(self):
        for V in [self.V, [[1, 2], [3, -1]], [[1, 2], [2, 4]], [[0, 0], [0, 0]]]:
            A, B = e.hessian_jets(V)
            self.assertEqual(A[0][1], B[0][0])
            self.assertEqual(A[1][1], B[0][1])
            d = e.response(V)
            self.assertEqual([list(v[:2]) for v in d['vectors']], self.u.a.transpose(V))
            self.assertEqual(d['marker'], self.u.det2(V)/F(2))
            self.assertEqual(d['C'], e.tensor(V))

    def test_common_positive_hessian_patch(self):
        u = self.u
        for V in [self.V, u.mm(self.O, self.V), [[1, 2], [3, -1]]]:
            self.assertLessEqual(e.inner(V, V), 25)
            for x, y in [(F(1, 30), 0), (F(-1, 30), 0), (0, F(1, 30)),
                         (0, F(-1, 30)), (F(1, 60), F(-1, 60))]:
                H = e.hessian(V, x, y)
                self.assertTrue(u.a.psd(u.a.add(H, u.a.scale(u.I, F(-1, 2)))))
                self.assertTrue(u.a.psd(u.a.add(u.a.scale(u.I, F(3, 2)), H, -1)))

    def test_analytic_budget_and_equal_branch_loss(self):
        for V in [self.V, [[1, 2], [3, -1]], [[2, 0], [0, 2]]]:
            s = e.state(V)
            v = e.velocity(V, 3, 7)
            vplus = self.u.a.scale(self.u.a.add(v, e.cofactor(v)), F(1, 2))
            vminus = self.u.a.scale(self.u.a.add(v, e.cofactor(v), -1), F(1, 2))
            Bdot = 2*e.inner(V, v)
            self.assertEqual(2*s['B']*Bdot, -12*(s['B']**2-4*s['d']**2))
            self.assertEqual(2*e.inner(s['plus'], vplus), Bdot/2)
            self.assertEqual(2*e.inner(s['minus'], vminus), Bdot/2)

    def test_exact_pair_orbit_and_record_recovery(self):
        step = e.pair_step(self.V, F(41, 50), F(3, 10), self.O)
        self.assertEqual(step['d'], 4)
        self.assertEqual(step['B'], F(881, 100))
        self.assertEqual(step['m'], F(819, 200))
        self.assertEqual(e.recover_rational(step['V'], self.O, step['m']), self.V)
        q = (step['B']**2-64)/F(17**2-64)
        self.assertTrue(0 < q < 1)  # exp(-4*kappa*time), finite positive time.
        s = e.state(step['V'])
        self.assertEqual(s['a_plus']+step['m'], F(25, 2))
        self.assertEqual(s['a_minus']+step['m'], F(9, 2))

    def test_negative_orientation_orbit(self):
        V = self.u.a.diag((4, -1))
        step = e.pair_step(V, F(3, 10), F(41, 50), self.O)
        self.assertEqual(step['d'], -4)
        self.assertEqual(e.recover_rational(step['V'], self.O, step['m']), V)

    def test_record_cannot_recover_a_lost_direction(self):
        with self.assertRaises(ValueError):
            e.recover_rational(self.u.a.diag((2, 2)), self.u.I, F(9, 2))

    def test_invalid_flow_and_pair_inputs_rejected(self):
        bad = [lambda: e.state([[0, 0], [0, 0]]), lambda: e.state([[1]]),
               lambda: e.velocity(self.V, -1), lambda: e.velocity(self.V, 0.5),
               lambda: e.pair_step(self.V, F(1, 2), F(1, 2)),
               lambda: e.pair_step(self.V, 0, 1),
               lambda: e.pair_step(self.V, 1, 1, self.u.K),
               lambda: e.pair_step([[1, 2], [2, 4]], 1, 1),
               lambda: e.rational_root(F(2))]
        for run in bad:
            with self.assertRaises((ValueError, TypeError)):
                run()

    def test_balanced_and_zero_relaxation_controls(self):
        u = self.u
        for V in [u.a.diag((2, 2)), u.a.diag((2, -2))]:
            self.assertEqual(e.velocity(V, 7, 3), u.a.scale(u.mm(u.R, V), 3))
            self.assertEqual(e.balances(V, 7, 3)['m_dot'], 0)
        self.assertEqual(e.velocity(self.V, 0, 2), u.a.scale(u.mm(u.R, self.V), 2))

    def test_euler_does_not_preserve_area(self):
        u = self.u
        v = e.velocity(self.V)
        W = u.a.add(self.V, u.a.scale(v, F(1, 10)))
        self.assertEqual(e.inner(e.cofactor(self.V), v), 0)
        self.assertNotEqual(u.det2(W), 4)

    def test_uniform_curvature_floor_is_rank_two(self):
        u = self.u
        for V in [self.V, [[1, 2], [3, -1]], u.mm(self.O, self.V)]:
            s, d = e.state(V), e.response(V)
            g = e.floor(abs(s['d']), s['B'])
            self.assertEqual(g['even'], d['beta_floor'])
            self.assertTrue(u.a.psd(u.a.add(u.a.scale(d['G'], F(1, 4)),
                                          u.a.scale(u.I, -g['even']))))
            self.assertFalse(u.a.psd(u.a.add(d['C'], u.a.scale(u.a.identity(3), -g['even']))))

    def test_rotating_generators_need_time_ordering(self):
        u = self.u
        C0 = e.tensor(self.V)
        C1 = e.tensor(u.mm(self.O, self.V))
        A, B = u.q.operator(2, C0), u.q.operator(2, C1)
        self.assertNotEqual(u.mm(A, B), u.mm(B, A))
        S, T = e.resolvent(C0, 2, F(1, 10)), e.resolvent(C1, 2, F(1, 10))
        self.assertNotEqual(u.mm(S, T), u.mm(T, S))

    def test_ordered_resolvent_common_reference_contraction(self):
        u = self.u
        h, rate = F(1, 10), e.floor(4, 17)['full']
        tensors = [e.tensor(self.V), e.tensor(u.mm(self.O, self.V)),
                   e.tensor(e.pair_step(self.V, F(41, 50), F(3, 10))['V'])]
        for degree in [1, 2, 3]:
            _, G = u.z.gram(degree)
            T = u.a.identity(len(G))
            for C in tensors:
                S = e.resolvent(C, degree, h)
                T = u.mm(S, T)
            bound = (1+h*rate)**(-2*len(tensors))
            deficit = u.a.add(u.a.scale(G, bound), u.mm(u.mm(u.a.transpose(T), G), T), -1)
            self.assertTrue(u.a.psd(deficit))

    def test_noise_source_product_identity(self):
        u = self.u
        C = e.tensor(u.mm(self.O, self.V))
        for f in [u.p.COORD[0], u.p.mul(u.p.COORD[0], u.p.COORD[2]), e.work_witness()[0]]:
            gamma = u.z.gamma(f, f, C)
            rhs = u.p.add(u.p.scale(u.p.mul(f, u.q.lap(f, C)), 2),
                          u.p.scale(gamma, -2))
            self.assertEqual(u.q.lap(u.p.mul(f, f), C), rhs)
            self.assertGreaterEqual(u.y.phi(gamma), 0)
            self.assertEqual(u.y.phi(u.q.lap(u.p.mul(f, f), C)), 0)

    def test_moving_work_prevents_false_energy_monotonicity(self):
        u = self.u
        f, C, Cdot, b = e.work_witness()
        self.assertEqual(Cdot, u.a.diag((F(-2400, 17), F(150, 17), 0)))
        self.assertEqual(u.y.phi(u.p.mul(f, f)), F(1, 12))
        self.assertEqual(b['energy'], F(1, 24))
        self.assertEqual(b['work'], F(25, 34))
        self.assertEqual(b['energy_dot'], F(283, 408))
        self.assertEqual(b['squared_norm_dot'], F(-1, 12))

    def test_work_formula_independent_energy_difference(self):
        u = self.u
        f, C, Cdot, b = e.work_witness()
        h = F(1, 1000)
        plus, minus = u.a.add(C, u.a.scale(Cdot, h)), u.a.add(C, u.a.scale(Cdot, -h))
        self.assertEqual((u.z.energy(f, plus)-u.z.energy(f, minus))/(2*h), b['work'])
        Lf = u.q.lap(f, C)
        fp, fm = u.p.add(f, u.p.scale(Lf, -h)), u.p.add(f, u.p.scale(Lf, h))
        self.assertEqual((u.z.energy(fp, C)-u.z.energy(fm, C))/(2*h), -b['dissipation'])

    def test_fixed_area_observers_choose_different_shapes(self):
        z = self.u.z
        full0, even0, _ = z.rates((F(0), F(2), F(2)))
        full1, even1, _ = z.rates((F(0), F(32, 25), F(25, 8)))
        self.assertEqual(full0, 1)
        self.assertEqual(full1, F(881, 800))
        self.assertGreater(full1, full0)
        self.assertLess(even1, even0)
        # At the exact crossing (b,c)=(1,3), full=1, even=1, d^2=12.
        self.assertEqual(z.rates((F(0), F(1), F(3)))[0], 1)
        for b in [F(1, 4), F(1, 2), F(3, 4), F(1), F(5, 4), F(3, 2)]:
            c = 3/b
            self.assertGreaterEqual(c, b)
            self.assertLessEqual(z.rates((F(0), b, c))[0], 1)

    def test_integrated_error_survival_gate(self):
        b = e.robust_budget(self.V, 5, F(1, 10))
        self.assertEqual(b['D'], F(699, 200))
        self.assertEqual(b['B'], F(2601, 100))
        self.assertGreater(b['full'], 0)
        for r, eps in [(4, 0), (5, -1), (5, 1)]:
            with self.assertRaises(ValueError):
                e.robust_budget(self.V, r, eps)

    def test_open_area_and_energy_sources(self):
        E = [[F(1, 7), F(2, 9)], [F(-3, 11), F(1, 5)]]
        b = e.balances(self.V, 3, 2, E)
        self.assertEqual(b['d_dot'], e.inner(e.cofactor(self.V), E))
        self.assertEqual(b['open_energy'], e.inner(self.V, E))
        self.assertGreaterEqual(b['m_dot'], 0)

    def test_budget_only_and_area_only_both_fail(self):
        # Rational samples of the two already-known YM54 obstruction families.
        for n in [2, 4, 16, 64]:
            losing = e.state([[1, 0], [0, F(1, n)]])
            escaping = e.state([[n, 0], [0, F(1, n)]])
            self.assertLess(losing['B'], 2)
            self.assertEqual(losing['d'], F(1, n))
            self.assertEqual(escaping['d'], 1)
            self.assertEqual(self.u.z.rates((F(0), F(1, 2*n*n), F(n*n, 2)))[1], F(1, 2*n*n))
        with self.assertRaises(ValueError):
            e.floor(0, 2)

    def test_native_scaling_keeps_clock_selection_explicit(self):
        u = self.u
        V = u.a.scale(self.V, 3)
        self.assertEqual(e.tensor(V), u.a.scale(e.tensor(self.V), 9))
        self.assertEqual(e.velocity(V, 2, 3), u.a.scale(e.velocity(self.V, 2, 3), 3))
        self.assertEqual(e.floor(36, 153)['full'], 9*e.floor(4, 17)['full'])


if __name__ == '__main__':
    raise SystemExit('Run verify.py with --publications-root so exact dependency pins are checked first')
