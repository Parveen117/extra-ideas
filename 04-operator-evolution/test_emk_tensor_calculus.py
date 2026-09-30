"""Focused checks of R7's tensor typing and load-bearing closure distinctions."""

import importlib.util
from itertools import product
from math import prod
import os
from pathlib import Path
import unittest

import emk_tensor_calculus as tc
from aghora_return import cayley_flow

Q = tc.Q
V = tc.Slot('V', 1)
DUAL = tc.Slot('V', -1)
ZERO2 = tc.scale(tc.I, 0)


def load_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sample(slots, start=1):
    return tc.Tensor(slots, tuple(Q(i, i+1) for i in
                                 range(start, start + prod(s.dimension for s in slots))))


def matrix_action(a, t):
    return tc.apply_slot(t, 0, a)


class EMKTensorCalculusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(os.environ['PUBLICATIONS_SOURCE_ROOT'])
        folder = root / 'papers/emk-ugd-algebra/certificates'
        cls.native = load_source('emk2_r7', folder/'emk2_native_carrier.py')
        cls.rtc = load_source('emkt1_r7', folder/'emkt1_master_tensor_and_time.py')
        cls.ordered = load_source('emkt2_r7', folder/'emkt2_time_ordered_transport.py')

    def test_01_consumes_native_basis_and_ordered_loop(self):
        for values in ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0),
                       (0, 0, 0, 1), (2, 3, 5, 7)):
            a, b, c, d = values
            expected = tc.add(tc.add(tc.scale(tc.I, a), tc.scale(tc.K, b)),
                              tc.add(tc.scale(tc.R, c), tc.scale(tc.RK, d)))
            self.assertEqual(expected, tc.matrix(self.native.block(*values)))
        self.assertEqual(tc.K, tc.matrix(self.rtc.cut_swap()))
        up = tc.matrix(self.ordered.CH_UP[0])
        lo = tc.matrix(self.ordered.CH_LO[0])
        a, b = Q(1, 3), Q(2, 5)
        u = tc.add(tc.I, tc.scale(up, a))
        w = tc.add(tc.I, tc.scale(lo, b))
        loop = tc.mul(tc.mul(u, w), tc.mul(tc.inverse(u), tc.inverse(w)))
        self.assertEqual(loop, tc.matrix(self.ordered.group_loop(a, up, b, lo)))
        # Same upstream order defect, now lifted to a genuine tensor type.
        slots = (V, DUAL, V)
        lifted = tc.induced_matrix(slots, {'V': loop})
        self.assertNotEqual(lifted, tc.identity(8))

    def test_02_separate_e_m_frames_need_a_pairing_bridge(self):
        e, ed = tc.Slot('e', 1), tc.Slot('e', -1)
        m, md = tc.Slot('m', 1), tc.Slot('m', -1)
        omega, v = tc.Tensor((ed,), (2, 3)), tc.Tensor((m,), (5, 7))
        bridge = tc.tensor_from_matrix(((1, 2), (0, 1)), e, md)
        def pairing(w, b, z):
            t = tc.tensor_product(tc.tensor_product(w, b), z)
            return tc.contract(tc.contract(t, 0, 1), 0, 1).data[0]
        changes = {'e': tc.matrix(((2, 1), (0, 1))),
                   'm': tc.matrix(((1, 0), (3, 2)))}
        self.assertEqual(pairing(omega, bridge, v), pairing(
            tc.transport(omega, changes), tc.transport(bridge, changes),
            tc.transport(v, changes)))
        with self.assertRaises(ValueError):
            tc.contract(tc.tensor_product(omega, v), 0, 1)

    def test_03_transport_composes_and_contraction_is_natural(self):
        a = {'V': tc.matrix(((2, 1), (1, 1)))}
        b = {'V': tc.matrix(((1, 0), (-2, 3)))}
        ab = {'V': tc.mul(b['V'], a['V'])}
        for slots in ((V,), (DUAL,), (V, DUAL), (V, DUAL, V, DUAL)):
            t = sample(slots)
            self.assertEqual(tc.transport(tc.transport(t, a), b), tc.transport(t, ab))
            if len(slots) >= 2:
                self.assertEqual(tc.contract(tc.transport(t, a), 0, 1),
                                 tc.transport(tc.contract(t, 0, 1), a))
        x, y = sample((V,)), sample((DUAL, V), 3)
        self.assertEqual(tc.transport(tc.tensor_product(x, y), a),
                         tc.tensor_product(tc.transport(x, a), tc.transport(y, a)))

    def test_04_tensor_transport_is_not_a_linear_KIR_algebra_map(self):
        slots = (V, V)
        lift = lambda a: tc.induced_matrix(slots, {'V': a})
        self.assertNotEqual(lift(tc.add(tc.I, tc.R)), tc.add(lift(tc.I), lift(tc.R)))
        self.assertEqual(tc.mul(lift(tc.R), lift(tc.R)), tc.identity(4))
        self.assertEqual(tc.mul(lift(tc.R), lift(tc.K)),
                         tc.mul(lift(tc.K), lift(tc.R)))
        self.assertNotEqual(tc.mul(lift(tc.R), lift(tc.R)), tc.scale(tc.identity(4), -1))

    def test_05_rank_parity_and_metric_return_blindness(self):
        for rank in range(1, 5):
            for signs in product((-1, 1), repeat=rank):
                slots = tuple(tc.Slot('V', s) for s in signs)
                size = 2**rank
                self.assertEqual(tc.induced_matrix(slots, {'V': tc.scale(tc.I, -1)}),
                                 tc.scale(tc.identity(size), (-1)**rank))
        p = tc.Move('0', 'a', (('V', tc.R),)).then(
            tc.Move('a', '2', (('V', tc.K),)))
        q = tc.Move('0', 'b', (('V', tc.K),)).then(
            tc.Move('b', '2', (('V', tc.R),)))
        loop = p.then(q.reverse())
        self.assertEqual(dict(loop.operators)['V'], tc.scale(tc.I, -1))
        audit = tc.return_report(loop, ((DUAL, DUAL), (V, DUAL), (V,)))
        self.assertEqual(audit['tensor_type_returns'], [True, True, False])
        self.assertFalse(audit['full_unledgered_return'])

    def test_06_generator_lift_preserves_brackets_but_not_primitive_square(self):
        slots = (V, V)
        lift = lambda a: tc.induced_matrix(slots, {'V': a}, infinitesimal=True)
        self.assertEqual(tc.commutator(lift(tc.R), lift(tc.K)),
                         lift(tc.commutator(tc.R, tc.K)))
        self.assertNotEqual(tc.mul(lift(tc.R), lift(tc.R)), tc.scale(tc.identity(4), -1))
        metric = tc.tensor_from_matrix(tc.I, DUAL, DUAL)
        self.assertEqual(tc.generator_action(metric, {'V': tc.R}), tc.tensor_scale(metric, 0))

    def test_07_seam_slot_decomposition_is_unique_and_covariant(self):
        t = sample((V, DUAL, V))
        change = tc.matrix(((2, 1), (1, 1)))
        kp = tc.mul(tc.mul(change, tc.K), tc.inverse(change))
        def blocks(value, k):
            out = []
            for signs in product((-1, 1), repeat=3):
                piece = value
                for axis, sign in enumerate(signs):
                    p = tc.scale(tc.add(tc.I, tc.scale(k, sign)), Q(1, 2))
                    piece = tc.apply_slot(piece, axis,
                                          p if value.slots[axis].variance == 1 else tc.transpose(p))
                out.append(piece)
            return out
        pieces, changed = blocks(t, tc.K), blocks(tc.transport(t, {'V': change}), kp)
        total = tc.tensor_scale(t, 0)
        for original, after in zip(pieces, changed):
            total = tc.tensor_add(total, original)
            self.assertEqual(tc.transport(original, {'V': change}), after)
        self.assertEqual(total, t)

    def test_08_exterior_forms_and_metric_variance_changes(self):
        t = sample((DUAL, DUAL))
        symmetric = tc.symmetrize(t)
        alternating = tc.symmetrize(t, alternating=True)
        self.assertEqual(tc.tensor_add(symmetric, alternating), t)
        a, b = tc.Tensor((DUAL,), (1, 2)), tc.Tensor((DUAL,), (3, 5))
        self.assertEqual(tc.wedge(a, b).data, (0, -1, 1, 0))
        self.assertEqual(tc.wedge(a, b), tc.tensor_scale(tc.wedge(b, a), -1))
        change = tc.matrix(((2, 1), (0, 1)))
        self.assertEqual(tc.transport(tc.wedge(a, b), {'V': change}),
                         tc.wedge(tc.transport(a, {'V': change}), tc.transport(b, {'V': change})))
        # The core is finite dimensional, not limited to the two-mode example.
        d3 = tc.Slot('W', -1, 3)
        forms = [tc.Tensor((d3,), tuple(int(i == j) for i in range(3))) for j in range(3)]
        self.assertEqual(tc.wedge(tc.wedge(forms[0], forms[1]), forms[2]).at((0, 1, 2)), 1)
        metric = tc.matrix(((2, 1), (1, 3)))
        v = sample((V, V))
        lowered = tc.change_variance(v, 0, metric)
        self.assertEqual(tc.change_variance(lowered, 0, metric), v)
        h_next = tc.mul(tc.mul(tc.transpose(tc.inverse(change)), metric), tc.inverse(change))
        self.assertEqual(tc.transport(lowered, {'V': change}),
                         tc.change_variance(tc.transport(v, {'V': change}), 0, h_next))

    def test_09_exact_finite_product_rule_and_ledger_coupling(self):
        operators = {'V': tc.matrix(((1, 2), (1, 3)))}
        x0, x1 = sample((V,)), sample((V,), 4)
        y0, y1 = sample((DUAL,)), sample((DUAL,), 6)
        dx, dy = tc.increment(x0, x1, operators), tc.increment(y0, y1, operators)
        expected = tc.tensor_add(tc.tensor_product(dx, tc.transport(y0, operators)),
                                 tc.tensor_product(x1, dy))
        actual = tc.increment(tc.tensor_product(x0, y0), tc.tensor_product(x1, y1), operators)
        self.assertEqual(actual, expected)
        lx, ly = tc.Tensor((V,), (1, 0)), tc.Tensor((DUAL,), (0, 2))
        lxy = tc.tensor_scale(actual, 0)
        dx, dy = tc.increment(x0, x1, operators, lx), tc.increment(y0, y1, operators, ly)
        expected = tc.tensor_add(tc.tensor_product(dx, tc.transport(y0, operators)),
                                 tc.tensor_product(x1, dy))
        coupling = tc.tensor_sub(tc.tensor_add(
            tc.tensor_product(lx, tc.transport(y0, operators)), tc.tensor_product(x1, ly)), lxy)
        self.assertNotEqual(coupling, tc.tensor_scale(coupling, 0))
        self.assertEqual(tc.increment(tc.tensor_product(x0, y0), tc.tensor_product(x1, y1),
                                     operators, lxy), tc.tensor_add(expected, coupling))

    def test_10_smooth_jet_covariance_and_noncoordinate_curvature(self):
        a = tc.add(tc.R, tc.scale(tc.K, Q(2, 3)))
        g, dg = tc.matrix(((2, 1), (1, 1))), tc.matrix(((1, 2), (0, 3)))
        slots = (V, DUAL, V)
        t, dt = sample(slots), sample(slots, 4)
        t_next = tc.transport(t, {'V': g})
        velocity = tc.mul(dg, tc.inverse(g))
        dt_next = tc.tensor_add(tc.generator_action(t_next, {'V': velocity}),
                                tc.transport(dt, {'V': g}))
        a_next = tc.sub(tc.mul(tc.mul(g, a), tc.inverse(g)), velocity)
        self.assertEqual(tc.covariant_derivative(t_next, dt_next, {'V': a_next}),
                         tc.transport(tc.covariant_derivative(t, dt, {'V': a}), {'V': g}))
        expected = tc.sub(tc.add(tc.sub(tc.K, tc.R), tc.commutator(a, tc.K)), tc.RK)
        self.assertEqual(tc.curvature(a, tc.K, tc.K, tc.R, tc.RK), expected)

    def test_11_curvature_lifts_and_bianchi_do_not_require_flatness(self):
        for slots in ((V,), (V, DUAL), (DUAL, DUAL, V)):
            a = tc.induced_matrix(slots, {'V': tc.R}, infinitesimal=True)
            b = tc.induced_matrix(slots, {'V': tc.K}, infinitesimal=True)
            f = tc.induced_matrix(slots, {'V': tc.curvature(tc.R, tc.K)}, infinitesimal=True)
            self.assertEqual(tc.commutator(a, b), f)
        edges = {'01': tc.add(tc.I, tc.matrix(((0, 1), (0, 0)))),
                 '12': tc.add(tc.I, tc.matrix(((0, 0), (1, 0)))),
                 '23': tc.R, '02': tc.I, '13': tc.K,
                 '03': tc.matrix(((2, 1), (1, 1)))}
        self.assertNotEqual(tc.triangle_curvature(edges['01'], edges['12'], edges['02']), ZERO2)
        self.assertEqual(tc.tetrahedron_bianchi(edges), ZERO2)
        for slots in ((V, DUAL), (V, V, V)):
            lifted = {name: tc.induced_matrix(slots, {'V': value}) for name, value in edges.items()}
            self.assertEqual(tc.tetrahedron_bianchi(lifted), tc.scale(tc.identity(2**len(slots)), 0))

    def test_12_exterior_covariant_square_is_curvature_and_lift_order_matters(self):
        u01 = tc.matrix(((1, 1), (0, 1)))
        u12 = tc.matrix(((1, 0), (1, 1)))
        u02 = tc.I
        s0, s1, s2 = (sample((V,), n) for n in (1, 4, 7))
        d01 = tc.tensor_sub(s1, tc.transport(s0, {'V': u01}))
        d12 = tc.tensor_sub(s2, tc.transport(s1, {'V': u12}))
        d02 = tc.tensor_sub(s2, tc.transport(s0, {'V': u02}))
        dd = tc.tensor_add(tc.tensor_sub(d12, d02), tc.transport(d01, {'V': u12}))
        f = tc.triangle_curvature(u01, u12, u02)
        self.assertEqual(dd, matrix_action(f, s0))
        self.assertNotEqual(dd, tc.tensor_scale(dd, 0))
        # Finite curvature is a difference AFTER transport is lifted.
        slots = (V, V)
        lifted_f = tc.triangle_curvature(*(tc.induced_matrix(slots, {'V': u})
                                         for u in (u01, u12, u02)))
        self.assertNotEqual(lifted_f, tc.induced_matrix(slots, {'V': f}))

    def test_13_full_continuous_RK_has_no_fixed_symmetric_metric(self):
        self.assertEqual(tc.invariant_symmetric_forms((tc.R,)), (tc.I,))
        self.assertEqual(tc.invariant_symmetric_forms((tc.K,)),
                         (tc.matrix(((-1, 0), (0, 1))),))
        self.assertEqual(tc.invariant_symmetric_forms((tc.R, tc.K)), ())
        # Discrete R and K ARE isometries of I; the obstruction concerns both flows.
        self.assertEqual(tc.metric_defect(tc.R, tc.I, tc.I), ZERO2)
        self.assertEqual(tc.metric_defect(tc.K, tc.I, tc.I), ZERO2)
        boost = cayley_flow(tc.K, Q(1, 3))
        self.assertNotEqual(tc.metric_defect(boost, tc.I, tc.I), ZERO2)
        self.assertEqual(tc.mul(tc.mul(tc.transpose(boost), tc.EPSILON), boost), tc.EPSILON)
        for u in (tc.R, tc.K, boost, tc.matrix(((2, 3), (5, 7)))):
            det = u[0][0]*u[1][1] - u[0][1]*u[1][0]
            self.assertEqual(tc.mul(tc.mul(tc.transpose(u), tc.EPSILON), u),
                             tc.scale(tc.EPSILON, det))
        for a in (tc.R, tc.K, tc.RK, tc.I):
            trace = a[0][0] + a[1][1]
            self.assertEqual(tc.add(tc.mul(tc.transpose(a), tc.EPSILON),
                                    tc.mul(tc.EPSILON, a)), tc.scale(tc.EPSILON, trace))

    def test_14_projection_loses_curvature_and_same_generator_cut_join_commute(self):
        p, q = tc.matrix(((1, 0), (0, 0))), tc.matrix(((0, 0), (0, 1)))
        compressed = tc.commutator(tc.mul(tc.mul(p, tc.R), p), tc.mul(tc.mul(p, tc.K), p))
        full = tc.mul(tc.mul(p, tc.commutator(tc.R, tc.K)), p)
        leakage = tc.sub(tc.mul(tc.mul(tc.mul(tc.mul(p, tc.R), q), tc.K), p),
                         tc.mul(tc.mul(tc.mul(tc.mul(p, tc.K), q), tc.R), p))
        self.assertEqual(tc.sub(full, compressed), leakage)
        self.assertEqual(compressed, ZERO2)
        self.assertNotEqual(full, ZERO2)
        for u in (tc.I, tc.R, tc.matrix(((2, 1), (1, 1)))):
            cut = tc.sub(u, tc.inverse(u))
            join = tc.add(u, tc.inverse(u))
            self.assertEqual(tc.commutator(cut, join), ZERO2)

    def test_15_path_integration_reversal_and_independent_sheet_memory(self):
        x = sample((V, DUAL))
        steps = []
        nodes = ('0', '1', '2', '0')
        for i, u in enumerate((tc.R, tc.K, tc.matrix(((2, 1), (1, 1))))):
            move = tc.Move(nodes[i], nodes[i+1], (('V', u),), i - 1)
            steps.append((move, sample(x.slots, i+2), sample(x.slots, i+5)))
        end, combined, correction = tc.integrate_path(x, steps)
        self.assertEqual(end, tc.tensor_sub(tc.transport(x, dict(combined.operators)), correction))
        self.assertEqual(combined.sheet, 0)
        self.assertEqual(dict(combined.then(combined.reverse()).operators)['V'], tc.I)
        c0, c1 = tc.matrix(((2, 1), (0, 1))), tc.matrix(((1, 0), (2, 1)))
        changed = steps[0][0].reframe({'V': c0}, {'V': c1}, source_sheet=5, target_sheet=8)
        self.assertEqual(changed.sheet, steps[0][0].sheet + 3)
        loop = tc.Move('0', '0', (('V', tc.I),), 2)
        audit = tc.return_report(loop, ((V,), (DUAL, DUAL)))
        self.assertEqual(audit['tensor_type_returns'], [True, True])
        self.assertFalse(audit['full_unledgered_return'])
        ledgered = tc.return_report(loop, ((V,),), sheet_ledger=2)
        self.assertTrue(ledgered['full_return_after_sheet_ledger'])
        self.assertFalse(ledgered['full_unledgered_return'])

    def test_16_refuses_untyped_operations_and_singular_changes(self):
        with self.assertRaises(TypeError):
            tc.Tensor((V,), (0.5, 1))
        with self.assertRaises(ValueError):
            tc.transport(sample((V,)), {'V': ZERO2})
        with self.assertRaises(ValueError):
            tc.contract(sample((V, V)), 0, 1)
        with self.assertRaises(ValueError):
            tc.change_variance(sample((V,)), 0, ZERO2)
        with self.assertRaises(ValueError):
            tc.symmetrize(sample((V, DUAL)))
        with self.assertRaises(TypeError):
            tc.Move('0', '1', (('V', tc.I),), Q(1, 2))
        with self.assertRaises(ValueError):
            tc.Move('0', '1', (('V', tc.I),)).then(tc.Move('2', '3', (('V', tc.I),)))
        with self.assertRaises(ValueError):
            tc.wedge(sample((V,)), sample((V,)))
        with self.assertRaises(ValueError):
            tc.Tensor((V, tc.Slot('V', -1, 3)), (1,)*6)


if __name__ == '__main__':
    unittest.main()
