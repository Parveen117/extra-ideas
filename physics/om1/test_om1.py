"""OM1: metric covariance, complete residuals and approximate-vacuum controls."""
from fractions import Fraction as F
import unittest

import om1_observer_memory as m


class ObserverMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_complete_packet(self):
        self.assertGreaterEqual(len(self.result['checks']), 22)
        for name, value in self.result['checks'].items():
            self.assertTrue(value, name)

    def test_approximate_vacuum_floor_is_not_second_eigenvalue(self):
        # h=diag(1,3,5), psi=(2,1,0)/sqrt(5); Q compression has levels 2.6,5.
        eta, variance = F(7, 5), F(16, 25)
        low, delta, fidelity = m.cut_bounds(eta, variance, F(3))
        self.assertEqual((low, delta, fidelity), (F(1), F(13, 5), F(4, 5)))
        self.assertLess(delta, 3)  # Using z itself as the hidden floor is false.
        mu = [variance*delta**k for k in range(3)]
        ret = m.return_from_blocks(m.blocks(mu, 1), F(1), delta)
        self.assertEqual(ret['lower'], eta-1)
        self.assertEqual(ret['upper'], eta-1)
        self.assertEqual(ret['residual'], 0)

    def test_nonorthogonal_hidden_coordinates_preserve_the_return(self):
        p = m.trial_packet(3)
        block = m.blocks(p['mu'], 2)
        energy = F(51867, 10000)
        expected = m.return_from_blocks(block, energy, p['delta'])
        r = [[F(1), F(2)], [F(0), F(3)]]
        rt = m.cm.transpose(r)
        transformed = {'norm': block['norm']}
        for name in ('s', 'h', 'h2'):
            transformed[name] = m.cm.matmul(rt, m.cm.matmul(block[name], r))
        for name in ('b', 'c'):
            transformed[name] = [row[0] for row in m.cm.matmul(rt, [[v] for v in block[name]])]
        actual = m.return_from_blocks(transformed, energy, p['delta'])
        for name in ('lower', 'upper', 'residual', 'bare_upper'):
            self.assertEqual(actual[name], expected[name])
        wrong = dict(block, s=[[F(1), F(0)], [F(0), F(1)]])
        with self.assertRaises(ValueError):
            m.return_from_blocks(wrong, energy, p['delta'])

    def test_full_core_residual_by_direct_projected_action(self):
        p = m.trial_packet(3)
        energy = F(51867, 10000)
        ret = m.return_bounds(p, energy, 2)
        psi, hpsi = p['vectors'][:2]
        source = m.cr.padd(hpsi, psi, -p['eta'])
        def hidden_action(v):
            hv = m.cm.scalar_h(v, 3, p['width'])
            return m.cr.padd(hv, psi, -p['inner'](psi, hv)/p['norm'])
        dsource = hidden_action(source)
        trial = m.cr.padd({e: ret['coefficients'][0]*v for e, v in source.items()},
                         dsource, ret['coefficients'][1])
        residual = m.cr.padd(m.cr.padd(source, hidden_action(trial), -1), trial, energy)
        self.assertEqual(p['inner'](psi, residual), 0)
        self.assertEqual(p['inner'](residual, residual)/p['norm'], ret['residual'])

    def test_schur_endpoints_have_the_required_exact_signs(self):
        for d in (3, 4):
            p = m.trial_packet(d)
            row = self.result['numbers'][str(d)]['return_refinements'][-1]
            low, high = F(row['ground_lower']), F(row['ground_upper'])
            self.assertGreater(p['eta']-low-m.return_bounds(p, low, 2)['upper'], 0)
            self.assertLess(p['eta']-high-m.return_bounds(p, high, 2)['lower'], 0)

    def test_unproved_or_insufficient_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            m.cut_bounds(F(3), F(1), F(3))
        with self.assertRaises(ValueError):
            m.cut_bounds(F(1), F(-1), F(3))
        p = m.trial_packet(4)
        with self.assertRaises(ValueError):
            m.return_bounds(p, p['delta'], 2)
        # A poor approximate vacuum may have a valid floor below its own eta.
        _, bad_floor, _ = m.cut_bounds(F(2), F(2), F(3))
        self.assertLess(bad_floor, 2)

    def test_cuts_and_continuum_claims_remain_distinct(self):
        boundary = self.result['claim_boundary']
        self.assertFalse(boundary['computed_reading_is_exact_vacuum'])
        self.assertEqual(boundary['CM1_67_reading_source_rank_26'], 'UNCHANGED_DIFFERENT_CUT')
        self.assertEqual(boundary['nonconstant_mode_matching'], 'OPEN')
        self.assertEqual(boundary['continuum_mass_gap'], 'OPEN')


if __name__ == '__main__':
    unittest.main()
