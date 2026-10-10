"""Independent parity, domain, full-return and metric controls for YC24."""
from fractions import Fraction as Q
from itertools import product
import hashlib
import json
import unittest

import sympy as sp

import yc24_parity_response_tower as y


class YC24Tests(unittest.TestCase):
    def test_return_words_at_odd_order_never_reach_the_vacuum_signature(self):
        for n in (1,3,5):
            for word in product(y.FACE_MASKS,repeat=n):
                mask = 0
                for face in word:
                    mask ^= face
                self.assertNotEqual(mask,0)
        for word in ((3,3),(3,6,12,9)):
            mask = 0
            for face in word:
                mask ^= face
            self.assertEqual(mask,0)

    def test_removed_empty_signature_has_floor_eight_not_zero(self):
        floors = dict(zip(y.B_MASKS,y.B_FLOORS))
        for ns in product(range(6),repeat=4):
            if not any(ns) or sum(ns)%2:
                continue
            mask = sum((n%2)<<j for j,n in enumerate(ns))
            energy = sum(n*(n+2) for n in ns)
            self.assertGreaterEqual(energy,floors.get(mask,6))
        self.assertEqual(floors[0],8)

    def test_actual_bridge_geometry_has_no_within_side_face(self):
        for mask in (*y.E_MASKS,*y.B_MASKS):
            targets = {mask^face for face in y.FACE_MASKS}
            self.assertEqual(targets,set(y.B_MASKS if mask in y.E_MASKS else y.E_MASKS))

    def test_nonsymmetric_four_face_ray_keeps_source_norm_and_full_bounds(self):
        row = y.response_bounds(Q(3,4),Q(1,3),xi=(1,Q(1,2),-Q(1,3),0))
        self.assertEqual(row['v0'],Q(49,144))
        self.assertGreater(row['alpha'],0)
        self.assertLess(row['rho'],1)
        self.assertLess(row['metric_cap'],y.response_bounds(1,Q(1,2))['metric_cap'])

    def test_domain_errors_are_not_reported_as_certified_towers(self):
        for kwargs in ({'t':Q(3,2)}, {'z':6}, {'z':Q(11,2)},
                       {'xi':(1,1,1,2)}, {'xi':(1,1)}, {'order':-1}):
            with self.assertRaises(ValueError):
                y.response_bounds(**kwargs)

    def test_source_free_limit_keeps_identity_metric_and_zero_tail(self):
        for kwargs in ({'t':0},{'xi':(0,0,0,0)}):
            row = y.response_bounds(**kwargs)
            self.assertEqual(row['return_cap'],0)
            self.assertEqual(row['metric_cap'],1)
            self.assertEqual(row['energy_tail'],0)
            self.assertEqual(row['metric_error'],0)

    def test_resolvent_leaves_an_individual_source_and_mixes_neutral_readings(self):
        row = y.finite_control()
        self.assertNotEqual(row['Z'][1,0],0)
        self.assertNotEqual(row['Z'][2,0],0)
        self.assertEqual(row['R'][:,0].dot(row['R'][:,1]),0)
        self.assertNotEqual(row['sigma'][0,1],0)

    def test_energy_tower_is_positive_but_metric_tail_can_be_indefinite(self):
        row = y.finite_control(order=0)
        self.assertTrue(y.y20.positive_semidefinite(row['sigma']-row['sigmaN']))
        diff = row['metric']-row['metricN']
        self.assertLess(sp.factor(diff.det()),0)

    def test_graph_metric_is_not_derivative_of_truncated_return(self):
        z = sp.symbols('z',real=True)
        row = y.finite_control(z=z,order=0)
        trial_derivative = sp.eye(2)+row['sigmaN'].diff(z)
        self.assertFalse(y.y20.equal(trial_derivative,row['metricN']))
        self.assertTrue(y.y20.equal(-row['F'].diff(z),row['metric']))

    def test_exact_geometric_remainder_retains_the_full_inverse(self):
        q = sp.Rational(9,16)
        row = y.finite_control(t=sp.Rational(3,4),z=sp.Rational(1,2),order=2)
        K, A, R = row['step'],row['A'],row['R']
        remainder = q**4*R.H*K**3*(sp.eye(2)-q*K).inv()*A.inv()*R
        self.assertTrue(y.y20.equal(row['sigma']-row['sigmaN'],remainder))

    def test_full_strength_tail_caps_decrease_strictly_with_order(self):
        rows = [y.response_bounds(order=n) for n in range(7)]
        for first,second in zip(rows,rows[1:]):
            self.assertGreater(first['energy_tail'],second['energy_tail'])
            self.assertGreater(first['metric_error'],second['metric_error'])
        self.assertLess(rows[3]['energy_tail'],Q(444,100000))
        self.assertLess(rows[3]['metric_error'],Q(302,100000))

    def test_smaller_couplings_preserve_both_block_spectral_reserves(self):
        for theta in (0,Q(1,2),1,Q(3,2),2):
            end = y.block_bounds(theta,1)
            for t in (0,Q(1,4),Q(2,3),1):
                row = y.block_bounds(theta,t)
                self.assertGreaterEqual(row['spectral_reserve'],end['spectral_reserve'])
                self.assertGreater(row['spectral_reserve'],0)
                self.assertEqual(row['full_gap'],Q(3,4))

    def test_lattice_transfer_changes_only_the_physical_reference_floor(self):
        for kind in ('tube','force'):
            before,after = y.y22.join_bounds(kind),y.join_bounds(kind)
            for key in ('full_factor_floor','incidence','beta','seed','contraction',
                        'relative_return','ball_reserve','external_cap','full_gap'):
                self.assertEqual(before[key],after[key])
            self.assertGreater(after['physical_gap'],before['physical_gap'])
        self.assertEqual(y.join_bounds('force')['physical_gap'],Q(6476,3125))

    def test_eight_sector_schur_comparison_matches_the_scalar_reserve(self):
        z = sp.Rational(12,5)
        full = y.scalar_hidden_comparison()-z*sp.eye(8)
        schur = full[:4,:4]-full[:4,4:]*full[4:,4:].inv()*full[4:,:4]
        row = y.response_bounds(1,Q(12,5))
        self.assertEqual((schur*sp.ones(4,1))[0],sp.Rational(row['alpha']))
        self.assertTrue(y.y20.positive_semidefinite(schur-sp.Rational(row['alpha'])*sp.eye(4)))

    def test_source_pins_and_scope(self):
        record = json.loads((y.HERE/'YC24_RESULT.json').read_text())
        for path,digest in record['source_pins'].items():
            self.assertEqual(hashlib.sha256((y.ROOT/path).read_bytes()).hexdigest(),digest)
        self.assertTrue(all(record['checks'].values()))
        self.assertTrue(record['claims']['all_hidden_harmonics_retained'])
        self.assertFalse(record['claims']['local_tower_ratio_is_RG_contraction'])
        self.assertFalse(record['claims']['continuum_mass_gap'])


if __name__ == '__main__':
    unittest.main()
