"""Full-space vacuum counts, noncommuting steps and scale-inference controls."""
from fractions import Fraction as F
import copy
import json
import unittest

import yc6_vacuum_step as v


def scale(a, value):
    return [[value*x for x in row] for row in a]


class YC6Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = v.jsonable(v.run())
        cls.records = json.loads((v.HERE/'YC6_ENDPOINTS.json').read_text())['endpoints']
        cls.rows = v.envelopes(cls.records)

    def test_frozen_packet_and_source_pins(self):
        self.assertEqual(self.result, json.loads((v.HERE/'YC6_RESULT.json').read_text()))
        self.assertTrue(all(self.result['checks'].values()))

    def test_free_vacuum_gap_and_all_reflection_blocks(self):
        self.assertEqual(set(v.m.resolved_blocks(0)), {(0,0,0), (0,1,1), (1,0,1), (1,1,0)})
        self.assertEqual(v.m.counts(0, 0, 8, 'lower', True)[0], 1)
        self.assertEqual(self.result['point_gap_enclosures']['0'], ['8', '8'])
        # A centre-odd bottom at 3 is not a vacuum-sector excitation.
        self.assertEqual(v.m.c.counts(1, 0, 3, False)[:2], (0, 1))
        self.assertGreater(8, 3)

    def test_ritz_upper_above_hidden_floor_does_not_use_an_inverse(self):
        record = self.records[-1]
        self.assertEqual(record['excited_upper_method'], 'ritz')
        z = F(record['levels'][1][1])
        self.assertGreater(z, 24)
        negative, zero, _ = v.m.c.counts(0, 14, z, False)
        self.assertGreaterEqual(negative+zero, 2)
        with self.assertRaises(ValueError):
            v.m.counts(0, 14, z, 'upper', True)

    def test_false_excited_lower_rejected(self):
        record = copy.deepcopy(self.records[16])
        record['levels'][1][0] = record['levels'][1][1]
        with self.assertRaises(AssertionError):
            v.validate_endpoint(record)

    def test_domain_shape_and_method_gates(self):
        for key, value in (('theta', '15'), ('excited_upper_method', 'discard_hidden'), ('levels', [['0', '1']])):
            record = copy.deepcopy(self.records[0])
            record[key] = value
            with self.assertRaises(ValueError):
                v.validate_endpoint(record)
        with self.assertRaises(AssertionError):
            v.envelopes(self.records[:-1])
        with self.assertRaises(AssertionError):
            v.envelopes(self.records[::-1])

    def test_old_floor_is_carried_forward_without_claiming_comparison_monotonicity(self):
        final = self.rows[-1]
        self.assertEqual(final['floor_source_theta'], '25/2')
        self.assertGreater(F(final['second_floor']), F(self.records[-1]['levels'][1][0]))
        at_source = next(r for r in self.records if r['theta'] == final['floor_source_theta'])
        self.assertEqual(F(final['second_floor']), F(at_source['levels'][1][0]))
        self.assertLess(F(final['floor_source_theta']), F(final['theta']))
        # A later endpoint may improve a later floor; it cannot alter earlier rows.
        changed = copy.deepcopy(self.records)
        changed[-1]['levels'][1][0] = '23'
        self.assertEqual(v.envelopes(changed)[:-1], self.rows[:-1])

    def test_all_cells_are_covered_and_bound_is_uniform(self):
        cells = self.result['coupling_cells']
        self.assertEqual(len(cells), 112)
        self.assertEqual([F(c['left']) for c in cells], [F(n,8) for n in range(112)])
        self.assertEqual([F(c['right']) for c in cells], [F(n+1,8) for n in range(112)])
        self.assertEqual(min(F(c['gap_lower']) for c in cells), F(2113,2000))
        self.assertTrue(all(F(c['gap_lower']) >= 1 for c in cells))

    def test_actual_vacuum_gap_is_not_globally_nondecreasing(self):
        self.assertLess(F(self.result['point_gap_enclosures']['2'][1]), 8)
        # CM2's positive asymptotic coefficient is retained separately.
        cm2 = json.loads((v.ROOT/'physics/cm2/CM2_RESULT.json').read_text())
        self.assertTrue(cm2['all_pass'])

    def test_step_interval_division_preserves_all_four_corners(self):
        gaps = {F(r['theta']): list(map(F, r['gap'])) for r in self.rows}
        for step in self.result['fourfold_coupling_steps']:
            theta = F(step['theta'])
            lower, upper = map(F, step['gap_ratio'])
            for numerator in gaps[4*theta]:
                for denominator in gaps[theta]:
                    self.assertLessEqual(lower, numerator/denominator)
                    self.assertGreaterEqual(upper, numerator/denominator)
            self.assertEqual(list(map(F, step['normalized_ratio_cubed'])), [lower**3/4, upper**3/4])

    def test_exact_noncommuting_resolvent_step(self):
        d0, w = [[F(5), F(1)], [F(1), F(8)]], [[F(2), F(1)], [F(1), F(3)]]
        eye = [[F(1), F(0)], [F(0), F(1)]]
        self.assertNotEqual(v.mm(d0,w), v.mm(w,d0))
        r1 = v.inverse(v.add(v.add(d0,w),eye,-1))
        r4 = v.inverse(v.add(v.add(d0,w,4),eye,-1))
        self.assertEqual(v.add(r4,r1,-1), scale(v.mm(r4,v.mm(w,r1)),-3))
        self.assertEqual(v.m.c.inertia(v.add(r1,r4,-1))[0], 0)
        with self.assertRaises(ValueError):
            v.inverse([[1,2],[2,4]])

    def test_fourfold_schur_identity_retains_the_final_interaction(self):
        control = v.noncommuting_step_control()
        old, new = control['sigma_before'], control['sigma_after']
        moment = v.add(old,new,-1)
        cross = scale(moment,F(1,3))
        kinetic = [[F(3),F(0)],[F(0),F(7)]]
        potential = [[F(4),F(1)],[F(1),F(5)]]
        s1 = v.add(v.add(kinetic,potential),old,-1)
        s4 = v.add(v.add(kinetic,potential,4),new,-16)
        transported = v.add(v.add(v.add(s1,potential,3),old,-15),cross,48)
        self.assertEqual(s4, transported)
        self.assertNotEqual(s4, v.add(v.add(s1,potential,3),old,-15))

    def test_alternating_remainder_at_a_large_step(self):
        control = v.noncommuting_step_control()
        old, actual = control['sigma_before'], control['sigma_after']
        partial = [[F(0),F(0)],[F(0),F(0)]]
        for n, coefficient in enumerate(control['unsigned_taylor']):
            self.assertEqual(v.m.c.inertia(coefficient)[0], 0)
            partial = v.add(partial,coefficient,(-3)**n)
            signed_rest = scale(v.add(actual,partial,-1),(-1)**(n+1))
            self.assertEqual(v.m.c.inertia(signed_rest)[0], 0)
            # ||W||<=6, D0>=4, z=1; no small-step/convergent-series claim.
            remainder_bound = v.add(scale(old,6**(n+1)),signed_rest,-1)
            self.assertEqual(v.m.c.inertia(remainder_bound)[0], 0)

    def test_ground_state_dirichlet_transform_control(self):
        h = [[F(2),F(-1)],[F(-1),F(1,2)]]
        psi = [[F(1)],[F(2)]]
        self.assertEqual(v.mm(h,psi), [[0],[0]])
        self.assertEqual(v.m.c.inertia(h), (0,1,1))
        for f in ((F(-4),F(1)),(F(2),F(3)),(F(7,3),F(-2,5))):
            embedded = [[f[0]],[2*f[1]]]
            energy = v.mm(v.transpose(embedded),v.mm(h,embedded))[0][0]/5
            self.assertEqual(energy,F(2,5)*(f[0]-f[1])**2)
            mean = (f[0]+4*f[1])/5
            variance = ((f[0]-mean)**2+4*(f[1]-mean)**2)/5
            self.assertEqual(energy/variance,F(5,2))

    def test_constant_information_does_not_fix_a_dimensionless_gap(self):
        probabilities = (F(1,5),F(4,5))
        for n in (1,2,4,8):
            epsilon = F(1,n*n)
            rates = [[-2*epsilon,2*epsilon],[epsilon/2,-epsilon/2]]
            self.assertEqual(v.mm([list(probabilities)],rates), [[0,0]])
            gap = -sum(rates[i][i] for i in range(2))
            self.assertEqual(gap,F(5,2*n*n))
        self.assertEqual(probabilities,(F(1,5),F(4,5)))

    def test_normalized_vanishing_does_not_imply_absolute_vanishing(self):
        for n in (2,3,5):
            # theta=n^6: theta^(1/3)=n^2; absolute examples shrink/stay/grow.
            absolute = [F(1,n**6),F(1),F(n)]
            self.assertEqual([x/n**2 for x in absolute], [F(1,n**8),F(1,n**2),F(1,n)])
            self.assertLess(absolute[0],1)
            self.assertEqual(absolute[1],1)
            self.assertGreater(absolute[2],1)

    def test_cutoff_and_fixed_box_limits_have_opposite_behavior_in_the_control(self):
        # a=2^(-n^3), g^2=1/n^3: g^(2/3)=1/n, logarithmic UV running.
        for n in range(1,7):
            cell_rate = F(2**(n**3),n)
            next_rate = F(2**((n+1)**3),n+1)
            self.assertGreaterEqual(next_rate/cell_rate,64)
            self.assertLess(F(1,n+1),F(1,n))  # unmatched bare-g core at fixed L=1
        self.assertEqual(v.scale_controls()['compact_physical_g_power'],F(2,3))

    def test_refinement_step_and_coupling_step_are_different(self):
        x, beta_increment = F(5),F(1)
        theta, refined = x*x,(x+beta_increment)**2
        self.assertEqual((theta,refined),(25,36))
        self.assertNotEqual(refined,4*theta)
        # A fixed physical mass has a power law in a, even with running g.
        mass,a,b = F(7,3),F(1,16),F(2)
        self.assertEqual((a/b*mass)/(a*mass),1/b)

    def test_ym_line_shift_cancels_in_the_gap_only(self):
        theta = F(4)
        h0,h1 = F(7),F(14)
        a0,a1 = h0/4-3*theta/4,h1/4-3*theta/4
        self.assertEqual(a1-a0,(h1-h0)/4)
        self.assertNotEqual(h0,4*a0)
        self.assertEqual(h0,4*a0+3*theta)

    def test_claim_boundaries(self):
        boundary = self.result['claim_boundary']
        for key in ('one_site_only','vacuum_sector_only','full_hidden_return_retained','ground_state_transform_analytic_identity'):
            self.assertTrue(boundary[key])
        for key,value in boundary.items():
            if key not in ('one_site_only','vacuum_sector_only','full_hidden_return_retained','ground_state_transform_analytic_identity'):
                self.assertFalse(value,key)


if __name__ == '__main__':
    unittest.main()
