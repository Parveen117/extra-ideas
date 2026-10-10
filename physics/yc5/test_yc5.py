"""Independent algebra, full-source and noncommuting-return controls for YC5."""
from fractions import Fraction as F
import copy
import json
import unittest

import yc5_resolved_return as m


def difference(a, b):
    return [[x-y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def mm(a, b):
    return [[sum(a[i][t]*b[t][j] for t in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def quadratic(a, coeff):
    return sum(coeff[i]*a[i][j]*coeff[j] for i in range(len(coeff)) for j in range(len(coeff)))


def inverse2(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[F(a[1][1], det), F(-a[0][1], det)],
            [F(-a[1][0], det), F(a[0][0], det)]]


class YC5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()
        cls.records = json.loads((m.HERE/'YC5_ENDPOINTS.json').read_text())['endpoints']

    def test_packet_and_source_pins(self):
        self.assertEqual(self.result, json.loads((m.HERE/'YC5_RESULT.json').read_text()))
        self.assertTrue(all(self.result['checks'].values()))

    def test_fast_generator_matches_predecessor(self):
        for p in (m.y.potential(), m.c.shell((3, 1, 2))[-1][0],
                  m.y.mono(0, 0, 1, 4, 4, 5, 8)):
            self.assertEqual(m.generator(p), m.y.h0(p))

    def test_lagrange_projectors_at_all_nodes(self):
        nodes = (24, 32, 40, 48, 56)
        for target in nodes:
            coeff = m.lagrange_coefficients(nodes, target)
            for value in nodes:
                self.assertEqual(sum(c*value**i for i, c in enumerate(coeff)), int(value == target))

    def test_free_source_first_and_second_moments(self):
        for k in range(4):
            block = max(m.resolved_blocks(k).values(), key=lambda b: len(b['sources']))
            i = next(i for i, r in enumerate(block['sources']) if m.c.inner(r, r))
            source = block['sources'][i]
            derivative = m.y.h0(source)  # independent predecessor implementation
            self.assertEqual(m.c.inner(source, derivative), sum(e*w[i][i] for e, w in block['weights'].items()))
            self.assertEqual(m.c.inner(derivative, derivative), sum(e*e*w[i][i] for e, w in block['weights'].items()))

    def test_every_resolved_component_is_hidden(self):
        for k in range(4):
            for block in m.resolved_blocks(k).values():
                for parts in block['parts']:
                    for energy, p in parts.items():
                        self.assertGreaterEqual(energy, 3*k+24)
                        for vector in block['base']['vectors']:
                            self.assertEqual(m.c.inner(p, vector[0]), 0)

    def test_omitting_a_nonzero_source_component_breaks_reconstruction(self):
        block = m.resolved_blocks(0)[(0, 0, 0)]
        i = next(i for i, parts in enumerate(block['parts']) if parts)
        omitted_energy = next(iter(block['parts'][i]))
        kept = {}
        for energy, part in block['parts'][i].items():
            if energy != omitted_energy:
                kept = m.y.padd(kept, part)
        defect = m.y.padd(block['sources'][i], kept, -1)
        self.assertGreater(m.c.inner(defect, defect), 0)

    def test_resolver_rejects_a_false_hidden_floor(self):
        with self.assertRaises(AssertionError):
            m.resolve({m.y.ZERO: F(1)}, 24)

    def test_interaction_tensor_against_direct_polynomial_integral(self):
        for k in (0, 1):
            block = m.resolved_blocks(k)[tuple(int(i < k) for i in range(3))]
            coeff = [F((-1)**i, i+1) for i in range(len(block['sources']))]
            for z in (F(1, 3), F(2)):
                returned = {}
                for value, parts in zip(coeff, block['parts']):
                    for energy, p in parts.items():
                        returned = m.y.padd(returned, p, value/(energy-z))
                direct = m.c.inner(returned, m.y.pmul(m.y.potential(), returned))
                self.assertEqual(quadratic(m.interaction_return(k, z), coeff), direct)
                self.assertGreater(direct, 0)

    def test_interaction_tensor_against_dense_hidden_coordinates(self):
        for k in (0, 1):
            h, z = m.hidden_interaction(k), F(7, 3)
            v = [[value/(h['basis'][a][0]-z) for value in row] for a, row in enumerate(h['coords'])]
            dense = mm(list(map(list, zip(*v))), mm(h['potential'], v))
            self.assertEqual(m.interaction_return(k, z), dense)

    def test_hidden_source_span_is_not_assumed_closed_under_potential(self):
        for k in (0, 1):
            h = m.hidden_interaction(k)
            block = m.resolved_blocks(k)[tuple(int(i < k) for i in range(3))]
            image = m.y.pmul(m.y.potential(), h['basis'][-1][1])
            residual = image
            for p, norm, *_ in block['base']['vectors']:
                residual = m.y.padd(residual, p, -m.c.inner(image, p)/norm)
            for _, p, norm in h['basis']:
                residual = m.y.padd(residual, p, -m.c.inner(image, p)/norm)
            self.assertGreater(m.c.inner(residual, residual), 0)

    def test_tangent_and_secant_scalar_identities(self):
        for u in (F(0), F(1, 3), F(2), F(12)):
            s = u/3
            for t in (F(0), F(1, 5), F(2, 3), F(1)):
                x = t*u
                exact = 1/(1+x)
                tangent = (1+2*s-x)/(1+s)**2
                secant = 1-x/(1+u)
                self.assertEqual(exact-tangent, (x-s)**2/((1+x)*(1+s)**2))
                self.assertEqual(secant-exact, x*(u-x)/((1+u)*(1+x)))
                self.assertLessEqual(tangent, exact)
                self.assertLessEqual(exact, secant)

    def test_noncommuting_operator_return_bounds(self):
        d0, potential = [[F(2), 0], [0, F(5)]], [[F(2), F(1)], [F(1), F(2)]]
        self.assertNotEqual(mm(d0, potential), mm(potential, d0))
        free = inverse2(d0)
        interaction = mm(free, mm(potential, free))
        actual = inverse2([[4, 1], [1, 7]])
        lower = [[(3*free[i][j]-interaction[i][j])/4 for j in range(2)] for i in range(2)]
        upper = [[free[i][j]-interaction[i][j]/4 for j in range(2)] for i in range(2)]
        for matrix in (difference(actual, lower), difference(upper, actual),
                       difference(free, actual), difference(actual, inverse2([[8, 0], [0, 11]]))):
            self.assertEqual(m.c.inertia(matrix)[0], 0)
        wrong = [[3*value/4 for value in row] for row in free]  # dropping the actual interaction
        self.assertGreater(m.c.inertia(difference(actual, wrong))[0], 0)
        self.assertGreater(m.c.inertia(difference(actual, free))[0], 0)  # reversed inverse order

    def test_lower_return_refinement_dominates_yc4(self):
        for k in range(4):
            for block in m.resolved_blocks(k).values():
                old = m.c.pencil(block['base'], k, F(2), F(1), True)
                resolved = m.pencil(block, 2, 1, 'lower')
                dynamic = m.pencil(block, 2, 1, 'lower', True)
                self.assertEqual(m.c.inertia(difference(resolved, old))[0], 0)
                self.assertEqual(m.c.inertia(difference(dynamic, resolved))[0], 0)

    def test_upper_and_lower_forms_are_ordered(self):
        for k in range(4):
            for block in m.resolved_blocks(k).values():
                for theta, z in ((F(0), F(0)), (F(4), F(7)), (F(12), F(15))):
                    for dynamic in (False, True):
                        lower = m.pencil(block, theta, z, 'lower', dynamic)
                        upper = m.pencil(block, theta, z, 'upper', dynamic)
                        self.assertEqual(m.c.inertia(difference(upper, lower))[0], 0)

    def test_parameter_gates(self):
        for k in range(4):
            block = next(iter(m.resolved_blocks(k).values()))
            with self.assertRaises(ValueError):
                m.pencil(block, 1, 3*k+24, 'lower', True)
            with self.assertRaises(ValueError):
                m.pencil(block, -1, 0, 'lower', True)
            with self.assertRaises(ValueError):
                m.pencil(block, 1, 0, 'ritz', True)
        with self.assertRaises(ValueError):
            m.hidden_interaction(2)

    def test_bad_endpoint_and_missing_cells_rejected(self):
        record = copy.deepcopy(self.records[16])
        record['ground'][1][0] = record['ground'][1][1]
        with self.assertRaises(AssertionError):
            m.validate_endpoint(record)
        record = copy.deepcopy(self.records[-1])
        record['theta'] = '13'
        with self.assertRaises(ValueError):
            m.validate_endpoint(record)
        record = copy.deepcopy(self.records[0])
        record['second_lower'].pop()
        with self.assertRaises(ValueError):
            m.validate_endpoint(record)
        with self.assertRaises(AssertionError):
            m.validate_mesh(self.records[:-1])
        with self.assertRaises(AssertionError):
            m.validate_mesh(list(reversed(self.records)))

    def test_whole_window_and_first_sector_multiplicity(self):
        self.assertEqual(len(self.records), 105)
        self.assertEqual(len(self.result['coupling_cells']), 104)
        self.assertEqual(F(self.result['minimum_mesh_gap_lower']), F(103, 2000))
        for cell in self.result['coupling_cells']:
            self.assertGreaterEqual(F(cell['gap_lower']), F(1, 20))
            self.assertGreater(F(cell['other_centre_margin']), 0)
            self.assertGreater(F(cell['second_level_margin']), 0)
        self.assertEqual(self.result['first_excitation_multiplicity'], 3)

    def test_tighter_point_enclosures_without_larger_retained_cut(self):
        old = json.loads((m.ROOT/'physics/yc4/YC4_RESULT.json').read_text())
        for theta in ('2', '4', '6', '8'):
            lower, upper = map(F, self.result['point_gap_enclosures'][theta])
            old_lower, old_upper = map(F, old['point_gap_enclosures'][theta])
            self.assertGreater(lower, old_lower)
            self.assertLess(upper, old_upper)
        self.assertEqual([s['retained_rank'] for s in self.result['sectors']], [13, 21, 32, 23])
        self.assertEqual([s['source_span_rank'] for s in self.result['hidden_interaction_representations']], [27, 33])

    def test_no_physical_or_closure_upgrade(self):
        boundary = self.result['claim_boundary']
        for key in ('one_site_only', 'all_gauge_and_centre_sectors_covered',
                    'uncomputed_complement_floor_proved', 'full_source_resolved'):
            self.assertTrue(boundary[key])
        for key in ('hidden_interaction_span_assumed_invariant', 'retained_cut_enlarged_relative_to_yc4',
                    'floating_point_used_to_certify', 'weak_limit_absolute_splitting_proved',
                    'volume_uniform_or_4d_gap', 'formal_proof_assistant_verification'):
            self.assertFalse(boundary[key])


if __name__ == '__main__':
    unittest.main()
