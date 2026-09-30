"""Finite exact checks supporting, but not replacing, the R1 written proofs."""

import unittest
from fractions import Fraction as Q

from response_transport import Edge, ResponseNetwork, triangle


def multiply(left, right):
    """Independent dense matrix multiplication for the operator checks."""
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


class ReturnInvariantTests(unittest.TestCase):
    def setUp(self):
        self.net = triangle()
        self.cycle = ('A', 'B', 'C', 'A')
        self.tree = ('ab', 'bc')

    def test_exact_return_reversal_and_path_comparison(self):
        self.assertEqual(self.net.holonomy(self.cycle), Q(6, 7))
        self.assertEqual(self.net.holonomy(tuple(reversed(self.cycle))), Q(7, 6))
        self.assertEqual(self.net.transport(('A', 'B', 'C')) /
                         self.net.transport(('A', 'C')), Q(6, 7))
        self.assertEqual(self.net.transport(('A',)), Q(1))

    def test_vertex_rescaling_cancels_exactly(self):
        changed = self.net.rescale({'A': Q(5), 'B': Q(11), 'C': Q(13)})
        self.assertEqual([e.gain for e in changed.edges],
                         [Q(22, 5), Q(39, 11), Q(5, 91)])
        self.assertEqual(changed.holonomy(self.cycle), Q(6, 7))
        self.assertNotEqual(changed.transport(('A', 'B', 'C')),
                            self.net.transport(('A', 'B', 'C')))

    def test_tree_normalization_and_exact_reconstruction(self):
        normal, gauges, invariants = self.net.normal_form('A', self.tree)
        self.assertEqual([e.gain for e in normal.edges], [Q(1), Q(1), Q(6, 7)])
        self.assertEqual(invariants, {'ca': Q(6, 7)})
        recovered = normal.rescale({v: 1 / g for v, g in gauges.items()})
        self.assertEqual(recovered.edges, self.net.edges)

    def test_reference_orbit_has_identical_normal_form(self):
        for gauges in ({'A': 1, 'B': -2, 'C': 3},
                       {'A': Q(2, 3), 'B': Q(7, 5), 'C': Q(-11, 13)}):
            with self.subTest(gauges=gauges):
                changed = self.net.rescale(gauges)
                self.assertEqual(self.net.normal_form('A', self.tree)[0].edges,
                                 changed.normal_form('A', self.tree)[0].edges)

    def test_changing_tree_preserves_every_tested_return(self):
        for tree in (('ab', 'bc'), ('bc', 'ca'), ('ca', 'ab')):
            with self.subTest(tree=tree):
                normalized, _, chords = self.net.normal_form('A', tree)
                self.assertEqual(len(chords), self.net.cycle_rank)
                self.assertEqual(normalized.holonomy(self.cycle), Q(6, 7))

    def test_two_independent_cycles(self):
        net = ResponseNetwork(('A', 'B', 'C', 'D'), (
            Edge('ab', 'A', 'B', 2), Edge('bc', 'B', 'C', 3),
            Edge('cd', 'C', 'D', 5), Edge('ca', 'C', 'A', Q(1, 7)),
            Edge('da', 'D', 'A', Q(1, 11)),
        ))
        normal, _, values = net.normal_form('A', ('ab', 'bc', 'cd'))
        self.assertEqual(net.cycle_rank, 2)
        self.assertEqual(values, {'ca': Q(6, 7), 'da': Q(30, 11)})
        changed = net.rescale({'A': 3, 'B': 7, 'C': 11, 'D': 13})
        self.assertEqual(normal.edges, changed.normal_form('A', ('ab', 'bc', 'cd'))[0].edges)

    def test_exact_common_parameter_responses_are_flat(self):
        tangent = {'A': Q(2), 'B': Q(5), 'C': Q(11)}
        flat = ResponseNetwork(self.net.vertices, (
            Edge(e.name, e.tail, e.head, tangent[e.head] / tangent[e.tail])
            for e in self.net.edges
        ))
        self.assertEqual(flat.holonomy(self.cycle), Q(1))
        self.assertEqual(flat.normal_form('A', self.tree)[2], {'ca': Q(1)})

    def test_all_nonzero_return_values_are_admissible(self):
        for value in (Q(-1), Q(1, 3), Q(1), Q(6, 7), Q(5)):
            with self.subTest(value=value):
                net = ResponseNetwork(('A', 'B', 'C'), (
                    Edge('ab', 'A', 'B', 1), Edge('bc', 'B', 'C', 1),
                    Edge('ca', 'C', 'A', value),
                ))
                self.assertEqual(net.holonomy(self.cycle), value)
                self.assertEqual(net.normal_form('A', self.tree)[2], {'ca': value})

    def test_six_stage_source_chain_has_no_return_parameter(self):
        vertices = ('OM', 'NA', 'MA', 'SHI', 'VA', 'YA')
        weights = (2, 3, 5, 7, 11)
        net = ResponseNetwork(vertices, (
            Edge(str(i), vertices[i], vertices[i + 1], w)
            for i, w in enumerate(weights)
        ))
        normal, _, invariants = net.normal_form('OM')
        self.assertEqual(net.cycle_rank, 0)
        self.assertEqual(invariants, {})
        self.assertTrue(all(e.gain == 1 for e in normal.edges))

    def test_typed_operator_composition_realizes_the_return(self):
        ab = self.net.path_operator(('A', 'B'))
        bc = self.net.path_operator(('B', 'C'))
        ca = self.net.path_operator(('C', 'A'))
        product = multiply(ca, multiply(bc, ab))
        self.assertEqual(product, [[Q(6, 7), Q(0), Q(0)],
                                   [Q(0), Q(0), Q(0)],
                                   [Q(0), Q(0), Q(0)]])
        self.assertEqual(product, self.net.path_operator(self.cycle))

    def test_operator_reference_change_is_diagonal_conjugation(self):
        gauges = (Q(5), Q(11), Q(13))
        diag = [[gauges[i] if i == j else Q(0) for j in range(3)] for i in range(3)]
        inv = [[1 / gauges[i] if i == j else Q(0) for j in range(3)] for i in range(3)]
        changed = self.net.rescale(dict(zip(self.net.vertices, gauges)))
        for path in (('A', 'B'), ('A', 'B', 'C'), self.cycle):
            with self.subTest(path=path):
                expected = multiply(diag, multiply(self.net.path_operator(path), inv))
                self.assertEqual(expected, changed.path_operator(path))

    def test_constant_returns_require_zero_cycle_log_rate(self):
        rates = {'ab': Q(2), 'bc': Q(3), 'ca': Q(-5)}
        self.assertEqual(self.net.rate_potential(rates, 'A', self.tree),
                         {'A': Q(0), 'B': Q(2), 'C': Q(5)})
        self.assertIsNone(self.net.rate_potential({'ab': 2, 'bc': 3, 'ca': -4},
                                                'A', self.tree))
        physical_change = ResponseNetwork(self.net.vertices, (
            Edge('ab', 'A', 'B', 2), Edge('bc', 'B', 'C', 4),
            Edge('ca', 'C', 'A', Q(1, 7)),
        ))
        self.assertEqual(physical_change.holonomy(self.cycle), Q(8, 7))

    def test_nonreciprocal_or_undefined_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            Edge('zero', 'A', 'B', 0)
        with self.assertRaises(TypeError):
            Edge('approximate', 'A', 'B', 0.5)
        with self.assertRaises(ValueError):
            ResponseNetwork(('A', 'B', 'C'), (Edge('ab', 'A', 'B', 2),))
        with self.assertRaises(ValueError):
            self.net.rescale({'A': 0, 'B': 1, 'C': 1})
        with self.assertRaises(ValueError):
            self.net.normal_form('A', ('ab',))
        with self.assertRaises(ValueError):
            self.net.holonomy(('A', 'B'))


if __name__ == '__main__':
    unittest.main()
