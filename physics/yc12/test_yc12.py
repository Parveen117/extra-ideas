"""Independent completeness, source, hidden-interaction and geometry controls."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
import unittest
import sympy as sp
import yc12_correlated_cube_gap as y


class YC12Tests(unittest.TestCase):
    def test_every_small_leafless_support_is_a_face(self):
        edges, faces, _, _ = y.cube()
        leafless = []
        for size in range(1, 6):
            for subset in combinations(range(12), size):
                degree = Counter(v for e in subset for v in edges[e])
                if min(degree.values()) >= 2:
                    leafless.append(frozenset(subset))
        self.assertEqual(set(leafless), set(faces))
        self.assertEqual(len(leafless), 6)

    def test_hidden_skew_cycles_are_not_dropped(self):
        _, faces, adjacent, _ = y.cube()
        boundaries = {faces[p] ^ faces[q] for p, q in adjacent}
        all_six = set(y.cycle_supports(6))
        self.assertEqual(len(boundaries), 12)
        self.assertTrue(boundaries < all_six)
        self.assertEqual(len(all_six-boundaries), 4)

    def test_per_link_harmonics_distinguish_all_equal_energy_sources(self):
        _, faces, adjacent, opposite = y.cube()
        signatures = {
            18: [tuple(int(e in faces[p] ^ faces[q]) for e in range(12)) for p, q in adjacent],
            26: [tuple(int(e in faces[p])+int(e in faces[q]) for e in range(12)) for p, q in adjacent],
            24: [tuple(int(e in faces[p] | faces[q]) for e in range(12)) for p, q in opposite],
            32: [tuple(2*int(e in face) for e in range(12)) for face in faces]}
        for energy, rows in signatures.items():
            self.assertEqual(len(set(rows)), len(rows))
            self.assertTrue(all(sum(n*(n+2) for n in row) == energy for row in rows))

    def test_direct_low_face_moments_reconstruct_full_source_gram(self):
        # Proper subsets of the six faces are independent Haar traces.
        # E[W^2]=1/4, E[W^4]=1/8. Here at most four faces occur.
        moment = {0: Q(1), 1: Q(0), 2: Q(1, 4), 3: Q(0), 4: Q(1, 8)}
        def mean(indices):
            value = Q(1)
            for count in Counter(indices).values():
                value *= moment[count]
            return value
        grams = y.source_grams()
        for p in range(6):
            for q in range(6):
                # r_p=-2 sum_a W_p W_a+1/2; all face projections vanish.
                direct = 4*sum(mean((p, a, q, b)) for a in range(6) for b in range(6))-Q(1, 4)
                self.assertEqual(direct, sum(m[p+1][q+1] for m in grams.values()))

    def test_link_parity_is_a_real_cube_transformation(self):
        mask = y.centre_flip()
        for face in y.cube()[1]:
            self.assertEqual(sum((mask >> e) & 1 for e in face) % 2, 1)
        # A single face flip is impossible on this closed surface.
        faces = y.cube()[1]
        for candidate in range(1 << 12):
            signs = [sum((candidate >> e) & 1 for e in face) % 2 for face in faces]
            self.assertEqual(sum(signs) % 2, 0)

    def test_parity_envelope_keeps_hidden_interaction_outside_source(self):
        D = sp.diag(19, 21, 26)
        T = sp.Matrix([[0, 0, 2], [0, 0, 3], [2, 3, 0]])
        source = sp.Matrix([[1, 0], [0, 1], [0, 0]])
        J = sp.diag(1, 1, -1)
        self.assertEqual(J*T*J, -T)
        self.assertNotEqual(T*source, sp.zeros(3, 2))
        for theta in (sp.Rational(1, 2), sp.Rational(1)):
            for zeta in (-1, 0, 10):
                A = D-zeta*sp.eye(3)
                bare = source.T*A.inv()*source
                actual = source.T*(A-theta*T).inv()*source
                upper = bare/(1-(6*theta/(18-zeta))**2)
                for difference in (actual-bare, upper-actual):
                    matrix = [[Q(x) for x in row] for row in difference.tolist()]
                    self.assertEqual(y.y4.inertia(matrix)[0], 0)
                self.assertGreater(actual[0, 0], bare[0, 0])

    def test_subtracted_return_has_correct_schur_sign(self):
        A, D = sp.Matrix([[4, 1], [1, 7]]), sp.diag(19, 22)
        r = sp.Matrix([[1, 2], [0, 3]])
        z = sp.Rational(2)
        H = A.row_join(r.T).col_join(r.row_join(D))
        S = A-z*sp.eye(2)-r.T*(D-z*sp.eye(2)).inv()*r
        self.assertEqual((H-z*sp.eye(4)).det(), (D-z*sp.eye(2)).det()*S.det())
        wrong = A-z*sp.eye(2)+r.T*(D-z*sp.eye(2)).inv()*r
        self.assertNotEqual((H-z*sp.eye(4)).det(), (D-z*sp.eye(2)).det()*wrong.det())

    def test_inverse_requires_a_full_hidden_reserve(self):
        for theta, zeta in ((1, 12), (1, 18), (-1, 0)):
            with self.assertRaises(ValueError):
                y.return_matrix(theta, zeta, True)
        self.assertEqual(y.return_matrix(0, 0, False), y.return_matrix(0, 0, True))

    def test_spectral_enclosures_are_supported_by_exact_inertias(self):
        for theta in (Q(1, 4), Q(1)):
            for index in (0, 1):
                lo, hi = y.spectral_enclosure(theta, index)
                self.assertLessEqual(y.y4.inertia(y.pencil(theta, lo, True))[0], index)
                self.assertGreater(y.y4.inertia(y.pencil(theta, hi, False))[0], index)
                self.assertLessEqual(lo, hi)

    def test_full_cube_bound_is_monotone_over_the_declared_window(self):
        t = sp.symbols('t')
        k, h = sp.Rational(y.GRAD_F2), sp.Rational(y.HESS_F2)
        G = 2-sp.Rational(4, 3)*t-2*h*t**2-8*k*t**3-12*k**2*t**4
        self.assertTrue(all(c < 0 for c in sp.Poly(sp.diff(G, t), t).all_coeffs()))
        self.assertGreater(G.subs(t, 1), sp.Rational(1, 4))
        self.assertLess(G.subs(t, 1), 2-sp.Rational(4, 3)-2*h)
        for tbad in (-1, Q(1001, 1000)):
            with self.assertRaises(ValueError):
                y.cube_bounds(tbad)

    def test_boundary_charge_gap_cannot_be_replaced_by_singlet_gap(self):
        # One fundamental link is a valid full-factor excitation of energy 3.
        # It is charged at both endpoints; four links are needed for a singlet loop.
        self.assertEqual(1*(1+2), 3)
        self.assertEqual(4*1*(1+2), 12)
        self.assertEqual(y.join_bounds(0)['s'], 1)
        self.assertEqual(y.join_bounds(0)['g'], Q(1, 4))

    def test_tiling_matches_independent_ordinary_lattice_faces(self):
        shape = (4, 4, 6)
        t = y.tiling(shape)
        ordinary = y.y8.geometry(shape)
        direct = []
        for face in ordinary['faces']:
            direct.append(frozenset(t['edge_factors'][ordinary['edges'][e]] for e in face))
        self.assertEqual(Counter(direct), Counter(t['internal']+t['external']))
        self.assertEqual(set(t['edge_factors']), set(ordinary['edges']))

    def test_odd_collapsed_and_noninteger_tori_are_rejected(self):
        for shape in ((3, 4, 4), (2, 4, 4), (4, 4, 5), (4, 4), (True, 4, 4), (4.0, 4, 4)):
            with self.assertRaises(ValueError):
                y.tiling(shape)

    def test_interface_budget_counts_cube_incidence(self):
        t = y.tiling((4, 4, 4))
        self.assertEqual({n for f, n in t['incidence'].items() if f[0] == 'cube'}, {24})
        self.assertEqual({n for f, n in t['incidence'].items() if f[0] == 'bridge'}, {4})
        b = y.join_bounds(Q(1, 16384))
        self.assertLess(b['mapping'], b['r'])
        self.assertLess(b['contraction'], 1)
        self.assertEqual(b['actual_full_gap_lower'], Q(55, 256))
        with self.assertRaises(ValueError):
            y.join_bounds(Q(1, 16000))


if __name__ == '__main__':
    unittest.main()
