"""Focused controls for the spatial carrier and full hidden-space comparison."""
from fractions import Fraction as F
import copy
import json
import unittest
import sympy as sp
import yc7_two_cell_return as y


class TwoCellTests(unittest.TestCase):
    def test_haar_coordinates_and_sphere_constraint(self):
        self.assertEqual(y.inner(y.mono(0), y.mono(0)), F(1, 4))
        self.assertEqual(y.inner(y.mono(0, 1), y.mono(0, 1)), F(1, 24))
        self.assertEqual(y.inner(y.mono(0, 4), y.mono(0, 4)), F(1, 16))
        norm = y.add(y.sum_polys(y.mono(i, i) for i in range(4)), y.ONE, -1)
        self.assertEqual(y.inner(norm, norm), 0)
        with self.assertRaises(ValueError):
            y.moment((0,)*12)

    def test_plaquettes_against_direct_su2_matrices_and_finite_local_gauge(self):
        def matrix(q):
            a, b, c, d = map(sp.Rational, q)
            # Quaternion convention: e_i = -i sigma_i.
            return sp.Matrix([[a-sp.I*d, -c-sp.I*b], [c-sp.I*b, a+sp.I*d]])
        qs = [(F(3, 5), F(4, 5), 0, 0), (F(1, 2),)*4,
              (0, F(3, 5), F(4, 5), 0), (F(3, 5), 0, 0, F(4, 5)),
              (F(1, 2), F(-1, 2), F(1, 2), F(-1, 2)), (0, 0, F(3, 5), F(4, 5))]
        matrices = [matrix(q) for q in qs]

        def direct(ms):
            out = []
            for a, b, c in ((0, 2, 4), (1, 4, 2), (0, 3, 5), (1, 5, 3)):
                out.append(sp.simplify(1-sp.trace(ms[a]*ms[c]*ms[a].H*ms[b].H)/2))
            for b, c in ((2, 3), (4, 5)):
                out.append(sp.simplify(1-sp.trace(ms[b]*ms[c]*ms[b].H*ms[c].H)/2))
            return out

        values = tuple(a for q in qs for a in q)
        from math import prod
        polynomial = [sum(c*prod(a**e for a, e in zip(values, power)) for power, c in p.items())
                      for p in y.plaquettes()]
        self.assertEqual(direct(matrices), polynomial)
        gs = [matrix((F(1, 2),)*4), matrix((F(3, 5), 0, F(4, 5), 0))]
        transformed = [gs[a]*m*gs[b].H for m, (a, b) in zip(matrices, y.EDGES)]
        self.assertEqual(direct(transformed), polynomial)
        self.assertTrue(all(0 <= v <= 2 for v in polynomial))

    def test_global_centre_is_not_independent_link_evenness(self):
        v = y.sum_polys(y.plaquettes())
        changed = {e: c*(-1)**sum(e[8:12]) for e, c in v.items()}
        self.assertGreater(y.inner(y.add(changed, v, -1), y.add(changed, v, -1)), 0)
        self.assertTrue(y.centre_even(v))
        self.assertTrue(y.centre_even(y.basis()[1]))
        self.assertTrue(any(sum(e[8:12]) % 2 for e in y.basis()[1]))

    def test_complete_cut_contains_two_winding_products(self):
        shells = y.electric_shells()
        self.assertEqual(sum(n for _, e, n in shells if e == 6), 2)
        self.assertEqual(sum(n for _, e, n in shells if e == 8), 4)
        self.assertEqual(sum(n for _, e, n in shells), 7)
        self.assertEqual(y.singlet_multiplicity([1, 1, 1, 1]), 2)
        self.assertEqual(y.singlet_multiplicity([1, 0]), 0)

    def test_hidden_source_resolution_includes_interface_kinetic_cost(self):
        d = y.data()
        self.assertEqual({l: m[0][0] for l, m in d['weights'].items() if m[0][0]},
                         {14: F(3, 4), 16: F(7, 24)})
        self.assertEqual(d['M'][0][0], F(11, 2))
        self.assertEqual(d['B'][0][0], F(25, 24))

    def test_cross_source_cannot_be_dropped(self):
        d = y.data()
        self.assertEqual(d['B_cross'][1][3], F(1, 4))
        separate = y.inner(d['source_local'][1], d['source_local'][3])
        separate += y.inner(d['source_interface'][1], d['source_interface'][3])
        self.assertNotEqual(separate, d['B'][1][3])
        self.assertEqual(separate+d['B_cross'][1][3], d['B'][1][3])

    def test_ordered_schur_forms(self):
        for theta, z in ((F(1, 2), F(5)), (F(2), F(11))):
            lo, hi, ritz = [sp.Matrix(y.pencil(theta, z, method)) for method in ('lower', 'upper', 'ritz')]
            self.assertEqual(y.inertia((hi-lo).tolist())[0], 0)
            self.assertEqual(y.inertia((ritz-hi).tolist())[0], 0)

    def test_noncommuting_full_resolvent_sandwich(self):
        d0 = sp.diag(12, 17)
        w = sp.Matrix([[3, 1], [1, 2]])
        theta, z = sp.Rational(2), sp.Rational(11)
        source = sp.Matrix([[1, 2], [-1, 3]])
        exact = source.T*(d0+theta*w-z*sp.eye(2)).inv()*source
        upper = source.T*(d0-z*sp.eye(2)).inv()*source
        lower = source.T*(d0+(12*theta-z)*sp.eye(2)).inv()*source
        self.assertNotEqual(d0*w, w*d0)
        self.assertEqual(y.inertia((upper-exact).tolist())[0], 0)
        self.assertEqual(y.inertia((exact-lower).tolist())[0], 0)

    def test_spatial_join_changes_source_as_well_as_hidden_inverse(self):
        theta, h, z = sp.Rational(3, 2), sp.Rational(2, 5), sp.Rational(1)
        d = sp.Matrix([[15, 1], [1, 19]])
        w = sp.Matrix([[2, 1], [1, 3]])
        r = sp.Matrix([[1, 2], [3, -1]])
        s = sp.Matrix([[2, -1], [0, 1]])
        old = (d-z*sp.eye(2)).inv()
        new = (d+theta*h*w-z*sp.eye(2)).inv()
        exact = (r+h*s).T*new*(r+h*s)-r.T*old*r
        rhs = h*(s.T*new*r+r.T*new*s)+h*h*s.T*new*s-theta*h*r.T*new*w*old*r
        self.assertEqual(exact, rhs)
        self.assertNotEqual(exact, -theta*h*r.T*new*w*old*r)

    def test_signed_frames_keep_the_full_return(self):
        r, d = sp.Matrix([[1, 2], [3, -1]]), sp.Matrix([[15, 1], [1, 19]])
        up, uq = sp.Matrix([[0, -1], [-1, 0]]), sp.diag(-1, 1)
        sigma = r.T*(d-sp.eye(2)).inv()*r
        rp, dp = uq*r*up.T, uq*d*uq.T
        self.assertEqual(rp.T*(dp-sp.eye(2)).inv()*rp, up*sigma*up.T)

    def test_endpoint_enclosures_and_point_gap(self):
        rows = json.loads((y.HERE/'YC7_ENDPOINTS.json').read_text())['endpoints']
        for row in rows:
            y.validate_endpoint(row)
        last = rows[-1]
        self.assertEqual(F(last['second'][0])-F(last['ground'][1]), F(12753, 10000))

    def test_whole_interval_uses_ordered_energies_not_gap_monotonicity(self):
        rows = json.loads((y.HERE/'YC7_ENDPOINTS.json').read_text())['endpoints']
        reserves = [F(a['second'][0])-F(b['ground'][1]) for a, b in zip(rows, rows[1:])]
        self.assertEqual(min(reserves), F(12431, 10000))
        self.assertGreater(min(reserves), F(6, 5))

    def test_invalid_thresholds_and_false_endpoint_rejected(self):
        for args in ((-1, 0, 'lower'), (1, 12, 'lower'), (1, 13, 'upper'), (1, 0, 'bad')):
            with self.assertRaises(ValueError):
                y.pencil(*args)
        row = copy.deepcopy(json.loads((y.HERE/'YC7_ENDPOINTS.json').read_text())['endpoints'][0])
        row['second'] = ['7', '8']
        with self.assertRaises(AssertionError):
            y.validate_endpoint(row)


if __name__ == '__main__':
    unittest.main()
