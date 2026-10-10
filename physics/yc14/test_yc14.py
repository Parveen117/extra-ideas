"""Controls for transported metrics, source sectors and complete return bounds."""
from fractions import Fraction as Q
import unittest
import sympy as sp
import yc14_observer_reset_and_interface as y


class YC14Tests(unittest.TestCase):
    def test_metric_reset_preserves_norm_and_both_projections(self):
        family = y.reset_family()
        for u, t in ((Q(1), Q(2)), (Q(2), Q(0)), (Q(1, 3), Q(-2, 5))):
            values = {family['u']: sp.Rational(u), family['t']: sp.Rational(t)}
            C, G, cut = (family[key].subs(values) for key in ('reset', 'metric', 'cut'))
            for v in (sp.Matrix([2, -3]), sp.Matrix([1, 1])):
                plus, minus = (sp.eye(2)+cut)*v/2, (sp.eye(2)-cut)*v/2
                self.assertTrue(y.zero((C*v).dot(C*v)-(v.T*G*v)[0]))
                self.assertTrue(y.zero((plus.T*G*minus)[0]))
                self.assertGreater((v.T*G*v)[0], 0)

    def test_rotation_only_cannot_remove_nonzero_obliqueness(self):
        f = y.reset_family()
        cut = f['cut'].subs({f['u']: 2, f['t']: Q(3, 7)})
        for O in (sp.eye(2), y.rotation(Q(3, 5), Q(4, 5)), y.OBS):
            changed = O*cut*O.T
            self.assertNotEqual(changed, changed.T)

    def test_whitening_is_not_the_old_return_metric_whitening(self):
        f = y.reset_family()
        values = {f['u']: 2, f['t']: 0}
        C, G, B, cut = (f[key].subs(values) for key in ('reset', 'metric', 'cycle', 'cut'))
        self.assertEqual(C*cut*C.inv(), y.K)
        self.assertNotEqual(G*cut*G.inv(), y.K)
        self.assertEqual(G*G, B)
        self.assertEqual(C.det(), 1)
        self.assertNotEqual(C.T*C, sp.eye(2))

    def test_leaving_turn_untransported_falsely_erases_memory(self):
        f = y.reset_family()
        values = {f['u']: 2, f['t']: 0}
        C, B = (f[key].subs(values) for key in ('reset', 'cycle'))
        actual = (C*y.R*C.inv()*y.K)**2
        self.assertEqual(actual, C*B*C.inv())
        self.assertNotEqual(actual, (y.R*y.K)**2)
        self.assertEqual(actual.trace(), Q(257, 16))

    def test_flat_local_readings_can_retain_open_path_ratio_memory(self):
        frames, cuts, steps = y.memory_example()
        self.assertTrue(all(O*k*O.T == y.K for O, k in zip(frames, cuts)))
        ratios = [step[1, 0]/step[0, 0] for step in steps]
        self.assertEqual(ratios, [Q(4, 3), Q(16, 63)])
        self.assertEqual(steps[1]*steps[0], frames[-1])
        self.assertEqual(frames[-1].T*steps[1]*steps[0], sp.eye(2))

    def test_ground_transform_keeps_the_derivative_term(self):
        # A confining scalar carrier independently checks the exact algebra,
        # not a claim that its ground is a YM cube ground.
        x = sp.symbols('x', real=True)
        omega, f = sp.exp(-x*x), x**3+2*x+1
        H = lambda v: -sp.diff(v, x, 2)+(4*x*x-2)*v
        self.assertEqual(sp.simplify(H(omega)), 0)
        self.assertEqual(sp.simplify(H(omega*f)/omega+sp.diff(f, x, 2)-4*x*sp.diff(f, x)), 0)
        self.assertNotEqual(sp.simplify(H(omega*f)/omega), -sp.diff(f, x, 2))

    def test_source_vector_bound_uses_charged_cube_floor(self):
        two, four = y.source_cost(2), y.source_cost(4)
        self.assertEqual(two['source_floor'], Q(13, 2))
        self.assertEqual(two['inverse_vector_upper'], Q(1, 13))
        self.assertEqual(four['inverse_vector_upper'], Q(1, 24))
        self.assertEqual(two['weighted_seed_upper'], Q(4, 13))
        self.assertGreater(y.y12.cube_bounds(1)['actual_full_gap_lower'], y.GAP)
        self.assertNotEqual(two['source_floor'], 6+2*11)

    def test_inverse_is_not_division_by_source_mean(self):
        # Full positive operator with the same norm and first energy moment.
        # The returning vector leaves the one-dimensional source span.
        D, r = sp.Matrix([[12, 5], [5, 12]]), sp.Matrix([sp.Rational(1, 2), 0])
        value = (r.T*D.inv()*r)[0]
        self.assertEqual(r.dot(r), Q(1, 4))
        self.assertEqual((r.T*D*r)[0], 3)
        self.assertGreater(value, Q(1, 48))
        self.assertLess(value, Q(1, 26))
        self.assertNotEqual((D.inv()*r)[1], 0)
        self.assertLessEqual((D.inv()*r).dot(D.inv()*r), Q(1, 169))

    def test_bridge_signature_distinguishes_faces_on_periodic_seams(self):
        for shape in ((4, 4, 4), (4, 6, 8), (8, 4, 6)):
            row = y.geometry_certificate(shape)
            self.assertTrue(row['bridge_signatures_unique'])
            self.assertEqual(row['cube_profiles'], [(24, 0)])
            self.assertEqual(row['bridge_profiles'], [(2, 2)])

    def test_internal_centres_do_not_change_cube_wilson_potential(self):
        edges, faces, _, _ = y.y12.cube()
        for vertex in {v for edge in edges for v in edge}:
            affected = {i for i, edge in enumerate(edges) if vertex in edge}
            self.assertEqual(len(affected), 3)
            self.assertTrue(all(len(face & affected) in (0, 2) for face in faces))
            self.assertTrue(all(i in affected for i, edge in enumerate(edges) if edge[0] == vertex))

    def test_whole_window_not_just_endpoint_has_all_three_reserves(self):
        for i in range(9):
            row = y.join_bounds(y.ETA_MAX*i/8)
            self.assertGreater(row['mapping_reserve'], 0)
            self.assertLess(row['contraction'], 1)
            self.assertLess(row['relative_return'], 1)
            self.assertGreaterEqual(row['full_gap_lower'], Q(21, 100))
            self.assertEqual(row['full_gap_lower'], Q(1, 4)-Q(1024, 5)*row['eta'])
            self.assertEqual(row['nonlinear_minimum_support'], 1)

    def test_new_window_requires_exact_seed_not_coarse_map(self):
        row = y.join_bounds(y.ETA_MAX)
        coarse = 4*row['beta']/y.GAP*(1+2*y.RADIUS)*y.EXP_UPPER
        self.assertGreater(coarse, y.RADIUS)
        self.assertEqual(row['mapping_reserve'], Q(7, 166400))
        self.assertEqual(row['fixed_point_norm_upper'], Q(15, 1976))
        self.assertLess(row['fixed_point_norm_upper'], y.RADIUS)

    def test_old_endpoint_improves_without_changing_internal_cube_window(self):
        old = y.y12.join_bounds(Q(1, 16384))
        new = y.join_bounds(Q(1, 16384))
        self.assertEqual(new['full_gap_lower'], Q(19, 80))
        self.assertGreater(new['full_gap_lower'], Q(55, 256))
        self.assertEqual(new['factor_floor'], old['g'])
        with self.assertRaises(ValueError):
            y.y12.cube_bounds(Q(3, 2))

    def test_outside_proved_carriers_and_windows_rejected(self):
        for eta in (-1, y.ETA_MAX+Q(1, 10**8)):
            with self.assertRaises(ValueError):
                y.join_bounds(eta)
        for shape in ((2, 4, 4), (3, 4, 4), (4, 4)):
            with self.assertRaises(ValueError):
                y.geometry_certificate(shape)
        for bridges, gap in ((0, y.GAP), (3, y.GAP), (2, 0), (2, 11)):
            with self.assertRaises(ValueError):
                y.source_cost(bridges, gap)

    def test_frozen_input_pins_still_match(self):
        for path in ('physics/yc12/YC12_RESULT.json', 'physics/yc13/YC13_RESULT.json'):
            y.verify_predecessor(path)


if __name__ == '__main__':
    unittest.main()
