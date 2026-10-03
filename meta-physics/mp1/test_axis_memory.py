"""Behavioral and adversarial checks for MP-1; rational arithmetic only."""

from fractions import Fraction as Q
from itertools import product
import unittest

import axis_memory as am


class AxisMemoryTests(unittest.TestCase):
    def test_existing_emk_relations(self):
        self.assertEqual(am.mul(am.K, am.K), am.I)
        self.assertEqual(am.mul(am.R, am.R), am.scale(am.I, -1))
        self.assertEqual(am.mul(am.K, am.R), am.scale(am.RK, -1))

    def test_quadratic_observer_from_actual_conjugation(self):
        for k in range(-13, 14):
            p = am.rational_turn(Q(k, 7))
            x, y = am.axis_coefficients(p)
            q = am.axis_operator(p)
            self.assertEqual(q, am.add(am.scale(am.K, x), am.scale(am.RK, y)))
            self.assertEqual(am.mul(q, q), am.I)
            self.assertEqual(am.axis_operator(tuple(-v for v in p)), q)

    def test_exactly_two_sheets_in_declared_rational_family(self):
        turns = [am.rational_turn(Q(k, 5)) for k in range(-8, 9)] + [(-1, 0)]
        turns += [tuple(-x for x in p) for p in turns]
        for a, b in product(turns, repeat=2):
            self.assertEqual(am.axis_operator(a) == am.axis_operator(b),
                             a == b or a == tuple(-x for x in b))

    def test_labeled_axis_does_not_erase_axis_sign(self):
        self.assertNotEqual(am.axis_operator((1, 0)), am.axis_operator((0, 1)))

    def test_reference_probe_recovers_relative_sign(self):
        for v in [(1, 0), (2, 3), (Q(2, 3), Q(-1, 7))]:
            self.assertEqual(am.reference_intensity((1, 0), v), 4*am.dot(v, v))
            self.assertEqual(am.reference_intensity((-1, 0), v), 0)

    def test_winding_repetitions_and_signs(self):
        for n in range(-12, 13):
            r = am.polygon_report(am.diamond(n))
            self.assertEqual(r['winding'], n)
            self.assertEqual(r['relative_frame_sign'], -1 if n % 2 else 1)

    def test_concatenated_opposite_loops_cancel(self):
        a, b = am.diamond(3), am.diamond(-2)
        self.assertEqual(am.winding(a[:-1] + b), 1)
        self.assertEqual(am.winding(a[:-1] + tuple(reversed(a))), 0)

    def test_subdivision_and_repeated_vertices(self):
        p = am.diamond(1)
        refined = []
        for a, b in zip(p, p[1:]):
            refined += [a, a, tuple((Q(x)+y)/2 for x, y in zip(a, b))]
        refined += [p[-1]]
        self.assertEqual(am.winding(refined), 1)

    def test_half_open_ray_vertex_and_horizontal_edges(self):
        boundary = ((2, 0), (2, 1), (-2, 1), (-2, -1), (2, -1), (2, 0))
        touching = ((1, 0), (2, 0), (2, 1), (1, 1), (1, 0))
        self.assertEqual(am.winding(boundary), 1)
        self.assertEqual(am.winding(tuple(reversed(boundary))), -1)
        self.assertEqual(am.winding(touching), 0)

    def test_orientation_and_offset_controls(self):
        p = am.diamond(1)
        shear = tuple((2*x+3*y, y) for x, y in p)
        reflected = tuple((-x, y) for x, y in p)
        shifted = tuple((x+2, y) for x, y in p)
        self.assertEqual(am.winding(shear), 1)
        self.assertEqual(am.winding(reflected), -1)
        self.assertEqual(am.winding(shifted), 0)

    def test_exact_segment_minima(self):
        self.assertEqual(am.segment_clearance_squared((1, 0), (0, 1)), Q(1, 2))
        self.assertEqual(am.segment_clearance_squared((-1, 1), (1, 1)), 1)
        self.assertEqual(am.segment_clearance_squared((2, 1), (3, 1)), 5)
        self.assertEqual(am.segment_clearance_squared((2, 1), (2, 1)), 5)

    def test_zero_and_open_paths_rejected(self):
        for p in [(), ((1, 0),), ((1, 0), (0, 1)),
                  ((1, 0), (-1, 0), (1, 0)), ((0, 0), (0, 0))]:
            with self.assertRaises(ValueError):
                am.winding(p)

    def test_no_float_rounding(self):
        for bad in [0.5, True, '1/2']:
            with self.assertRaises(TypeError):
                am.rational_turn(bad)
        with self.assertRaises(ValueError):
            am.axis_operator((1, 1))

    def test_all_6561_vertex_perturbations_in_declared_grid(self):
        p = am.diamond(1)
        r = am.polygon_report(p, vertex_error=Q(1, 4))
        self.assertTrue(r['tube_certified'])
        for noise in product((-Q(1, 4), Q(0), Q(1, 4)), repeat=8):
            q = tuple((Q(a)+noise[2*j], Q(b)+noise[2*j+1])
                      for j, (a, b) in enumerate(p[:-1]))
            self.assertEqual(am.winding(q+(q[0],)), 1)

    def test_strict_boundary_rejected_and_zero_reached(self):
        p = am.diamond(1)
        r = am.polygon_report(p, vertex_error=Q(1, 2))
        self.assertEqual(r['tube_squared_bound'], r['clearance_squared'])
        self.assertFalse(r['tube_certified'])
        q = tuple((x-Q(1, 2), y-Q(1, 2)) for x, y in p)
        with self.assertRaises(ValueError):
            am.winding(q)

    def test_interpolation_budget_is_not_ignored(self):
        self.assertTrue(am.polygon_report(am.diamond(1), vertex_error=Q(1, 8),
                                         interpolation_error=Q(1, 8))['tube_certified'])
        self.assertFalse(am.polygon_report(am.diamond(1), vertex_error=Q(1, 8),
                                          interpolation_error=Q(3, 8))['tube_certified'])
        with self.assertRaises(ValueError):
            am.polygon_report(am.diamond(1), vertex_error=-1)

    def test_omitting_between_sample_path_hides_a_loop(self):
        self.assertEqual(am.winding(am.diamond(0)), 0)
        self.assertEqual(am.winding(am.diamond(1)), 1)
        self.assertEqual(am.diamond(0)[::len(am.diamond(0))-1],
                         am.diamond(1)[::len(am.diamond(1))-1])

    def test_phase_sample_recovery_with_branch_bound(self):
        for direction in [-1, 1]:
            p = tuple(Q(3, 5) + direction*Q(j, 9) for j in range(31))
            wrapped = tuple((2*x) % 1 for x in p)
            self.assertEqual(am.principal_frame_lift(wrapped, p[0]), p)

    def test_sharp_quarter_turn_tie(self):
        plus, minus = Q(1, 4), -Q(1, 4)
        self.assertEqual((2*plus) % 1, (2*minus) % 1)
        with self.assertRaises(ValueError):
            am.principal_frame_lift((0, Q(1, 2)))

    def test_hidden_half_turn_without_bound(self):
        actual = (Q(0), Q(1, 2))
        observed = tuple((2*x) % 1 for x in actual)
        self.assertEqual(am.principal_frame_lift(observed), (0, 0))
        self.assertNotEqual(am.principal_frame_lift(observed), actual)

    def test_phase_input_contract(self):
        for p in [(), (-1, 0), (0, 1)]:
            with self.assertRaises(ValueError):
                am.principal_frame_lift(p)
        with self.assertRaises(ValueError):
            am.principal_frame_lift((0,), Q(1, 4))

    def test_parity_does_not_recover_integer_history(self):
        for n in range(10):
            self.assertEqual(am.polygon_report(am.diamond(n))['relative_frame_sign'],
                             am.polygon_report(am.diamond(n+2))['relative_frame_sign'])
            self.assertNotEqual(am.winding(am.diamond(n)), am.winding(am.diamond(n+2)))
        self.assertEqual(am.required_history_labels(4), 9)

    def test_flat_patch_transport_telescopes(self):
        frames = [am.rotor(am.rational_turn(Q(k, 9))) for k in [0, 1, 3, 2, 0]]
        out = am.I
        for a, b in zip(frames, frames[1:]):
            out = am.mul(am.mul(b, am.inverse(a)), out)
        self.assertEqual(out, am.I)


if __name__ == '__main__':
    unittest.main()
