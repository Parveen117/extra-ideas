"""Independent moment checks, sign replay, controls and frozen DR2 evidence."""
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

import dr2_all_turn_gap as d2


def cart_add(p, q, a=1):
    out = dict(p)
    for e, x in q.items():
        out[e] = out.get(e, F(0))+a*x
    return {e: x for e, x in out.items() if x}


def cart_mul(p, q):
    out = {}
    for a, x in p.items():
        for b, y in q.items():
            e = tuple(i+j for i, j in zip(a, b))
            out[e] = out.get(e, F(0))+x*y
    return {e: x for e, x in out.items() if x}


def cart_derivative(p, axis):
    out = {}
    for e, x in p.items():
        if e[axis]:
            k = list(e)
            k[axis] -= 1
            out[tuple(k)] = x*e[axis]
    return out


def wick(p):
    """Independent Cartesian Gaussian moments, each variance = 1/2."""
    total = F(0)
    for e, x in p.items():
        if any(k % 2 for k in e):
            continue
        value = x
        for k in e:
            for i in range(1, k, 2):
                value *= F(i, 2)
        total += value
    return total


def cartesian_trial(d):
    dimension = 3*d

    def coordinate(i):
        e = [0]*dimension
        e[i] = 1
        return {tuple(e): F(1)}

    axes = [coordinate(i) for i in range(dimension)]
    q = {}
    for a in range(3):
        q = cart_add(q, cart_mul(axes[a], axes[a]))
        q = cart_add(q, cart_mul(axes[3+a], axes[3+a]), -1)
    potential = {}
    for i in range(d):
        for j in range(i):
            for a in range(3):
                for b in range(a):
                    minor = cart_add(cart_mul(axes[3*i+a], axes[3*j+b]),
                                     cart_mul(axes[3*i+b], axes[3*j+a]), -1)
                    potential = cart_add(potential, cart_mul(minor, minor))
    norm = wick(cart_mul(q, q))
    kinetic = F(0)
    for axis, x in enumerate(axes):
        derivative = cart_add(cart_derivative(q, axis), cart_mul(x, q), -1)
        kinetic += wick(cart_mul(derivative, derivative))/2
    energy = kinetic+wick(cart_mul(potential, cart_mul(q, q)))/F(2*(d-1))
    return norm, kinetic/norm, energy/norm, q


class DR2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = d2.run()

    def test_certificate_and_source_pins(self):
        self.assertTrue(self.result['all_pass'])
        self.assertTrue(all(self.result['checks'].values()))
        frozen = json.loads((Path(__file__).parent/'DR2_RESULT.json').read_text())
        self.assertEqual(frozen, dict(self.result, source_sha256=d2.source_pins()))

    def test_symbolic_generator_and_moments_against_cr1(self):
        cr = d2.dr.cr1
        exponents = [(a, b, c) for a in range(5) for b in range(3) for c in range(2)
                     if a+2*b+3*c <= 6]
        for d in (3, 4, 7, 23, 101):
            free = cr.free_means(d, F(1))
            for e in exponents:
                self.assertEqual({k: d2.evaluate(v, d) for k, v in d2.gen_mono(e).items()}, cr.gen_mono(e, d))
                self.assertEqual(d2.evaluate(d2.mono_mean(e), d), free({e: F(1)}))

    def test_symbolic_source_action_against_cr1(self):
        cr = d2.dr.cr1
        source = d2.source_numerator()
        for d in (3, 4, 11, 101):
            n = d-1
            p = {e: d2.evaluate(v, d) for e, v in source.items()}
            hp = cr.padd({}, cr.gen(p, d), -n)
            hp = cr.padd(hp, {e: n*(2*cr.degree(e)+F(3*d, 2))*v for e, v in p.items()})
            hp = cr.padd(hp, cr.pmul(p, {d2.E1: -F(n, 2), d2.E2: F(1, 2)}))
            actual = {e: d2.evaluate(v, d) for e, v in d2.scaled_action(source).items()}
            self.assertEqual(actual, hp)

    def test_turning_trial_by_cartesian_wick_and_exact_swap(self):
        for d in (3, 4, 6):
            norm, kinetic, energy, q = cartesian_trial(d)
            self.assertEqual(norm, 3)
            self.assertEqual(kinetic, F(3*d, 4)+1)
            self.assertEqual(energy, d2.excitation_upper_scaled(d))
            swapped = {e[3:6]+e[:3]+e[6:]: v for e, v in q.items()}
            self.assertEqual(swapped, {e: -v for e, v in q.items()})

    def test_nodal_ode_and_wrong_ground_interpretation_control(self):
        for m, j in ((3, 1), (4, 2), (11, 7), (101, 19)):
            cs = d2.laguerre(m, j)
            alpha = F(2*m-4, 3)
            self.assertEqual(d2.laguerre_residual(cs, alpha, j), ())
            damaged = list(cs)
            damaged[-1] += F(1, 10**9)
            self.assertNotEqual(d2.laguerre_residual(damaged, alpha, j), ())
        # For A_3, s*exp(-s) gives a ground upper 5/2. The nodal j=1
        # floor exceeds it, so treating a nodal test as a ground floor is false.
        self.assertGreater(d2.radial_floor_cube(3, 1), F(5, 2)**3)
        self.assertLess(d2.radial_floor_cube(3, 0), F(5, 2)**3)

    def test_outward_roots_and_uniform_margin(self):
        for x in (F(0), F(8), F(1, 7), F(999999999), F(2, 10**50)):
            lo, hi = d2.root_floor(x), d2.root_upper(x)
            self.assertLessEqual(lo**3, x)
            self.assertGreaterEqual(hi**3, x)
            self.assertLessEqual(hi-lo, F(1, 10**12))
        for d in (11, 12, 100, 10**9):
            self.assertGreater(d2.gap_coefficient(d), F(1, 64))
        self.assertGreater(d2.gap_coefficient(9), 0)
        # The shifted polynomial identities in run(), not these samples,
        # certify the unbounded dimension ranges.

    def test_actual_source_is_needed_for_this_count_margin(self):
        for d in (9, 11, 100, 10**6):
            count_minus_bare_gaussian = -F(17, 32*(d-1))
            self.assertLess(count_minus_bare_gaussian, 0)
            self.assertEqual(count_minus_bare_gaussian+d2.gain(d), d2.gap_coefficient(d))
            self.assertGreater(d2.gap_coefficient(d), 0)

    def test_invalid_domains_are_rejected(self):
        for d in (2, -1, F(7, 2), True):
            with self.assertRaises(ValueError):
                d2.gain(d)
        for m, j in ((2, 0), (3, -1), (F(7, 2), 1), (3, F(1, 2))):
            with self.assertRaises(ValueError):
                d2.radial_floor_cube(m, j)
        with self.assertRaises(ValueError):
            d2.gap_coefficient(8)
        with self.assertRaises(ValueError):
            d2.small_count_scaled(11, d2.frozen_floors())
        with self.assertRaises(ValueError):
            d2.root_floor(-1)


if __name__ == '__main__':
    unittest.main()
