"""CM2 independent local calculus, sector and complete-error controls."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import unittest

import cm2_compact_core_matching as m


def mono(e, c):
    result = F(1)
    for power, x in zip(e, c):
        result *= x**power
    return result


def derivative(e, c, axes):
    e = list(e)
    factor = 1
    for axis in axes:
        factor *= e[axis]
        if not factor:
            return F(0)
        e[axis] -= 1
    return factor*mono(e, c)


def spherical_conjugation(e, q, c):
    """Differentiate f=(1-q*r^2)^(1/4)*c^e in the original sphere operator.

    Return J^(1/2)*(-Delta_S3/2)*f and the claimed flattened expression.
    The factor q implements the coordinate dilation of the sphere.
    """
    y = sum(x*x for x in c)
    den = 1-q*y
    p = mono(e, c)
    grad = [derivative(e, c, (i,)) for i in range(3)]
    hess = [[derivative(e, c, (i, j)) for j in range(3)] for i in range(3)]
    log_first = [-q*x/(2*den) for x in c]
    factor_second = [[-q*F(i == j)/(2*den)-3*q*q*c[i]*c[j]/(4*den**2)
                      for j in range(3)] for i in range(3)]
    B = [[F(i == j)-q*c[i]*c[j] for j in range(3)] for i in range(3)]
    original = -sum(B[i][j]*(hess[i][j]+log_first[i]*grad[j]
                            +log_first[j]*grad[i]+factor_second[i][j]*p)
                    for i in range(3) for j in range(3))/2
    original += F(3, 2)*q*sum(c[i]*(grad[i]+log_first[i]*p) for i in range(3))
    Q = m.geometry(q, c)[2]
    flat = -sum(B[i][j]*hess[i][j] for i in range(3) for j in range(3))/2
    flat += 2*q*sum(c[i]*grad[i] for i in range(3))+Q*p
    return original, flat


class CM2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_packet_and_source_pins(self):
        self.assertTrue(self.result['all_pass'])
        self.assertTrue(all(self.result['checks'].values()))
        stored = json.loads((Path(__file__).parent/'CM2_RESULT.json').read_text())
        self.assertEqual(stored, dict(self.result, source_sha256=m.source_pins()))

    def test_frozen_pair_against_full_gram_action(self):
        data = m.pair_data()
        basis, s, h = m.cr.ladder(2, data['width'], 16)
        vector = [data['p'].get(e, F(0)) for e in basis]
        image = [data['hp'].get(e, F(0)) for e in basis]
        self.assertEqual(m.gc.quad(s, vector), data['norm'])
        self.assertEqual(m.gc.quad(h, vector)/data['norm'], data['eta'])
        self.assertEqual(m.gc.quad(s, image)/data['norm'], data['second'])
        self.assertGreater(data['line3'], m.RHO)
        self.assertEqual(len(data['p']), 64)

    def test_original_spherical_generator_after_half_density(self):
        for q, c in ((F(1, 7), (F(1, 2), F(1, 3), F(2, 5))),
                     (F(2, 3), (F(0), F(2, 5), F(1, 3))),
                     (F(0), (F(1), F(2), F(3)))):
            for e in product(range(5), repeat=3):
                if sum(e) <= 4:
                    original, flattened = spherical_conjugation(e, q, c)
                    self.assertEqual(original, flattened)
        q, c = F(1, 7), (F(1, 2), F(1, 3), F(2, 5))
        actual, _ = spherical_conjugation((0, 0, 0), q, c)
        Q = m.geometry(q, c)[2]
        self.assertEqual(actual, Q)
        self.assertGreater(Q, 0)
        self.assertNotEqual(actual, 0, 'dropping the scalar metric correction must fail')
        self.assertNotEqual(actual, -Q, 'reversing its sign must fail')

    def test_exact_ims_energy_and_omission_control(self):
        B, _, _ = m.geometry(F(1, 7), (F(1, 3), F(1, 2), F(2, 5)))
        norm2 = lambda v: sum(v[i]*B[i][j]*v[j] for i in range(3) for j in range(3))
        for z in (F(1, 3), F(1), F(7, 2)):
            c, s = (1-z*z)/(1+z*z), 2*z/(1+z*z)
            value, grad, angle_grad = F(3, 2), (F(1), F(-2), F(3)), (F(1, 5), F(2, 7), F(-1, 4))
            first = [c*grad[i]-value*s*angle_grad[i] for i in range(3)]
            second = [s*grad[i]+value*c*angle_grad[i] for i in range(3)]
            self.assertEqual(norm2(first)+norm2(second), norm2(grad)+value**2*norm2(angle_grad))
            self.assertGreater(norm2(first)+norm2(second), norm2(grad))
        # Flat metric would overstate the kinetic lower along the radial direction.
        radial = (F(1, 3), F(1, 2), F(2, 5))
        self.assertLess(norm2(radial), sum(x*x for x in radial))

    def test_centre_intertwiner_and_global_count_control(self):
        signs = list(product((-1, 1), repeat=3))
        for sigma in signs:
            for sheet in signs:
                for axis in range(3):
                    moved = list(sheet)
                    moved[axis] *= -1
                    self.assertEqual(m.character(sigma, tuple(moved)),
                                     sigma[axis]*m.character(sigma, sheet))
        # Each sector has one level below 3; the direct sum has eight, not one.
        toy_sector = (F(2), F(4))
        self.assertEqual(sum(e < 3 for e in toy_sector), 1)
        combined = sorted(toy_sector*8)
        self.assertEqual(sum(e < 3 for e in combined), 8)
        self.assertEqual(combined[1]-combined[0], 0)
        self.assertEqual(toy_sector[1]-toy_sector[0], 2)

    def test_uniform_window_monotonicity_and_complete_cost(self):
        for row in self.result['weak_windows']:
            start = row['theta_start']
            radius, width = F(row['localization_radius']), F(row['localization_width'])
            initial = m.cb.compact_upper(start, F(6), F(1, 4), 28, F(1, 500000))['upper']
            for theta in (start, 8*start):
                local = m.localization(theta, radius, width)
                upper = m.cb.compact_upper(theta, F(6), F(1, 4), 28, F(1, 500000))['upper']
                self.assertLessEqual(upper, initial)
                coefficient = local['excitation_coefficient']-2*upper
                self.assertGreaterEqual(coefficient, F(row['advertised_coefficient']))
                self.assertLess(local['excitation_coefficient'], 2*m.RHO)
                self.assertGreater(local['subtracted_constant'], 0)
                w = m.cb.root_floor(F(theta), 3)
                inner = (1-local['reach']**2/w)*m.RHO+F(9, 4)/w
                outer = radius-3/w
                self.assertGreaterEqual(outer, inner)

    def test_compact_normalization_and_claim_boundary(self):
        for epsilon in (F(1, 3), F(2, 7), F(1, 100)):
            theta = epsilon**(-6)
            w = epsilon**(-2)
            self.assertEqual(theta*epsilon**4, w)
            # H = -Delta + 2 theta X corresponds to 2w*(-Delta/2 + X).
            self.assertEqual(2*theta*epsilon**4, 2*w)
        boundary = self.result['claim_boundary']
        self.assertFalse(boundary['volume_uniformity'])
        self.assertFalse(boundary['continuum_mass_gap'])
        self.assertFalse(boundary['norm_resolvent_convergence'])
        self.assertIn('-> 0', boundary['gap_across_all_centre_sectors'])

    def test_invalid_chart_and_uncertified_windows_refused(self):
        for args in ((F(-1), (F(0),)*3), (F(1), (F(1), F(0), F(0)))):
            with self.assertRaises(ValueError):
                m.geometry(*args)
        for args in ((0, F(6), F(3)), (1, F(6), F(3)), (10**10, F(5), F(4)),
                     (10**10, F(6), F(-1))):
            with self.assertRaises(ValueError):
                m.localization(*args)
        with self.assertRaises(ValueError):
            m.weak_window(2*10**9, F(31, 10), F(1), F(878))
        with self.assertRaises(ValueError):
            m.character((1, 0, 1), (1, 1, 1))


if __name__ == '__main__':
    unittest.main()
