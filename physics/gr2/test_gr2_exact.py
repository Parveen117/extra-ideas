"""GR2: stdlib only (one float route for the once-around turn)."""
from fractions import Fraction as F
import math
import unittest
from unittest.mock import patch

import gr2_thesis_core as y


def sqrt_spd(a):
    """Square root of a symmetric positive 2x2 matrix, closed form."""
    s = math.sqrt(a[0][0]*a[1][1]-a[0][1]*a[1][0])
    t = math.sqrt(a[0][0]+a[1][1]+2*s)
    return [[(a[0][0]+s)/t, a[0][1]/t], [a[1][0]/t, (a[1][1]+s)/t]]


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def inv(a):
    d = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[a[1][1]/d, -a[0][1]/d], [-a[1][0]/d, a[0][0]/d]]


class ThesisCoreTests(unittest.TestCase):
    def test_recoverable_memory_is_removed_by_one_frame_change(self):
        rows = y.equivalence_control()
        self.assertTrue(all(r['memory_after_frame_change'] == '0' for r in rows))

    def test_the_wrong_frame_change_does_not_remove_it(self):
        real = y.ms1.boostK
        with patch.object(y.ms1, 'boostK', lambda b: real(b*b if b < 1 else b)):
            with self.assertRaises(ValueError):
                y.equivalence_control()

    def test_the_source_survives_every_frame_change(self):
        rows = y.source_control()
        self.assertEqual(rows[0]['unrecoverable'], '1225/144')
        self.assertGreater(len(rows[0]['recoverable_by_frame']), 1)

    def test_fluxes_add(self):
        rows = y.linearity_control()
        self.assertEqual(rows[0]['flux_per_unit'], '-1/2')

    def test_once_around_against_a_chain_of_responses(self):
        # responses H(phi) = cosh(eta) + sinh(eta)(cos(phi) K + sin(phi) RK), cosh(eta) = 1/N = 5/4
        ch, sh = 1.25, 0.75
        steps = 6000
        Hs = [[[ch+sh*math.cos(p), sh*math.sin(p)], [sh*math.sin(p), ch-sh*math.cos(p)]]
              for p in (2*math.pi*k/steps for k in range(steps))]
        f = sqrt_spd(Hs[0])
        f0 = [row[:] for row in f]
        for h in Hs[1:]+[Hs[0]]:
            fi = inv(f)
            a = mul(mul([[fi[0][0], fi[1][0]], [fi[0][1], fi[1][1]]], h), fi)
            f = mul(sqrt_spd(a), f)
        rot = mul(f, inv(f0))
        angle = abs(math.atan2(rot[1][0]-rot[0][1], rot[0][0]+rot[1][1]))
        self.assertAlmostEqual(angle, math.pi*(ch-1), places=5)          # probe-plane turn pi (1/N - 1)
        self.assertEqual(y.around_control()[0]['turn_over_2pi'], '1/4')  # a direction turns twice that

    def test_difference_from_the_cone_value(self):
        row = y.around_control()[0]
        self.assertEqual(row['difference'], '1/20')


if __name__ == '__main__':
    unittest.main()
