"""PH2 exact part: stdlib only."""
from fractions import Fraction as F
import unittest

import ph2_strain_cycle_rotation as y


class ExactTests(unittest.TestCase):
    def test_two_pure_strains_leave_a_rotation(self):
        row = y.exact_witness()
        self.assertEqual(row['tan_theta'], '288/541')
        self.assertEqual(row['aligned_tan_theta'], '0')

    def test_polar_formula_rejects_a_wrong_angle(self):
        g = [[F(2), F(1)], [F(-1), F(3)]]
        num, den = y.polar_tangent(g)
        self.assertEqual(den*g[0][1]+num*g[1][1], -num*g[0][0]+den*g[1][0])
        num += 1
        self.assertNotEqual(den*g[0][1]+num*g[1][1], -num*g[0][0]+den*g[1][0])


if __name__ == '__main__':
    unittest.main()
