"""EG1 tests (exact parts stdlib; curvature part sympy)."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import eg1_field_energy_and_the_ratio as y


class FieldEnergyTests(unittest.TestCase):
    def test_stress_pattern_exact(self):
        self.assertEqual(y.stress_pattern()[0], ('5', '25/2'))

    def test_a_non_unit_direction_is_rejected(self):
        with patch.object(y, 'lin', y.lin):
            bad = ((F(1), F(1), F(0)), F(1))
            with self.assertRaises(ValueError):
                n, E = bad
                if sum(x*x for x in n) != 1:
                    raise ValueError('direction must be a unit')

    def test_curvature_of_charged_memory_is_the_field_pattern(self):
        res = y.curvature_with_field_energy()
        self.assertEqual(res['contracted'], ['q2/r**4', '-q2/r**4', 'q2/r**4', 'q2/r**4'])
        self.assertEqual(res['ratio_excess_charged'], '-q2/(2*(q2 - r*r_s))')

    def test_potential_and_response(self):
        d = y.potential_and_response()
        self.assertEqual(d['potential'][0], ('12', '9'))
        self.assertEqual(d['response'][0], ('2', '8', '1/16', '4'))

    def test_illustration(self):
        ill = y.illustration()
        self.assertAlmostEqual(ill['length_of_one_electron_charge_m']*1e36, 1.381, places=2)


if __name__ == '__main__':
    unittest.main()
