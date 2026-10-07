"""OB1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import ob1_one_block as y


class OneBlockTests(unittest.TestCase):
    def test_algebra(self):
        self.assertEqual(y.algebra_control()['volume'], 'iota = C1 C2 C3')

    def test_field_equations_as_grades(self):
        row = y.field_control()
        self.assertTrue(row['source_real'] and row['source_conserved'])

    def test_a_third_cut_without_iota_breaks_the_algebra(self):
        broken = [y.C[0], y.C[1], y.bconst([[(0, 0), (-1, 0)], [(1, 0), (0, 0)]])]     # R instead of iota R
        with patch.object(y, 'C', broken):
            with self.assertRaises(ValueError):
                y.algebra_control()

    def test_wave_and_its_wrong_sheet(self):
        self.assertTrue(y.wave_control()['null'])

    def test_matter_reading_is_conserved_and_null(self):
        self.assertTrue(y.matter_control()['null'])

    def test_force(self):
        rows = y.force_control()
        self.assertEqual(rows[0]['power'], '-3')
        self.assertEqual(rows[0]['force'], ['-3', '7', '1'])
        self.assertEqual(rows[3]['force'], ['12', '0', '0'])

    def test_a_non_reading_potential_is_rejected(self):
        real = y.potentials

        def complex_potential():
            phi, A = real()
            return y.padd(phi, y.pscale((0, 1), y.var(1))), A
        with patch.object(y, 'potentials', complex_potential):
            with self.assertRaises(ValueError):
                y.field_control()


if __name__ == '__main__':
    unittest.main()
