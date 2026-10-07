"""PR4: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import pr4_source_is_memory as y

G = F(3, 2)


class SourceTests(unittest.TestCase):
    def test_currents_and_residue(self):
        self.assertEqual(y.current_control(G)['residue'], 'd_t j - d_x n = -2 g sigma')

    def test_residue_with_a_split_coupling_is_rejected(self):
        real = y.rhs

        def split(psi, g):                    # g K instead of g R
            p1, p2 = psi
            return (y.fadd(y.d_x(p1), y.scale(g, p2)), y.fadd(y.scale(F(-1), y.d_x(p2)), y.scale(g, p1)))
        with patch.object(y, 'rhs', split):
            with self.assertRaises(ValueError):
                y.current_control(G)
        self.assertIs(y.rhs, real)

    def test_frame_residue(self):
        self.assertTrue(y.frame_control(G))

    def test_memory_is_proper_density(self):
        rows = y.invariance_control()
        self.assertEqual(rows[0]['memory'], rows[0]['sigma_over_n_squared'])
        self.assertEqual(rows[-1]['memory'], '0')

    def test_uniqueness_and_sheet_parity(self):
        self.assertEqual(y.uniqueness_control()['invariant_entries'], ['12', '21'])
        self.assertEqual(y.sheet_control(G)['tau'], 'even')


if __name__ == '__main__':
    unittest.main()
