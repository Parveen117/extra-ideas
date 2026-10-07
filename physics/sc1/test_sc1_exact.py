"""SC1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import sc1_local_scale as y


class LocalScaleTests(unittest.TestCase):
    def test_scale_potential_changes_the_count_and_phase_keeps_it(self):
        row = y.balance_control()
        self.assertEqual(row['phase_potential'], 'count kept')

    def test_a_phase_potential_that_is_not_self_dagger_is_caught(self):
        real = y.potential_block
        calls = {'n': 0}

        def twisted(w0, w):
            calls['n'] += 1
            out = real(w0, w)
            if calls['n'] == 2:                               # the phase potential
                out = y.ob.badd(out, y.ob.bscale(y.pc(0, 1), y.ONE))
            return out
        with patch.object(y, 'potential_block', twisted):
            with self.assertRaises(ValueError):
                y.balance_control()

    def test_exact_potential_has_no_field(self):
        row = y.exact_control()
        self.assertEqual(row['exact_potential_field'], '0')

    def test_histories(self):
        row = y.history_control()
        self.assertEqual(row['count_ratio'], '9')
        self.assertEqual(row['invariant_ratio'], '81')

    def test_scale_is_not_a_unit_block(self):
        M = y.in1.scale(F(3), y.in1.ONE)
        self.assertNotEqual(y.in1.det(M), y.in1.O)


if __name__ == '__main__':
    unittest.main()
