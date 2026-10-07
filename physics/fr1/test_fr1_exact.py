"""FR1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import fr1_the_frame_cuts as y


class FrameTests(unittest.TestCase):
    def test_one_invariant_many_readings(self):
        rows = y.cut_control()
        self.assertEqual(rows[0]['readings'][0], ['1', '0'])
        self.assertEqual(rows[0]['readings'][3], ['25/16', '-9/16'])
        self.assertFalse(rows[1]['cut_dependent'])
        self.assertTrue(rows[3]['cut_dependent'])

    def test_a_non_cut_is_rejected(self):
        with self.assertRaises(ValueError):
            y.channels(y.block(F(1), F(0), F(0), F(0)), F(3, 5), F(3, 5))

    def test_only_the_scalar_survives_a_frame_change(self):
        row = y.frame_control()
        self.assertNotEqual(row['turn_component_before'], row['turn_component_after'])

    def test_sheets(self):
        forms = [r['form'] for r in y.sheets_control()]
        self.assertEqual(forms, ['1', '1', '1', '1', '-1', '-1', '0'])

    def test_counting_of_cuts(self):
        row = y.counting_control()
        self.assertEqual((row['real_block_cuts'], row['cut_complex_cuts'], row['linear_fourth'], row['conjugate_family']), (2, 3, 0, 2))

    def test_without_the_turn_there_is_no_third_cut(self):
        # replacing iota R by K (not anticommuting with K) must be caught
        real = y.kron
        calls = {'n': 0}

        def fake(a, b):
            calls['n'] += 1
            if a == y.R:
                return real(y.K, b)
            return real(a, b)
        with patch.object(y, 'kron', fake):
            with self.assertRaises(ValueError):
                y.counting_control()


if __name__ == '__main__':
    unittest.main()
