"""IN1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import in1_the_invariant as y


class InvariantTests(unittest.TestCase):
    def test_single_readings_are_null_and_complete_with_three_cuts(self):
        rows = y.single_control()
        self.assertEqual(rows[0]['readings'], ['8', '6', '0'])
        self.assertTrue(rows[0]['two_cuts_complete'])
        self.assertFalse(rows[2]['two_cuts_complete'])

    def test_dropping_the_third_cut_breaks_the_tensor(self):
        with patch.object(y, 'C3', y.C2):
            with self.assertRaises(ValueError):
                y.readings(y.tensor(y.STATES[2]))

    def test_frame_independence(self):
        rows = y.frame_control()
        self.assertEqual([r['invariant'] for r in rows], ['0', '0', '0', '808/9'])

    def test_a_frame_change_of_determinant_two_is_rejected(self):
        with patch.object(y, 'FRAMES', [[[y.c(2), y.Z], [y.Z, y.c(1)]]]):
            with self.assertRaises(ValueError):
                y.frame_control()

    def test_records(self):
        rows = y.record_control()
        self.assertEqual(rows[0]['invariant_at_half'], '10')
        self.assertTrue(rows[-1]['parallel'])

    def test_wrong_share_law_is_rejected(self):
        real = y.form
        with patch.object(y, 'form', lambda n, r: real(n, r)+1):
            with self.assertRaises(ValueError):
                y.record_control()


if __name__ == '__main__':
    unittest.main()
