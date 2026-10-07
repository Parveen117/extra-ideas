"""DM1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import dm1_one_more_dimension as y


class OneMoreDimensionTests(unittest.TestCase):
    def test_hypotenuse_and_legs(self):
        rows = y.pythagoras_control()
        self.assertEqual(rows[1]['seen_speed_squared'], '9/25')
        self.assertEqual(rows[1]['unseen_squared'], '16/25')
        self.assertEqual(rows[0]['unseen_squared'], '0')

    def test_third_reading_is_the_two_cut_invariant(self):
        rows = y.third_reading_control()
        self.assertEqual(rows[0]['third_reading'], '4')
        self.assertNotEqual(rows[0]['after_a_cut_complex_frame_change'], '4')

    def test_a_cut_complex_frame_is_not_a_two_cut_frame(self):
        with patch.object(y, 'REAL_FRAMES', [y.COMPLEX_FRAME]):
            with self.assertRaises(ValueError):
                y.third_reading_control()

    def test_mass_is_the_third_momentum(self):
        rows = y.mass_control()
        self.assertEqual(rows[0]['flip_rate'], '-3')
        self.assertEqual({r['sheet'] for r in rows}, {'upper', 'lower'})

    def test_without_iota_the_third_term_is_not_a_turn(self):
        with patch.object(y, 'C3', y.RMAT):
            with self.assertRaises(ValueError):
                y.mass_control()

    def test_record_is_single_one_level_up(self):
        rows = y.doubled_control()
        self.assertEqual(rows[0]['unrecoverable'], '148')
        self.assertTrue(rows[2]['single'])


if __name__ == '__main__':
    unittest.main()
