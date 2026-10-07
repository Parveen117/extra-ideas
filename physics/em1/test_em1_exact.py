"""EM1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import em1_wave_cut as y


class WaveCutTests(unittest.TestCase):
    def test_identity_and_single_wave(self):
        rows = y.identity_control()
        self.assertEqual(rows[0]['unrecoverable_squared'], '41/4')
        self.assertTrue(rows[2]['single'])

    def test_frame_change_keeps_the_invariant_and_moves_the_readings(self):
        rows = y.frame_control()
        self.assertEqual(rows[0]['E_squared'], ['5', '13/8'])
        self.assertEqual(rows[0]['E2_minus_B2'], '-5')
        self.assertEqual(rows[3]['B_squared'], ['0', '9/16'])

    def test_a_real_turn_is_not_the_frame_change(self):
        def real_turn(Fv, ch, sh):
            fy = y.cadd(y.cmul(y.c(ch), Fv[1]), y.cmul(y.c(sh), Fv[2]))
            fz = y.cadd(y.cmul(y.c(ch), Fv[2]), y.cmul(y.c(-sh), Fv[1]))
            return [Fv[0], fy, fz]
        with patch.object(y, 'boost', real_turn):
            with self.assertRaises(ValueError):
                y.frame_control()

    def test_two_waves(self):
        row = y.two_wave_control()
        self.assertEqual(row['counter_running']['unrecoverable_squared'], '144')
        self.assertEqual(row['parallel_unrecoverable'], '0')

    def test_train_energy_over_frequency(self):
        rows = y.doppler_control()
        self.assertEqual(rows[0]['doppler'], '2')
        self.assertTrue(all(r['train_energy_over_frequency'] == 'unchanged' for r in rows))


if __name__ == '__main__':
    unittest.main()
