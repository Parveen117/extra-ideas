"""PR2: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import pr2_two_sheets_pair_observer as y


class TwoSheetTests(unittest.TestCase):
    def test_form_and_sheets(self):
        rows = y.form_control()
        self.assertEqual([r['kind'] for r in rows], ['timelike', 'timelike', 'timelike', 'spacelike'])

    def test_a_non_boost_is_rejected(self):
        with self.assertRaises(ValueError):
            y.boost((F(1), F(0)), F(5, 4), F(1))

    def test_energy_momentum(self):
        row = y.energy_momentum_control(F(3, 2))[0]
        self.assertEqual(row['energy'], '51/16')
        self.assertEqual(F(row['energy'])**2-F(row['momentum_squared']), F(9, 4))

    def test_balanced_pair_is_frame_independent(self):
        row = y.pair_control(F(3, 5), F(4, 5), F(1, 2))
        self.assertTrue(all(v == ['0', '0', '0'] for v in row['odd_part_by_frame']))
        self.assertEqual(row['memory'], '16/25')

    def test_single_reading_is_frame_dependent(self):
        row = y.pair_control(F(3, 5), F(4, 5), F(1))
        self.assertEqual(row['memory'], '0')
        self.assertGreater(len({tuple(v) for v in row['odd_part_by_frame']}), 1)

    def test_one_sided_frame_change_breaks_the_scalars(self):
        with patch.object(y, 'act', lambda g, X: y.mm(g[0], X)):
            with self.assertRaises(ValueError):
                y.pair_control(F(3, 5), F(4, 5), F(3, 4))

    def test_three_sectors(self):
        self.assertEqual([y.sector_control(a, F(1))['square'] for a in (F(1, 2), F(1), F(2))], ['-3/4', '0', '3'])


if __name__ == '__main__':
    unittest.main()
