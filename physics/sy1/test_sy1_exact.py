"""SY1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import sy1_one_structure as y


class OneStructureTests(unittest.TestCase):
    def test_emk1_identities_on_the_grid(self):
        self.assertEqual(y.relations_control(), 1296)

    def test_unit_blocks_in_three_sectors(self):
        rows = y.unit_control()
        self.assertEqual(sorted({r['sector'] for r in rows}), ['circular', 'dual', 'split'])
        self.assertTrue(all((F(r['defect']) > 0) == (r['sector'] == 'circular') or r['sector'] == 'dual' for r in rows))

    def test_off_quadric_block_is_rejected(self):
        with patch.dict(y.UNITS, {'circular': [y.block(F(3, 5), F(3, 5), F(1), F(1, 5))]}):
            with self.assertRaises(ValueError):
                y.unit_control()

    def test_sector_mislabel_is_rejected(self):
        with patch.dict(y.UNITS, {'circular': y.UNITS['split']}):
            with self.assertRaises(ValueError):
                y.unit_control()

    def test_composition_law(self):
        self.assertEqual(y.composition_control(), 121)

    def test_stage_laws_as_cases(self):
        cases = y.cases_control()
        self.assertEqual(cases['PH2 rotation tangent'], '26/127')
        self.assertEqual(cases['MS1 state: speed, memory'], ['35/37', '144/1369'])
        self.assertEqual(cases['CL1 mass^2 = clock curvature = det(I - coin)'], '4/5')

    def test_a_wrong_quarter_turn_breaks_everything(self):
        with patch.object(y, 'R', y.K):
            with self.assertRaises(ValueError):
                y.relations_control()


if __name__ == '__main__':
    unittest.main()
