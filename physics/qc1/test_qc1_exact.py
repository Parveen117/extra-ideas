"""QC1 exact part: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import qc1_reading_content_closure as y


class ExactTests(unittest.TestCase):
    def test_pair_bracket_and_graph(self):
        self.assertEqual(len(y.pair_controls()), 3)

    def test_wrong_share_leaves_the_graph(self):
        real = y.scal
        with patch.object(y, 'scal', lambda c, a: real(F(1, 3), a) if c == F(1, 2) else real(c, a)):
            with self.assertRaises(ValueError):
                y.pair_controls()

    def test_content_spectrum(self):
        self.assertEqual(y.content_control(F(3, 5), F(4, 5), 2), [-2, 0, 2])
        self.assertEqual(y.content_control(F(5, 13), F(12, 13), 5), [-5, -3, -1, 1, 3, 5])

    def test_a_non_rotation_breaks_the_content_law(self):
        with self.assertRaises(ValueError):
            y.content_control(F(3, 5), F(3, 5), 2)

    def test_weight_ladder(self):
        lad = y.ladder(2, 3)
        self.assertEqual([c['height'] for c in lad['circles']], ['2', '3', '4'])
        self.assertEqual(lad['circles'][0]['klein_radius_squared'], '3/4')
        self.assertEqual(y.ladder(1, 2)['index_k'], '1/2')


if __name__ == '__main__':
    unittest.main()
