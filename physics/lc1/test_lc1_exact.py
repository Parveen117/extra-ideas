"""LC1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import lc1_lost_is_counted as y


class LostIsCountedTests(unittest.TestCase):
    def test_all(self):
        res = y.run()
        self.assertEqual(res['elements'][0], ('2', '17/8', '9/16'))
        self.assertEqual(res['ladder'][0], (1, '1/2', ['3/2', '5/2', '7/2'], ['1', '2', '3']))
        self.assertEqual(res['ladder'][1][3], ['1/2', '1', '3/2'])
        self.assertEqual(res['gravity'][3], ('1/3', '1/2', '1'))

    def test_the_root_is_needed(self):
        """the lost part of H itself is not the return; that of its root is."""
        M = y.element(F(2))
        H = y.ln.mul(M, M)
        self.assertNotEqual(y.ln.lost(H)[0], (H[0][0] - 1)/2)
        self.assertEqual(y.ln.lost(M)[0], (H[0][0] - 1)/2)

    def test_wrong_lost_part_is_rejected(self):
        with patch.object(y.ln, 'lost', lambda H: [H[0][1], H[0][1]]):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
