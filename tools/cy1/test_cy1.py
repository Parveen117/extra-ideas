"""CY1: stdlib only."""
from fractions import Fraction as F
import unittest

import cy1_three_reading_cyclic_form as m


class CyclicFormTests(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)

    def test_by_hand_in_two_dimensions(self):
        a, b, c = [F(1), F(0)], [F(0), F(1)], [F(0), F(0)]
        t = m.cyclic(a, b, c)
        self.assertEqual(t, [[0, 1], [0, 0]])
        self.assertEqual(m.madd(t, m.mscale(m.transpose(t), -1)), m.wedge(a, b))

    def test_a_non_common_factor_changes_the_alternating_part(self):
        p = (F(2), F(3))
        base = [m.dlog_monomial(e, p) for e in ((1, 0), (0, 1), (0, 0))]
        moved = [m.dlog_monomial(e, p) for e in ((1, 0), (0, 1), (2, 5))]      # only K rescaled
        alt = lambda t: m.madd(t, m.mscale(m.transpose(t), -1))
        self.assertNotEqual(alt(m.cyclic(*base)), alt(m.cyclic(*moved)))


if __name__ == '__main__':
    unittest.main()
