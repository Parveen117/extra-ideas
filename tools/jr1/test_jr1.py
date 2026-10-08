"""JR1: stdlib only."""
from fractions import Fraction as F
import unittest

import jr1_jet_reading_of_the_change_operator as m


class JetReadingTests(unittest.TestCase):
    def test_all(self):
        for k, v in m.run().items():
            self.assertTrue(v, k)

    def test_iota_squares_to_minus_one(self):
        self.assertEqual(m.IOTA*m.IOTA, m.G(-1))

    def test_double_pole_by_hand(self):
        # f = a t^-2 + b t^-1 :  reading = b + a c ,  c = iota (lam - 1)
        a, b, lam = m.G(2, 1), m.G(-1, 3), m.G(F(1, 2), 2)
        f = {-2: a, -1: b}
        self.assertEqual(m.reading(f, m.c_of(lam)), b + a*m.IOTA*(lam - 1))
        self.assertEqual(m.reading(m.change(f), m.c_of(lam)), lam*m.reading(f, m.c_of(lam)))

    def test_regular_part_is_not_read(self):
        f = {-2: m.G(1), -1: m.G(4)}
        h = dict(f)
        h.update({0: m.G(9), 3: m.G(0, 5)})
        self.assertEqual(m.reading(f, m.G(1, 1)), m.reading(h, m.G(1, 1)))

    def test_singular_recovery_is_refused(self):
        self.assertIsNone(m.solve([[m.G(1), m.G(2)], [m.G(2), m.G(4)]], [m.G(1), m.G(0)]))


if __name__ == '__main__':
    unittest.main()
