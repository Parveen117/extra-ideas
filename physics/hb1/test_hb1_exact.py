"""HB1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import hb1_unit_of_the_frame as y


class UnitOfTheFrameTests(unittest.TestCase):
    def test_commutators(self):
        self.assertEqual(y.commutator_control()['cyclic_sizes'], [3, 4, 5])

    def test_dropping_the_wrap_breaks_the_cyclic_identity(self):
        real = y.even
        # pretend R = C on cyclic marks: must fail
        def no_wrap(a, cyclic=False):
            out = real(a, cyclic)
            if cyclic:
                Q = len(a)
                out = [y.cadd(x, y.cmul(y.c(F(Q, 2)), w)) for x, w in zip(out, [a[Q-1]]+[y.Z]*(Q-2)+[a[0]])]
            return out
        with patch.object(y, 'even', no_wrap):
            with self.assertRaises(ValueError):
                y.commutator_control()

    def test_inequality_and_sharp_instance(self):
        row = y.inequality_control()
        self.assertEqual(row['sharp_instance']['R'], '-1')
        self.assertEqual(row['states'][-1]['ratio'], '15/14')

    def test_unit_is_tick_even_part_times_overlap(self):
        rows = y.phase_control()
        first = rows[0]
        self.assertEqual(first['unit'], '9/20')
        quarter = [r for r in rows if r['tick'] == ['0', '1']]
        self.assertTrue(all(r['unit'] == '0' for r in quarter))
        flat = [r for r in rows if r['tick'] == ['1', '0'] and r['profile'] == 'box 12']
        self.assertEqual(flat[0]['unit'], '11/12')

    def test_a_tick_that_is_not_a_turn_breaks_the_reading(self):
        real = y.cmul
        with patch.object(y, 'conj', lambda a: a):
            with self.assertRaises(ValueError):
                y.phase_control()
        self.assertIs(y.cmul, real)


if __name__ == '__main__':
    unittest.main()
