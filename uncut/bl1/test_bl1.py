"""BL1: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import bl1_balance_of_shared_information as y


class BalanceTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['cycles']['exp_I0_after_k_cycles_w_3_4'][:2], ['17/8', '257/32'])

    def test_self_dagger_cut_adds_nothing(self):
        kap = sp.Rational(3, 5)*y.K + sp.Rational(4, 5)*y.S
        self.assertEqual((y.R*kap)**2, y.I2)

    def test_channels_in_a_turned_frame(self):
        M = sp.Matrix([[3, 1], [2, 2]])
        W = sp.Rational(3, 5)*y.I2 + sp.Rational(4, 5)*y.R
        Mt = W*M*W.inv()
        s1, l1 = y.channels(M)
        s2, l2 = y.channels(Mt)
        self.assertEqual(s1 + l1, s2 + l2)
        self.assertNotEqual(s1, s2)

    def test_wrong_channel_split_is_rejected(self):
        def wrong(M, across=False):
            a, b, c, d = y.parts(M)
            return (a*a, d*d - c*c)
        with patch.object(y, 'channels', wrong):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
