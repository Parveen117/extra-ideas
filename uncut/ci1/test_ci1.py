"""CI1: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import ci1_shared_information as y


class SharedInformationTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['arrow']['witness_exp_2I0'], ['289/244', '5/4', '16/11'])
        self.assertEqual(res['sectors']['exp_2I0'], {'cut': '25/17', 'cone': '1', 'turn': '25/33'})

    def test_I0_does_not_depend_on_the_frame_and_I_cut_does(self):
        L = sp.Matrix([[3, 1], [1, 2]])
        turn = sp.Rational(3, 5)*y.I2 + sp.Rational(4, 5)*y.R
        boost = sp.Rational(5, 4)*y.I2 + sp.Rational(3, 4)*y.K
        for W in (turn, boost, turn*boost):
            M = W*L*W.inv()
            self.assertEqual(y.e2I_0(M), y.e2I_0(L))
        self.assertNotEqual(y.e2I_cut(turn*L*turn.inv()), y.e2I_cut(L))

    def test_zero_cases(self):
        self.assertEqual(y.e2I_cut(sp.Matrix([[3, 0], [0, 2]])), 1)          # no coupling seen by this cut
        self.assertEqual(y.e2I_0(sp.Matrix([[2, 0], [0, 2]])), 1)            # isotropic: nothing shared in any cut

    def test_opposite_sense_lowers_it(self):
        a, b = sp.Rational(5, 4), sp.Rational(13, 12)                        # cosh values with sinh 3/4 and 5/12
        sa, sb = sp.Rational(3, 4), sp.Rational(5, 12)
        self.assertLess(a*b - sa*sb, a*b)                                    # cosh(a - b)
        self.assertGreater(a*b + sa*sb, a*b)                                 # cosh(a + b)

    def test_cycle_of_tt1(self):
        kap = sp.Rational(5, 4)*y.K + sp.Rational(3, 4)*y.R                  # w = 3/4
        one = (y.R*kap)**2
        two = one*one
        self.assertEqual(y.e2I_0(one), sp.Rational(17, 8)**2)
        self.assertGreater(y.e2I_0(two), y.e2I_0(one)**2)

    def test_reading_the_other_channel_is_rejected(self):
        def wrong(M):
            a, b, c, d = y.parts(M)
            return (a*a - c*c)/M.det()
        with patch.object(y, 'e2I_cut', wrong):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
