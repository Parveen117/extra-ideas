"""FC1: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import fc1_what_the_first_cut_leaves_out as y


class FirstCutTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['witness']['moving'], 'E = 5, p = 3, W = 3 ; E^2 - p^2 = 16 ; det = 25')

    def test_opposite_turn_parts_make_a_flat_pair(self):
        r1 = 4*y.I2 + 3*y.R
        r2 = 4*y.I2 - 3*y.R
        tot = r1 + r2
        self.assertEqual(y.parts(tot), (8, 0, 0, 0))
        self.assertEqual(tot.det(), 64)
        self.assertEqual(r1.det() + r2.det(), 50)

    def test_conjugation_is_a_different_action(self):
        rho = 4*y.I2 + 1*y.K + 3*y.R
        Bst = sp.Matrix([[2, 0], [0, sp.Rational(1, 2)]])
        conj = Bst*rho*Bst.inv()
        cong = Bst*rho*y.dagger(Bst)
        self.assertEqual(y.parts(conj)[0], 4)                 # conjugation keeps the scalar part ...
        self.assertNotEqual(y.parts(conj)[3], 3)              # ... and changes the turn-part
        self.assertEqual(y.parts(cong)[3], 3)                 # the change of frame keeps the turn-part
        self.assertNotEqual(y.parts(cong)[0], 4)

    def test_with_the_ordinary_transpose_removed_the_stage_fails(self):
        with patch.object(y, 'dagger', lambda M: M):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
