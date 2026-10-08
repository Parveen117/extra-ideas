"""CO1: sympy."""
import unittest
from unittest.mock import patch

import sympy as sp

import co1_expansion_as_fall as y


class ExpansionAsFallTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['law_value'], '-6 h^2 in both frames')
        self.assertEqual(res['pressure'], '3*kappa*n_0**2/(16*A**6)')

    def test_fall_frame_has_flat_slices_and_one_time(self):
        E = y.fall_frame()
        self.assertEqual(E[1:, 1:], sp.eye(3))
        self.assertEqual(E[0, 0], 1)

    def test_it_is_one_frame_in_two_coordinate_systems(self):
        # every invariant agrees separately, not only the selected combination
        Ic = y.tp.invariants(y.tp.structure_in(y.X, y.comoving_frame()))
        If = y.tp.invariants(y.tp.structure_in(y.X, y.fall_frame()))
        for u, v in zip(Ic, If):
            self.assertEqual(sp.simplify(u - v), 0)
        self.assertEqual([sp.nsimplify(sp.simplify(u/y.h**2)) for u in If], [6, 3, 9])

    def test_a_wrong_law_is_rejected(self):
        def other(E):
            I1, I2, I3 = y.tp.invariants(y.tp.structure_in(y.X, E))
            return sp.simplify(I1 + I2 + I3)
        with patch.object(y, 'law', other):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
