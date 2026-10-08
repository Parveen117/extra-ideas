"""RB1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import rb1_return_boost_boundary_memory as y


class ReturnBoostBoundaryTests(unittest.TestCase):
    def test_run_and_pinned_values(self):
        res = y.run()
        self.assertEqual(res['uniform'], {'x': '2/3', 'q': '6/5', 'unit_norm': '5/9', 'next_level': '12/13'})
        self.assertAlmostEqual(res['inward']['width'], 0.12266457577602966, places=15)
        self.assertLess(res['outward']['width_400'], 1e-6)

    def test_algebra(self):
        self.assertEqual(y.mul(y.R, y.R), y.scal(-1, y.ONE))
        self.assertEqual(y.mul(y.K, y.K), y.ONE)
        self.assertEqual(y.mul(y.L, y.L), y.ONE)
        self.assertEqual(y.mul(y.R, y.K), y.scal(-1, y.L))

    def test_finite_chain_matches_r1_recursion(self):
        q, d = F(6, 5), F(0)
        for k in range(1, 9):
            d = q / (1 + q*d)
            self.assertEqual(y.chain_return([q]*k), d)
        lo, hi = sorted((y.chain_return([q]*40), y.chain_return([q]*41)))
        self.assertTrue(lo < F(2, 3) < hi and hi - lo < F(1, 10**12))

    def test_directional_amplitudes_are_not_seen(self):
        Ft = y.el(1, F(1, 4))
        self.assertEqual(y.append_cell(Ft, F(3), F(3), F(1), F(1)), y.append_cell(Ft, F(3), F(3), F(7, 2), F(1, 9)))

    def test_unpaired_cell_leaves_the_line(self):
        self.assertNotEqual(y.append_cell(y.el(1, F(1, 4)), F(3), F(2))[0], 1)

    def test_horizon_closure_changes_the_outer_reading(self):
        xs = [F(1, 3)]
        for _ in range(7):
            xs.append(y.tower(xs[-1]))
        c = y.cells_of_profile(xs)
        wall, null = y.chain_return(c, F(0)), y.chain_return(c, F(1))
        self.assertGreater(wall - null, F(6, 100))
        self.assertLess(abs(null - F(1, 3)), F(1, 10**30))

    def test_flat_end_closure_does_not(self):
        c = y.cells_of_profile([F(1, j + 3) for j in range(402)])
        self.assertLess(abs(y.chain_return(c, F(0)) - y.chain_return(c, F(5))), F(1, 10**6))

    def test_wrong_tower_law_is_rejected(self):
        with patch.object(y, 'tower', lambda x: (x + 1)/2):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
