"""TS1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import ts1_two_sheet_closure as y
rb = y.rb


class TwoSheetClosureTests(unittest.TestCase):
    def test_run_and_pinned_values(self):
        res = y.run()
        self.assertAlmostEqual(res['width_one_sheet'], 0.12266457577602966, places=12)
        self.assertAlmostEqual(res['width_two_sheets']['0'], 0.02319519064044732, places=12)
        self.assertLess(res['width_two_sheets']['20'], 1e-40)

    def test_any_closure_behind_the_second_sheet_gives_the_point_mass_value(self):
        cells = rb.cells_of_profile(y.two_sheet_profile(F(1, 3), 6, 40))
        for tail in (F(0), F(1), F(7), F(10**9)):
            self.assertLess(abs(rb.chain_return(cells, tail) - F(1, 3)), F(1, 10**60))

    def test_one_sheet_closures_differ(self):
        cells = rb.cells_of_profile(y.inward(F(1, 3), 6))
        self.assertGreater(abs(rb.chain_return(cells, F(0)) - rb.chain_return(cells, F(7))), F(1, 10))
        odd = rb.cells_of_profile(y.inward(F(1, 3), 7))
        self.assertGreater(rb.chain_return(odd, F(0)), F(39, 100))     # wall value depends on the parity of the cell count
        self.assertLess(rb.chain_return(cells, F(0)), F(28, 100))

    def test_other_starting_radius(self):
        cells = rb.cells_of_profile(y.inward(F(1, 10), 7) + y.inward(F(1, 10), 7)[::-1] + y.flat_end(60, start=11))
        self.assertLess(abs(rb.chain_return(cells, F(0)) - F(1, 10)), F(1, 10**40))

    def test_a_tower_without_the_mirror_is_rejected(self):
        with patch.object(rb, 'tower', lambda x: (3*x + 1)/4):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
