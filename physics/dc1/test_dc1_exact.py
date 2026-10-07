"""DC1: stdlib only."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import dc1_memory_between_masses as y


class MemoryBetweenMassesTests(unittest.TestCase):
    def test_exchange_fixes_the_relative_unit(self):
        rows = y.unit_control()
        self.assertEqual(rows[0]['bracket'], '17/2')

    def test_a_scalar_readout_cannot_fix_it(self):
        with patch.object(y, 'trace', lambda a: y.c(0)):
            with self.assertRaises(ValueError):
                y.unit_control()

    def test_pair_memory_law(self):
        row = y.memory_control()
        self.assertEqual(row['cases'], 27)
        self.assertEqual(row['isolated_or_no_phase_memory'], '0')
        self.assertEqual(row['example']['memory'], '165888/1953125')

    def test_phase_on_every_branch_pair_leaves_no_memory(self):
        def global_phase(a, b, zeta):
            return [[y.cmul(y.cmul(y.c(a[i]), y.c(b[j])), zeta) for j in range(2)] for i in range(2)]
        with patch.object(y, 'pair_state', global_phase):
            with self.assertRaises(ValueError):
                y.memory_control()

    def test_illustration(self):
        ill = y.phase_illustration()
        self.assertAlmostEqual(ill['clock_factor_phase_rate'], ill['direct_phase_rate'], places=9)
        self.assertAlmostEqual(ill['owner_example']['rate_scale_per_s']*1e4, 6.33, places=1)
        self.assertLess(ill['owner_example']['rows'][1]['pair_memory'], ill['owner_example']['rows'][1]['exponential_law_deficit']/50)


if __name__ == '__main__':
    unittest.main()
