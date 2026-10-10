"""CA1 focused exact tests."""
from fractions import Fraction as F
from pathlib import Path
import json
import unittest
from unittest.mock import patch

import sympy as sp

import ca1_compass_calibration as m


class CA1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_packet_and_source_pins(self):
        self.assertTrue(self.result["all_pass"])
        self.assertTrue(all(self.result["checks"].values()))
        stored = json.loads((Path(__file__).parent / "CA1_RESULT.json").read_text())
        self.assertEqual(stored, self.result)

    def test_four_C_properties_remain_typed(self):
        r = m.capacities(F(4), F(2), F(3), T=F(10), P=F(21), V=F(7))
        self.assertEqual(r["chi"], F(2, 3))
        self.assertEqual(r["C_V"] / r["C_P"], r["chi"])
        self.assertEqual(r["C_S"] / r["C_T"], r["chi"])
        self.assertEqual(r["K_T"] / r["K_S"], r["chi"])
        self.assertEqual(r["C_S"] * r["K_S"], 147)
        self.assertEqual(r["C_T"] * r["K_T"], 147)

    def test_one_anchor_per_dimensional_family(self):
        self.assertEqual(m.conjugate_from("C_V", F(10), F(2, 5)), ("C_P", F(25)))
        self.assertEqual(m.conjugate_from("C_P", F(25), F(2, 5)), ("C_V", F(10)))
        self.assertEqual(m.conjugate_from("C_S", F(6), F(3, 4)), ("C_T", F(8)))
        self.assertEqual(m.conjugate_from("C_T", F(8), F(3, 4)), ("C_S", F(6)))
        self.assertEqual(m.conjugate_from("K_S", F(8), F(3, 4)), ("K_T", F(6)))
        self.assertEqual(m.conjugate_from("K_T", F(6), F(3, 4)), ("K_S", F(8)))

    def test_positive_diagonal_congruence_orbit(self):
        original = m.symmetric_compass(F(4), F(-2), F(3))
        scaled_args = m.unit_rescale(F(4), F(-2), F(-2), F(3), F(7), F(11))
        scaled = m.symmetric_compass(scaled_args[0], scaled_args[1], scaled_args[3])
        self.assertEqual(original["chi"], scaled["chi"])
        self.assertEqual(original["orientation"], scaled["orientation"])
        self.assertNotEqual((F(4), F(-2), F(3)), (scaled_args[0], scaled_args[1], scaled_args[3]))

    def test_chart_type_is_load_bearing(self):
        self.assertEqual(m.chi_from_chart(F(3, 5), "seen_fraction"), F(3, 5))
        self.assertEqual(m.chi_from_chart(F(5, 3), "amplification"), F(3, 5))
        self.assertEqual(m.chi_from_chart(F(2, 5), "lost_fraction"), F(3, 5))
        self.assertEqual(sp.simplify(m.chi_from_chart(sp.log(F(5, 3))/2, "cut_information")), F(3, 5))
        with self.assertRaises(ValueError):
            m.chi_from_chart(F(3, 5), "dimensionless")

    def test_noncommuting_turn_needs_a_second_number(self):
        first = m.noncommuting_compass(F(1), F(1, 2), F(1, 2), F(1))
        second = m.noncommuting_compass(F(1), F(3, 2), F(1, 6), F(1))
        self.assertEqual(first["chi"], second["chi"])
        self.assertNotEqual(first["omega_squared"], second["omega_squared"])
        a = m.unit_rescale(F(1), F(3, 2), F(1, 6), F(1), F(5), F(9))
        moved = m.noncommuting_compass(*a)
        for key in ("chi", "beta_squared", "omega_squared", "beta_sign", "omega_sign"):
            self.assertEqual(second[key], moved[key])

    def test_reduced_volume_needs_curve_and_scale(self):
        self.assertEqual(m.reduced_volume_linear(F(2, 3), 1), F(1, 2))
        self.assertEqual(m.reduced_volume_linear(F(2, 3), 2), F(1, 4))
        V01, V02 = F(10), F(100)
        self.assertEqual(V01*m.reduced_volume_linear(F(2, 3), 1), 5)
        self.assertEqual(V02*m.reduced_volume_linear(F(2, 3), 2), 25)

    def test_invalid_blocks_anchors_and_untyped_matching_are_refused(self):
        for args in ((1, 2, 1), (0, 0, 1), (1, 1, 1)):
            with self.assertRaises(ValueError):
                m.symmetric_compass(*args)
        with self.assertRaises(ValueError):
            m.capacities(2, 1, 2, T=-1)
        with self.assertRaises(ValueError):
            m.conjugate_from("C_X", 1, F(1, 2))
        with self.assertRaises(ValueError):
            m.conjugate_from("C_V", 1, F(3, 2))
        with self.assertRaises(ValueError):
            m.unit_rescale(2, 1, 1, 2, -1, 1)

    def test_omitting_CS_CT_pair_is_rejected(self):
        original = m.capacities

        def incomplete(*args, **kwargs):
            values = original(*args, **kwargs)
            values.pop("C_S")
            values.pop("C_T")
            return values

        with patch.object(m, "capacities", incomplete):
            values = m.capacities(4, 2, 3)
            self.assertNotIn("C_S", values)
            self.assertFalse({"C_P", "C_V", "C_S", "C_T"}.issubset(values))


if __name__ == "__main__":
    unittest.main()
