"""CA2 exact focused tests."""
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

import ca2_symmetry_compass_atlas as m


class CA2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_packet_and_pins(self):
        self.assertTrue(self.result["all_pass"])
        self.assertTrue(all(self.result["checks"].values()))
        stored = json.loads((Path(__file__).parent / "CA2_RESULT.json").read_text())
        self.assertEqual(stored, self.result)

    def test_electromagnetic_diagram_and_ratios(self):
        em = m.electromagnetic(F(4), F(1), F(3), F(8), F(9))
        self.assertEqual(em["pre_observation_diagram"], ("E", "B", "D", "H"))
        self.assertEqual(em["post_observation_chart"], ("B", "E", "D", "H"))
        self.assertEqual(em["chi_EM"], F(11, 12))
        self.assertEqual(em["C_B"]/em["C_H"], em["chi_EM"])
        self.assertEqual(em["C_D"]/em["C_E"], em["chi_EM"])
        self.assertEqual(em["epsilon_B"]/em["epsilon_H"], em["chi_EM"])
        self.assertEqual(em["mu_D"]/em["mu_E"], em["chi_EM"])

    def test_vacuum_and_magnetoelectric_separate(self):
        vacuum = m.electromagnetic(F(2), 0, F(7), F(3), F(5))
        coupled = m.electromagnetic(F(2), 1, F(7), F(3), F(5))
        self.assertEqual(vacuum["chi_EM"], 1)
        self.assertEqual(coupled["chi_EM"], F(13, 14))
        self.assertNotEqual(vacuum["chi_EM"], coupled["chi_EM"])

    def test_thermodynamic_instance_retains_four_definitions(self):
        generic = m.scale_coefficients(F(4), F(2), F(3), F(10), F(21))
        thermo = m.ca1.capacities(F(4), F(2), F(3), F(10), F(21), F(7))
        self.assertEqual(generic["x_at_y"], thermo["C_V"])
        self.assertEqual(generic["x_at_q"], thermo["C_P"])
        self.assertEqual(generic["y_at_x"], thermo["C_S"])
        self.assertEqual(generic["y_at_p"], thermo["C_T"])

    def test_source_pair_is_part_of_the_address(self):
        G = [[F(4), F(1), F(0)], [F(1), F(3), F(1)], [F(0), F(1), F(2)]]
        self.assertEqual(m.principal_compass(G, 0, 1), F(11, 12))
        self.assertEqual(m.principal_compass(G, 1, 2), F(5, 6))
        self.assertNotEqual(m.principal_compass(G, 0, 1), m.principal_compass(G, 1, 2))

    def test_covariance_and_elastic_instances(self):
        cov = m.covariance_compass(F(9), F(3), F(4))
        self.assertEqual(cov["chi_cov"], F(3, 4))
        self.assertEqual(cov["squared_correlation"], F(1, 4))
        elastic = m.elastic(F(6), F(2), F(5), F(3), F(4))
        self.assertEqual(elastic["chi_elastic"], F(13, 15))

    def test_nonreciprocal_response_retains_turn(self):
        turn = m.ca1.noncommuting_compass(F(4), F(-1, 2), F(3, 2), F(2))
        self.assertGreater(turn["chi"], 1)
        self.assertGreater(turn["omega_squared"], 0)

    def test_invalid_adapters_are_refused(self):
        for args in ((1, 2, 1), (0, 0, 1), (1, 1, 1)):
            with self.assertRaises(ValueError):
                m.scale_coefficients(*args)
        with self.assertRaises(ValueError):
            m.scale_coefficients(2, 1, 2, p_scale=-1)
        with self.assertRaises(ValueError):
            m.principal_compass([[1, 2], [3, 4]], 0, 1)
        with self.assertRaises(ValueError):
            m.principal_compass([[1, 0], [0, 1]], 0, 0)


if __name__ == "__main__":
    unittest.main()
