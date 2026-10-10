"""YC2 exact compact-centre response tests."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import unittest

import yc2_compact_centre_response as m


class YC2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()
        cls.data = m.haar_energy_data()

    def test_packet_and_source_pins(self):
        self.assertTrue(self.result["all_pass"])
        self.assertTrue(all(self.result["checks"].values()))
        stored = json.loads((Path(__file__).parent / "YC2_RESULT.json").read_text())
        self.assertEqual(stored, self.result)

    def test_cross_and_gram_forms_agree(self):
        for edge in ((0, 1), (0, 2), (1, 2)):
            self.assertEqual(m.energy_polynomial(*edge), m.gram_energy_polynomial(*edge))
            self.assertEqual(m.link_parities(m.energy_polynomial(*edge)), (0, 0, 0))

    def test_exact_haar_energy_moments(self):
        self.assertEqual(self.result["haar_energy_response"]["centre_odd_loop_covariance"], [
            ["1/4", "0", "0"], ["0", "1/4", "0"], ["0", "0", "1/4"]])
        self.assertEqual(self.data["means"], [F(3, 8)]*3)
        self.assertEqual(self.data["covariance"], [
            [F(13, 192), F(1, 64), F(1, 64)],
            [F(1, 64), F(13, 192), F(1, 64)],
            [F(1, 64), F(1, 64), F(13, 192)],
        ])

    def test_shared_link_compass_and_modes(self):
        self.assertEqual(self.data["pair_chi"], F(160, 169))
        self.assertEqual(self.data["pair_lost"], F(9, 169))
        self.assertEqual(self.data["symmetric_eigenvalue"], F(19, 192))
        self.assertEqual(self.data["anisotropy_eigenvalue"], F(5, 96))
        self.assertGreater(self.data["anisotropy_eigenvalue"], 0)

    def test_total_wilson_core_channel(self):
        self.assertEqual(self.data["mean_Q"], F(9, 8))
        self.assertEqual(self.data["variance_Q"], F(19, 64))
        self.assertEqual(self.data["normalized_determinant"], F(1900, 2197))

    def test_character_fourier_transform_is_invertible(self):
        signs = list(product((-1, 1), repeat=3))
        table = [[m.character(sigma, sheet) for sheet in signs] for sigma in signs]
        for i, row in enumerate(table):
            for j, other in enumerate(table):
                self.assertEqual(sum(a*b for a, b in zip(row, other)), 8 if i == j else 0)

    def test_observation_pushes_the_typed_chart(self):
        before, after = m.observation_chart("odd")
        self.assertEqual(before, ("<c_i>", "kappa_j", "kappa_i", "<c_j>"))
        self.assertEqual(after, ("kappa_j", "<c_i>", "kappa_i", "<c_j>"))
        before, after = m.observation_chart("energy")
        self.assertEqual(after, (before[1], before[0], before[2], before[3]))

    def test_claim_boundary_selects_twisted_heat_trace(self):
        boundary = self.result["claim_boundary"]
        self.assertTrue(boundary["centre_odd_and_even_response_blocks_separate_at_the_symmetric_seam"])
        self.assertFalse(boundary["centre_even_Wilson_potential_alone_can_identify_a_centre_character"])
        self.assertTrue(boundary["twisted_heat_trace_is_required_for_sector_splitting"])
        self.assertFalse(boundary["absolute_CM2_sector_splitting_rate_computed"])

    def test_invalid_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            m.energy_polynomial(0, 0)
        with self.assertRaises(ValueError):
            m.energy_polynomial(0, 3)
        with self.assertRaises(ValueError):
            m.sphere_moment((2, -2, 0))
        with self.assertRaises(ValueError):
            m.permute_links(m.energy_polynomial(0, 1), (0, 0, 2))
        with self.assertRaises(ValueError):
            m.character((1, 1, 0), (1, 1, 1))
        with self.assertRaises(ValueError):
            m.observation_chart("unknown")


if __name__ == "__main__":
    unittest.main()
