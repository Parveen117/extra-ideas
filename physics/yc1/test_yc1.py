"""YC1 exact centre-symmetry, response-compass and degree-filter tests."""
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

import mpmath as mp

import yc1_centre_compass_potential as m


class YC1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_packet_and_source_pins(self):
        self.assertTrue(self.result["all_pass"])
        self.assertTrue(all(self.result["checks"].values()))
        stored = json.loads((Path(__file__).parent / "YC1_RESULT.json").read_text())
        self.assertEqual(stored, self.result)

    def test_centre_hessian_and_compass(self):
        hessian, chi = m.centre_hessian()
        self.assertEqual(hessian, ((F(1, 4), F(0)), (F(0), F(1, 16))))
        self.assertEqual(chi, 1)

    def test_st_scalar_and_diagram_centre_series(self):
        s = m.response_series(12)
        self.assertEqual(s["st"][2], F(1, 4))
        self.assertEqual(s["st"][4], F(-1, 96))
        self.assertEqual(s["centre"][:4], [F(0)]*4)
        self.assertEqual(s["centre"][4], F(1, 384))
        self.assertEqual(s["centre"][6], F(-1, 4608))

    def test_compass_loss_begins_quadratically_off_the_seam(self):
        s = m.response_series(12)
        self.assertEqual(s["lost"][:5], [F(0), F(0), F(1, 4), F(0), F(-1, 16)])
        self.assertEqual(s["chi"][:5], [F(1), F(0), F(-1, 4), F(0), F(1, 16)])

    def test_bessel_response_is_potential_derivative(self):
        mp.mp.dps = 50
        for kappa in (mp.mpf("0.2"), mp.mpf("1.5"), mp.mpf("12")):
            self.assertAlmostEqual(m.face_response(kappa), mp.diff(m.face_potential, kappa), places=40)

    def test_rw1_global_bound_implies_monotone_centre_cut(self):
        mp.mp.dps = 50
        for kappa in (mp.mpf("0.001"), mp.mpf("0.5"), mp.mpf("5"), mp.mpf("100")):
            r = m.face_response(kappa)
            mid_prime = (4*r-kappa*(1-r*r))/2
            self.assertGreater(mid_prime, 0)
            self.assertLess(mid_prime, r/2)

    def test_global_st_scalar_sandwich(self):
        mp.mp.dps = 50
        for kappa in (mp.mpf("0.01"), mp.mpf("0.3"), mp.mpf("2"), mp.mpf("30")):
            psi = m.face_potential(kappa)
            st = kappa*m.face_response(kappa)
            mid = m.centre_cut(kappa)
            self.assertGreater(mid, 0)
            self.assertLess(mid, psi/2)
            self.assertGreater(st, psi)
            self.assertLess(st, 2*psi)

    def test_source_and_observation_positions_are_typed(self):
        tvsp = self.result["tvsp"]
        self.assertEqual(tvsp["pre_observation"], "(T,V,S,P) = (Psi_kappa,h,kappa,Psi_h)")
        self.assertEqual(tvsp["post_observation"], "(V,T,S,P) = (h,Psi_kappa,kappa,Psi_h)")

    def test_claim_boundary_rejects_unproved_physical_identifications(self):
        boundary = self.result["claim_boundary"]
        self.assertTrue(boundary["gauge_invariant_two_source_response_potential_constructed"])
        self.assertTrue(boundary["first_retained_face_term_and_TC1_core_share_quartic_degree"])
        self.assertFalse(boundary["one_face_kappa_identified_with_CM2_theta_or_continuum_coupling"])
        self.assertFalse(boundary["one_face_quartic_coefficient_identified_with_TC1_e2_core"])
        self.assertFalse(boundary["mass_gap_or_volume_uniformity_proved"])

    def test_invalid_series_and_moment_inputs_are_refused(self):
        with self.assertRaises(ValueError):
            m.response_series(3)
        with self.assertRaises(ValueError):
            m.haar_moment(-1)
        with self.assertRaises(ValueError):
            m.sdiv([F(1)], [F(0)], 3)


if __name__ == "__main__":
    unittest.main()
