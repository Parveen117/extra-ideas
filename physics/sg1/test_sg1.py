"""Independent ODE, angular-selection and fail-closed checks for SG1."""
from fractions import Fraction as F
import math
import unittest

import sg1_singlet_gap as m


def rk4_shoot(energy, radius=12.0, step=1/2048):
    # A floating cross-check only, independent of the certified Taylor series.
    y, v, zeros = 0.0, 1.0, 0
    def rhs(x, a, b):
        return b, (x-energy)*a
    for i in range(round(radius/step)):
        x = i*step
        k1 = rhs(x, y, v)
        k2 = rhs(x+step/2, y+step*k1[0]/2, v+step*k1[1]/2)
        k3 = rhs(x+step/2, y+step*k2[0]/2, v+step*k2[1]/2)
        k4 = rhs(x+step, y+step*k3[0], v+step*k3[1])
        yn = y+step*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        vn = v+step*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        zeros += int(y*yn < 0)
        y, v = yn, vn
    return y, v, zeros


class SingletGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = m.run()

    def test_all_outward_checks(self):
        self.assertGreaterEqual(len(self.result['checks']), 18)
        for name, passed in self.result['checks'].items():
            self.assertTrue(passed, name)

    def test_shooting_with_independent_runge_kutta(self):
        for row in self.result['radial_certificates']:
            y, v, zeros = rk4_shoot(float(F(row['energy'])))
            self.assertEqual(zeros, row['dirichlet_negative_count'])
            self.assertEqual([1 if a > 0 else -1 for a in (y, v)], row['endpoint_signs'])
            for value, enclosure in zip((y, v), row['endpoint_intervals']):
                mid = sum(float(F(x)) for x in enclosure)/2
                self.assertLess(abs(value-mid)/abs(mid), 1e-5)

    def test_sturm_counts_survive_different_mesh_and_order(self):
        for row in self.result['radial_certificates'][:2]:
            other = m.radial_count(F(row['energy']), step=F(1, 16), order=28, bits=144)
            self.assertEqual(other['negative_count'], row['negative_count'])
            for old, new in zip(row['endpoint_intervals'], other['endpoint_intervals']):
                self.assertLessEqual(max(F(old[0]), F(new[0])), min(F(old[1]), F(new[1])))

    def test_neumann_count_needs_endpoint_correction(self):
        # E=2.3 on [0,2.5] has no Dirichlet eigenvalue below it, but has one
        # Dirichlet-Neumann eigenvalue. A root count alone would be wrong.
        row = m.radial_count(F(23, 10), radius=F(5, 2))
        self.assertEqual(row['dirichlet_negative_count'], 0)
        self.assertEqual(row['negative_count'], 1)
        self.assertEqual(row['endpoint_signs'], [1, -1])

    def test_outward_rounding_for_both_signs(self):
        for a in [(F(-7, 11), F(5, 13)), (F(-1, 3), F(-1, 7)), (F(2, 9), F(2, 7))]:
            lo, hi = m.round_out(a, 8)
            self.assertLessEqual(lo, a[0])
            self.assertGreaterEqual(hi, a[1])
            self.assertLess(a[0]-lo, F(1, 256))
            self.assertLess(hi-a[1], F(1, 256))

    def test_angular_trial_energy_from_factorial_moments(self):
        # u=r^2 exp(-a r): integrate derivative and centrifugal terms separately.
        a = F(27, 25)
        moment = lambda n: F(math.factorial(n))/(2*a)**(n+1)
        kinetic = (4*moment(2)-4*a*moment(3)+a*a*moment(4))/moment(4)
        centrifugal = 2*moment(2)/moment(4)
        potential = moment(5)/moment(4)
        self.assertEqual(kinetic+centrifugal+potential, a*a+F(5, 2)/a)

    def test_gauge_selection_examples_include_higher_spins(self):
        self.assertNotIn(0, m.coupled_spins([0, 73, 0, 0]))
        self.assertNotIn(0, m.coupled_spins([1, 2, 0]))
        self.assertIn(0, m.coupled_spins([2, 2, 0]))
        self.assertIn(0, m.coupled_spins([1, 1, 1]))
        self.assertEqual(m.coupled_spins([0, 0, 0]), {0})

    def test_unresolved_or_invalid_certificate_fails_closed(self):
        with self.assertRaises(ValueError):
            m.radial_count(F(3), radius=F(2))
        with self.assertRaises(ValueError):
            m.radial_count(F(23381, 10000), order=2)
        with self.assertRaises(ValueError):
            m.taylor_cell(F(0), F(1), (m.iv(0), m.iv(1)), step=F(1))

    def test_exact_gap_and_scope(self):
        self.assertEqual(F(self.result['numbers']['3']['core_gap_lower']), F(1671, 5000))
        self.assertEqual(F(self.result['numbers']['4']['core_gap_lower']), F(3, 1000))
        self.assertEqual(self.result['claim_boundary']['continuum_mass_gap'], 'OPEN')
        self.assertEqual(self.result['claim_boundary']['lowest_excited_sector_identity'], 'OPEN')


if __name__ == '__main__':
    unittest.main()
