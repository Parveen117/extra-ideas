"""Independent arithmetic, metric, derivative and false-transfer controls for CB1."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import unittest

import cb1_compact_bridge as c


def decimal(x):
    return Decimal(x.numerator)/Decimal(x.denominator)


class CompactBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = c.run()

    def test_whole_packet(self):
        self.assertEqual(len(self.result['checks']), 18)
        self.assertTrue(all(self.result['checks'].values()))

    def test_heat_bounds_against_independent_exponentials(self):
        with localcontext() as ctx:
            ctx.prec = 90
            for even in (False, True):
                packet = c.heat_deviation_upper(even)
                time = decimal(packet['time'])
                ns = range(2, 202, 2) if even else range(1, 201)
                direct = sum((Decimal((n+1)**2)*(-time*n*(n+2)).exp() for n in ns), Decimal(0))
                self.assertLess(direct, decimal(packet['deviation_upper']))
                self.assertLess(decimal(packet['deviation_upper'])-direct, Decimal('1e-10'))

    def test_sharp_two_state_contraction_and_missing_barrier_control(self):
        # K=[[1,b],[b,1]]: its nonvacuum/vacuum eigenvalue ratio is exactly (1-b)/(1+b).
        for beta in (F(1, 1000000), F(1, 5), F(7, 8)):
            contraction = (1-beta)/(1+beta)
            with localcontext() as ctx:
                ctx.prec = 60
                self.assertGreaterEqual(-decimal(contraction).ln(), 2*decimal(beta))
            # Every putative beta larger than the actual kernel ratio overstates the gap.
            wrong_beta = (1+beta)/2
            self.assertGreater(contraction, (1-wrong_beta)/(1+wrong_beta))
        # Nonconstant vacuum: eigenvalue ratio of [[a,b],[b,d]] obeys the same bound.
        for a, b, d in ((F(4), F(1), F(2)), (F(3), F(1, 2), F(1))):
            beta = min(a, b, d)/max(a, b, d)
            discriminant = (a-d)**2+4*b*b
            self.assertGreaterEqual(discriminant, beta*beta*(a+d)**2)

    def test_weighted_kinetic_from_pointwise_gradient_polynomial(self):
        p = c.flat_trial()
        poly, square, w = p['p'], p['square'], p['width']
        # L=Delta/2 gives |grad p|^2=L(p^2)-2 p L(p).
        gradient = c.cr.padd(c.cr.gen(square, 3), c.cr.pmul(poly, c.cr.gen(poly, 3)), -2)
        euler = {e: c.cr.degree(e)*v for e, v in poly.items()}
        density = {e: v/2 for e, v in gradient.items()}
        density = c.cr.padd(density, c.cr.pmul(poly, euler), -2*w)
        density = c.cr.padd(density, c.cr.pmul(c.cm.E1, square), w*w/2)
        self.assertEqual(p['mean'](density)/p['norm'], p['kinetic'])
        weighted = {(a+1, b, cc): v for (a, b, cc), v in density.items()}
        self.assertEqual(p['mean'](weighted)/p['norm'], p['weighted_kinetic'])
        self.assertNotEqual(p['mean'](weighted)/p['norm'], p['weighted_kinetic']-F(9, 2))

    def test_compact_metric_and_haar_are_not_flat(self):
        u = [F(1, 3), F(1, 4), F(1, 5)]
        radius2 = sum(x*x for x in u)
        metric = [[F(i == j)+u[i]*u[j]/(1-radius2) for j in range(3)] for i in range(3)]
        inverse = [[F(i == j)-u[i]*u[j] for j in range(3)] for i in range(3)]
        self.assertEqual(c.cm.matmul(metric, inverse), [[F(i == j) for j in range(3)] for i in range(3)])
        determinant = (metric[0][0]*(metric[1][1]*metric[2][2]-metric[1][2]*metric[2][1])
            -metric[0][1]*(metric[1][0]*metric[2][2]-metric[1][2]*metric[2][0])
            +metric[0][2]*(metric[1][0]*metric[2][1]-metric[1][1]*metric[2][0]))
        self.assertEqual(determinant, 1/(1-radius2))
        self.assertNotEqual(determinant, 1)

    def test_invalid_lift_and_incomplete_comparison_rejected(self):
        with self.assertRaises(ValueError):
            c.compact_upper(1, F(5), F(1), 20, F(1, 100))
        with self.assertRaises(ValueError):
            c.compact_upper(10**7, F(1, 10), F(1, 4), 20, F(1, 100))
        with self.assertRaises(ValueError):
            c.weak_row(10**3, F(39, 100), F(5), F(1, 4), 20, F(1, 6000), F(36, 100))
        with self.assertRaises(ValueError):
            c.return_reserve(3, F(1))
        with self.assertRaises(ValueError):
            c.sin_ratio_lower(F(2))

    def test_lift_improves_with_coupling_but_global_reserve_decays_with_volume(self):
        a = c.compact_upper(2*10**6, F(19, 4), F(1, 4), 18, F(442, 10**6))
        b = c.compact_upper(10**7, F(19, 4), F(1, 4), 18, F(442, 10**6))
        self.assertLess(b['upper'], a['upper'])
        self.assertGreater(a['upper'], c.flat_trial()['eta'])
        self.assertLess(c.return_reserve(100, F(1, 4)), c.return_reserve(3, F(1, 4)))
        self.assertFalse(self.result['claim_boundary']['volume_uniformity'])
        self.assertFalse(self.result['claim_boundary']['continuum_mass_gap'])

    def test_exact_root_enclosures(self):
        for x in (F(0), F(2), F(1, 2000000), F(2000000)):
            for power in (2, 3):
                lo, hi = c.root_floor(x, power), c.root_upper(x, power)
                self.assertLessEqual(lo**power, x)
                self.assertGreaterEqual(hi**power, x)
                self.assertLessEqual(hi-lo, F(1, 10**12))


if __name__ == '__main__':
    unittest.main()
