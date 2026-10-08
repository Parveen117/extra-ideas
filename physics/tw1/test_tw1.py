"""TW1: stdlib only."""
from fractions import Fraction as F
import math
import unittest

import tw1_turn_block_walk as m


def floats(u):
    """plain floating sums: a, b, K+, K- at Exp(-pi u); the mirror description below u = 0.3 (cross-check only)"""
    if u >= 0.3:
        q = math.exp(-math.pi*u)
        a = 1 + 2*sum(q**(n*n) for n in range(1, 60))
        b = 1 + 2*sum((-1)**n*q**(n*n) for n in range(1, 60))
        kp = sum(n*n*q**(n*n) for n in range(1, 60))
        km = sum((-1)**(n - 1)*n*n*q**(n*n) for n in range(1, 60))
        return a, b, kp, km
    v = 1/u
    t3 = 1 + 2*sum(math.exp(-math.pi*v*n*n) for n in range(1, 30) if math.pi*v*n*n < 700)
    kpv = sum(n*n*math.exp(-math.pi*v*n*n) for n in range(1, 30) if math.pi*v*n*n < 700)
    t2 = 2*sum(math.exp(-math.pi*v*(j + .5)**2) for j in range(30) if math.pi*v*(j + .5)**2 < 700)
    mm = sum((j + .5)**2*math.exp(-math.pi*v*(j + .5)**2) for j in range(30) if math.pi*v*(j + .5)**2 < 700)
    return t3/math.sqrt(u), t2/math.sqrt(u), u**-1.5*t3/(4*math.pi) - u**-2.5*kpv, u**-2.5*mm - u**-1.5*t2/(4*math.pi)


class TurnBlockWalkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()

    def test_all(self):
        self.assertGreaterEqual(len(self.out), 33)
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_first_terms_by_hand(self):
        a, b, kp, km, od = m.series(30)
        self.assertEqual(kp, {1: 1, 4: 4, 9: 9, 16: 16, 25: 25})
        self.assertEqual(km, {1: 1, 4: -4, 9: 9, 16: -16, 25: 25})
        self.assertEqual(od, {1: 1, 9: 9, 25: 25})
        self.assertEqual(a, {0: 1, 1: 2, 4: 2, 9: 2, 16: 2, 25: 2})
        # up to q^5 :  (1 + 2q + 2q^4)(q + 4q^4) = q + 2q^2 + 4q^4 + 10q^5 ,  (1 - 2q + 2q^4)(q - 4q^4) = q - 2q^2 - 4q^4 + 10q^5
        left = m.sadd(m.smul(a, kp, 5), m.smul(b, km, 5), -1)
        self.assertEqual(left, {2: 4, 4: 8})                              # = 4 (1 + 2q^2) q^2 : 4 a(q^2) K+(q^2)
        for order in (4, 8, 30):
            a, b, kp, km, _ = m.series(order)
            self.assertEqual(m.sscale(m.smul(m.dil(a, 2, order), m.dil(kp, 2, order), order), 4),
                             m.sadd(m.smul(a, kp, order), m.smul(b, km, order), -1))

    def test_lost_parts_by_hand(self):
        # L_0^2/16 = q + 4q^3 + 6q^5 + ... ,  L_1^2/16 = q^2 + 4q^6 ,  L_2^2/16 = q^4 ,  b = 1 - 2q + 2q^4
        # (q + 2q^2 + 4q^3 + 4q^4)(1 - 2q + 2q^4) = q - 4q^4 + ...   = K- to order 4
        a, b, kp, km, _ = m.series(4)
        four = lambda x: m.smul(m.smul(x, x, 4), m.smul(x, x, 4), 4)
        l0 = m.sadd(four(a), four(b), -1)
        self.assertEqual(l0, {1: 16, 3: 64})
        total = m.sadd(m.sadd(l0, m.dil(l0, 2, 4), 2), m.dil(l0, 4, 4), 4)
        self.assertEqual(total, {1: 16, 2: 32, 3: 64, 4: 64})
        self.assertEqual(m.smul(b, total, 4), m.sscale(km, 16))

    def test_lost_parts_against_floats(self):
        for u in (0.08, 0.2, 0.5, 1.0):
            lost = []
            for n in range(12):
                a, b, _, _ = floats(u*2**n)
                lost.append(2**n*(a**4 - b**4))
            a, b, kp, km = floats(u)
            self.assertAlmostEqual(16*km/(b*sum(lost)), 1, places=9)
            self.assertAlmostEqual(16*kp/(a*(lost[0] - sum(lost[1:]))), 1, places=6)

    def test_rest_is_an_upper_bound(self):
        for q in (F(1, 5), F(1, 2), F(7, 10)):
            for power in (0, 2):
                more = sum(n**power*q**(n*n) for n in range(m.TERMS + 1, m.TERMS + 12))
                self.assertGreater(m.rest(q, power), more)
                self.assertLess(m.rest(q, power), 2*more)

    def test_loss(self):
        self.assertGreater(m.loss(F(1, 10), F(5)), 1)
        self.assertLess(m.loss(F(1, 10), F(5)), F(1001, 1000))        # delta_1 = 2 (1/200)^2 5 (1 + ...) : tiny
        self.assertGreater(m.loss(F(2, 5), F(5)), m.loss(F(1, 5), F(5)))
        with self.assertRaises(AssertionError):
            m.loss(F(1), F(5))

    def test_half_step_solves_the_mean_step(self):
        r = (F(2, 5), F(2, 5))
        x, _ = m.half_step(r, (F(3), F(3)))
        for end, sign in ((x[0], -1), (x[1], 1)):
            self.assertLess(abs(2*end/(1 + end*end) - F(4, 25)), F(1, 10**100))
            self.assertGreaterEqual(sign*(2*end/(1 + end*end) - F(4, 25)), 0)

    def test_box_against_floats(self):
        for q1, q2 in ((F(1, 4), F(26, 100)), (F(48, 100), F(49, 100))):
            r, beta = m.box(q1, q2)
            for q in (q1, (q1 + q2)/2, q2):
                a, b, kp, km = floats(-math.log(float(q))/math.pi)
                self.assertTrue(float(r[0]) - 1e-12 <= b/a <= float(r[1]) + 1e-12)
                self.assertTrue(float(beta[0]) - 1e-9 <= km*a/(kp*b) <= float(beta[1]) + 1e-9)

    def test_inequality_against_floats(self):
        # the two sides of the product form on a fine grid of the strong side
        for i in range(20, 700):
            q = i/1000
            k = lambda x: sum(n*n*x**(n*n - 1) for n in range(1, 40))
            o = lambda x: sum(n*n*x**(n*n - 1) for n in range(1, 40, 2))
            self.assertGreater(k(q*q)**2*k(q**4)**2, o(q)**2*o(q*q)*k(q**8))
        # the inequality itself and beta on the weak side
        for i in range(5, 115):
            u = i/1000
            a, b, kp, km = floats(u)
            _, _, kp2, km2 = floats(2*u)
            rho, rho2 = km/kp, km2/kp2
            self.assertGreater(rho2*(1 + rho), 2*math.sqrt(rho))
            self.assertGreater(rho*a/b, 9.9)
        a, b, kp, km = floats(-math.log(0.7)/math.pi)
        self.assertGreater(km*a/(kp*b), self.num['smallest certified beta for q >= 7/10'])

    def test_step_law_against_floats(self):
        for u in (0.05, 0.11, 0.2, 0.4, 0.8):
            a, b, kp, km = floats(u)
            a2, b2, kp2, km2 = floats(2*u)
            self.assertAlmostEqual(4*a2*kp2/(a*kp - b*km), 1, places=9)
            self.assertAlmostEqual(4*b2*km2/(a*km - b*kp), 1, places=6)


if __name__ == '__main__':
    unittest.main()
