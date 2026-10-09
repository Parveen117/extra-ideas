"""HL1: stdlib only."""
from fractions import Fraction as F
import math
import random
import unittest

import hl1_helical_walk as m


class HelicalWalkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()

    def test_all(self):
        self.assertGreaterEqual(len(self.out), 30)
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_digits_by_hand(self):
        self.assertEqual(m.digit((1, 1)), (0, 1, (1, 0)))
        self.assertEqual(m.digit((2, 0)), (3, 2, (1, 0)))                 # 2 = iota^3 (1 + iota)^2
        self.assertEqual(m.digit((3, 0)), (2, 0, (-3, 0)))                # 3 = iota^2 (-3) ; -3 = 1 mod 4
        self.assertEqual(m.digit((2, 1)), (3, 0, (-1, 2)))                # 2 + iota = iota^3 (-1 + 2 iota)
        self.assertEqual(m.digit((0, 4)), (3, 4, (1, 0)))                 # 4 iota = iota^3 (1 + iota)^4 , (1 + iota)^4 = -4
        self.assertEqual(m.times((1, 1), (1, 1)), (0, 2))
        self.assertTrue(m.is_primary((1, 0)) and m.is_primary((-3, 0)) and m.is_primary((-1, 2)))
        self.assertFalse(m.is_primary((3, 0)) or m.is_primary((1, 2)) or m.is_primary((2, 1)))

    def test_mean(self):
        one = m.mean(m.iv(1), m.iv(1))
        self.assertTrue(one[0] <= 1 <= one[1] and one[1] - one[0] < F(1, 10**100))
        g = m.mean(m.iv(24), m.iv(6))
        self.assertAlmostEqual(float(g[0]), 13.458171481725615, places=12)   # the mean of 24 and 6
        self.assertLess(g[1] - g[0], F(1, 10**100))
        self.assertTrue(6 < g[0] and g[1] < 24)

    def test_clock_against_floats(self):
        for u in (0.3, 0.7, 1.0, 2.5):
            q = math.exp(-math.pi*u)
            a = 1 + 2*sum(q**(n*n) for n in range(1, 40))
            l = (2*sum(q**((n + .5)**2) for n in range(40)))**2
            x, y = a*a, l
            for _ in range(40):
                x, y = (x + y)/2, math.sqrt(x*y)
            self.assertAlmostEqual(u*x, 1, places=12)

    def test_sectors_against_floats(self):
        q = 0.3
        b = 1 + 2*sum((-1)**n*q**(n*n) for n in range(1, 30))
        km = sum((-1)**(n - 1)*n*n*q**(n*n) for n in range(1, 30))
        series = sum(2**m.split2(n)[0]*m.sigma(m.split2(n)[1])*q**n for n in range(1, 200))
        self.assertAlmostEqual(km/b, series, places=12)
        self.assertEqual([m.split2(n) for n in (1, 2, 12, 40)], [(0, 1), (1, 1), (2, 3), (3, 5)])
        self.assertEqual([m.sigma(n) for n in (1, 3, 9, 15)], [1, 4, 13, 24])

    def test_sphere_moments(self):
        self.assertEqual(m.sphere_moment((2, 0, 0, 0)), F(1, 4))
        self.assertEqual(m.sphere_moment((4, 0, 0, 0)), F(1, 8))
        self.assertEqual(m.sphere_moment((2, 2, 0, 0)), F(1, 24))
        self.assertEqual(m.sphere_moment((1, 1, 0, 0)), 0)
        self.assertEqual(4*F(1, 8) + 12*F(1, 24), 1)                      # mean of (x0^2 + ... + x3^2)^2
        self.assertEqual(m.cheb_u(2), [F(-1), F(0), F(4)])                # U_2 = 4c^2 - 1
        self.assertEqual(m.x_moment(0), 1)

    def test_commutator_by_hand(self):
        i, j = (0, 1, 0, 0), (0, 0, 1, 0)
        comm = m.qmul(m.qmul(i, j), m.qmul(m.qconj(i), m.qconj(j)))
        self.assertEqual(comm, (-1, 0, 0, 0))                             # i j = -j i : the twisted pair
        self.assertEqual(1 - 2*m.cross_sq(i, j), -1)

    def test_three_turns_against_random_blocks(self):
        rnd = random.Random(7)

        def block():
            while True:
                v = [rnd.gauss(0, 1) for _ in range(4)]
                n = math.sqrt(sum(x*x for x in v))
                if n > 1e-9:
                    return tuple(x/n for x in v)

        def c(u, v):
            return m.qmul(m.qmul(u, v), m.qmul(m.qconj(u), m.qconj(v)))[0]
        n, t1, t2, t3 = 60000, 0.0, 0.0, 0.0
        for _ in range(n):
            u1, u2, u3 = block(), block(), block()
            c12, c23, c31 = c(u1, u2), c(u2, u3), c(u3, u1)
            t1 += c12
            t2 += c12*c23
            t3 += c12*c23*c31
        self.assertAlmostEqual(t1/n, 1/4, delta=0.01)
        self.assertAlmostEqual(t2/n, 1/8, delta=0.01)
        self.assertAlmostEqual(t3/n, 5/72, delta=0.01)

    def test_residue_against_floats(self):
        u = 0.25
        v = 1/u
        t3 = 1 + 2*sum(math.exp(-math.pi*v*n*n) for n in range(1, 10))
        t4 = 1 + 2*sum((-1)**n*math.exp(-math.pi*v*n*n) for n in range(1, 10))
        kpv = sum(n*n*math.exp(-math.pi*v*n*n) for n in range(1, 10))
        delta = (1 - 4*math.pi*v*kpv/t3)/t4**4 - 1                         # the mirror description of the residue
        row = [r for r in self.num['u, residue, first sheet term, power of u it is below'] if r[0] == '1/4'][0]
        self.assertAlmostEqual(row[1]/delta, 1, places=9)


if __name__ == '__main__':
    unittest.main()
