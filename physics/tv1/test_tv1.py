"""TV1: stdlib only."""
from fractions import Fraction as F
from math import comb, cos, exp, log, pi, sin, sqrt
import unittest

import tv1_valley_of_turns as m


class ValleyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()

    def test_all(self):
        self.assertGreaterEqual(len(self.out), 13)
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_closed_paths_by_hand(self):
        self.assertEqual(m.closed_paths(1, 16)[:6], [comb(2*k, k) for k in range(6)])
        self.assertEqual(m.closed_paths(2, 16)[:6], [comb(2*k, k)**2 for k in range(6)])
        self.assertEqual(m.closed_paths(3, 16)[:3], [1, 6, 90])        # 90 = 6 x 6 + 6 x 4 x 2 + 6: back-and-forth pairs
        self.assertEqual(m.closed_paths_by_walking(2, 4), 36)

    def test_determinant_by_hand(self):
        # d = 3 at x = (1, 2, 3): [[1(2+3), -2], [-2, 2(1+3)]] = [[5, -2], [-2, 8]], determinant 36 = 1.2.3.(1+2+3)
        det = m.pdet(m.graph_form(3), 3)
        self.assertEqual(sum(c*m._prod([F(1), F(2), F(3)], e) for e, c in det.items()), 36)
        self.assertEqual(m.perm_sign((1, 0, 2)), -1)
        self.assertEqual(m.perm_sign((1, 2, 0)), 1)

    def test_split_by_hand(self):
        n = (F(0), F(0), F(1))
        u, v = (F(1), F(0), F(2)), (F(0), F(1), F(3))                  # p = 2, 3 ; t = (1,0,0), (0,1,0)
        x = m.dot(m.cross(u, v), m.cross(u, v))                        # u x v = (-2, -3, 1) : 14
        self.assertEqual(x, 14)
        layer = m.sub(m.scale(F(2), (F(0), F(1), F(0))), m.scale(F(3), (F(1), F(0), F(0))))
        self.assertEqual(m.dot(layer, layer) + 1, 14)                  # 13 from the layer, 1 from the core

    def test_path_expansion_against_a_grid(self):
        # mean of 1/(1 + sum sin^2 a_i), three turns: smooth, so a grid mean is accurate
        g = 24
        grid = sum(1/(1 + sin(pi*(i + .5)/g)**2 + sin(pi*(j + .5)/g)**2 + sin(pi*(k + .5)/g)**2)
                   for i in range(g) for j in range(g) for k in range(g))/g**3
        c3 = m.closed_paths(3)
        series = (2/5)*sum((3/5)**(2*k)*c3[k]/6**(2*k) for k in range(120))    # 1 + (3/2)(1 - y) = (5/2)(1 - (3/5) y)
        self.assertAlmostEqual(grid, series, places=10)

    def test_against_the_classical_values(self):
        # cross-check only (classical, not used): Watson's integral through the pair at heat time sqrt 6,
        # and Log 2/pi^2 per doubling for four turns
        q = exp(-pi*sqrt(6))
        s = (1 + 2*sum(q**(k*k) for k in range(1, 8)))**2
        watson = (18 + 12*sqrt(2) - 10*sqrt(3) - 7*sqrt(6))*s*s
        self.assertAlmostEqual(watson, 0.5054620197, places=9)
        last = self.num['three turns: partial sums at 2m = 16 .. 1024'][-1]
        tail = (2/3)*2*(3/(4*pi))**1.5*2/sqrt(512)                      # rest of the sum beyond 2m = 1024
        self.assertAlmostEqual(last + tail, 2*watson, places=4)
        self.assertLess(last, 2*watson)
        steps = self.num['four turns: added per doubling']
        self.assertAlmostEqual(steps[-1], log(2)/pi**2, places=3)


if __name__ == '__main__':
    unittest.main()
