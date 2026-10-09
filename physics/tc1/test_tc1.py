"""TC1: stdlib only."""
from fractions import Fraction as F
from math import exp, log, sqrt
import random
import unittest

import tc1_core_of_turns as m


class CoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = m.run()

    def test_all(self):
        self.assertGreaterEqual(len(self.out), 19)
        for k, v in self.out.items():
            self.assertTrue(v, k)

    def test_weight_by_hand(self):
        two = [[F(1), F(0)], [F(0), F(1)], [F(0), F(0)]]               # two directions at right angles: |c1 x c2|^2 = 1
        self.assertEqual(m.weight(two), 1)
        self.assertEqual(m.e2_of(m.gram(two)), 1)                       # G = diag(1, 1, 0)
        same = [[F(2), F(4)], [F(0), F(0)], [F(1), F(2)]]              # parallel: the valley
        self.assertEqual(m.weight(same), 0)
        self.assertEqual(m.det([[F(2), F(1)], [F(1), F(3)]]), 5)
        r = m.rotation(3, random.Random(3))
        self.assertEqual(m.mat_mul(r, m.transpose(r)), [[F(int(i == j)) for j in range(3)] for i in range(3)])

    def test_free_means_against_random_matrices(self):
        rnd = random.Random(11)
        n, tot = 20000, 0.0
        for _ in range(n):
            c = [[rnd.gauss(0, sqrt(.5)) for _ in range(4)] for _ in range(3)]
            cols = list(zip(*c))
            tot += sum(m.dot(m.cross(cols[i], cols[j]), m.cross(cols[i], cols[j]))
                       for i in range(4) for j in range(i + 1, 4))
        self.assertAlmostEqual(tot/n, 9, delta=0.2)                     # four turns: mean of e2 is 9

    def test_ray_series_by_hand(self):
        s = m.ray_series()
        self.assertEqual(s[-1], F(1, 16))        # x1^2 . t :  1!/(2 x1)^2 . 1/(4 x1) = 1/(16 x1)
        self.assertEqual(sorted(s), [-7, -5, -3, -1])

    def test_logarithm_against_a_grid(self):
        # K(4) - K(2) by a plain grid (cross-check only): it must lie in [(Log 2)/16 - delta(2), (Log 2)/16]
        n1, n2 = 40, 160
        total = 0.0
        for i in range(n1):
            x1 = 2 + (i + .5)*2/n1
            top = min(x1, 8/x1)                                         # beyond this the weight is below e^-16
            h = top/n2
            inner = 0.0
            for j in range(n2):
                x2 = (j + .5)*h
                hk = x2/n2
                for k in range(n2):
                    x3 = (k + .5)*hk
                    inner += (x1 - x2)*(x1 - x3)*(x2 - x3)*exp(-2*(x1*x2 + x1*x3 + x2*x3))*h*hk
            total += inner*2/n1
        upper = log(2)/16
        lower = upper - self.num['delta(Lambda), Lambda = 2 .. 32'][0]
        self.assertTrue(lower - 1e-4 <= total <= upper + 1e-4, (lower, total, upper))

    def test_scale_against_a_search(self):
        for d in (3, 4):
            a, b = 3*d/4, 1.5*d*(d - 1)/2
            best = min(a/(l*l) + b*l**4 for l in [0.3 + 0.0005*i for i in range(3000)])
            self.assertAlmostEqual(best, (729/256*d**3*(d - 1))**(1/3), places=5)


if __name__ == '__main__':
    unittest.main()
