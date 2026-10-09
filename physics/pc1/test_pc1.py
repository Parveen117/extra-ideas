"""Tests for PC1.  python -m unittest test_pc1"""
from fractions import Fraction as F
from math import sqrt
import itertools
import json
import os
import random
import unittest

import pc1_pair_clusters as pc


def jacobi_levels(a, sweeps=60):
    """Eigenvalues of a symmetric matrix by plane rotations (floats)."""
    n = len(a)
    a = [row[:] for row in a]
    for _ in range(sweeps):
        off = sum(a[i][j]**2 for i in range(n) for j in range(i))
        if off < 1e-22:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-15:
                    continue
                t = (a[q][q] - a[p][p])/(2*a[p][q])
                t = (1 if t >= 0 else -1)/(abs(t) + sqrt(t*t + 1))
                c = 1/sqrt(t*t + 1)
                s = t*c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c*akp - s*akq, s*akp + c*akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c*apk - s*aqk, s*apk + c*aqk
    return sorted(a[i][i] for i in range(n))


class TestPC1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = pc.run()

    def test_all(self):
        bad = [k for k, v in self.out['detail'].items() if not v]
        self.assertEqual(bad, [])
        self.assertEqual(self.out['checks'], 15)

    def test_cluster_inequality_on_small_models(self):
        # pair operator on 3 x 3 points with a symmetric lowest reading g (rate e) and everything else at
        # least rho ; three turns: H = sum over pairs ; its two lowest levels against the bound
        rnd = random.Random(3)
        m = 3
        for trial in range(4):
            g = [[0.0]*m for _ in range(m)]
            for i in range(m):
                for j in range(i, m):
                    g[i][j] = g[j][i] = rnd.uniform(-0.3, 0.3) + (2.0 if i == j == 0 else 0.0)
            norm = sqrt(sum(x*x for row in g for x in row))
            g = [[x/norm for x in row] for row in g]
            e, rho = 1.0, 2.0 + rnd.random()
            flat = [g[i][j] for i in range(m) for j in range(m)]
            pair = [[rho*(1 if p == q else 0) - (rho - e)*flat[p]*flat[q] for q in range(m*m)] for p in range(m*m)]
            pts = list(itertools.product(range(m), repeat=3))
            big = [[0.0]*len(pts) for _ in pts]
            for a, x in enumerate(pts):
                for b, y in enumerate(pts):
                    for (i, j, k) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
                        if x[k] == y[k]:
                            big[a][b] += pair[x[i]*m + x[j]][y[i]*m + y[j]]
            levels = jacobi_levels(big)
            red = [[sum(g[i][k]*g[j][k] for k in range(m)) for j in range(m)] for i in range(m)]
            r = sorted(jacobi_levels(red), reverse=True)
            self.assertAlmostEqual(sum(r), 1.0, places=9)
            # this model has exactly the levels 3 rho - (rho - e)(1 + 2 r_n), ... ; the bound is attained
            self.assertAlmostEqual(levels[0], 3*rho - (rho - e)*(1 + 2*r[0]), places=7)
            self.assertGreaterEqual(levels[1] + 1e-9, 3*rho - (rho - e)*(1 + 2*r[1]))
            self.assertGreaterEqual(levels[1] + 1e-9, 2*rho + e - 2*(1 - r[0])*(rho - e))

    def test_direction_average_of_the_product_reading(self):
        # two turns: the mean of |c_1|^2 |c_2|^2 f equals the mean of (e1^2/8 + e2/2) f for f of e1, e2
        tc = pc.cr1.tc1
        e1, e2, e3 = tc.entry_polys(2)
        var = lambda k: {tuple(1 if i == k else 0 for i in range(6)): F(1)}
        sq = lambda cols: tc.padd(tc.padd(tc.pmul(var(cols[0]), var(cols[0])), tc.pmul(var(cols[1]), var(cols[1]))),
                                  tc.pmul(var(cols[2]), var(cols[2])))
        # entries are ordered row by row: (a, i) -> a*2 + i
        n1, n2 = sq((0, 2, 4)), sq((1, 3, 5))
        self.assertEqual(tc.padd(n1, n2), e1)
        lhs = tc.pmul(n1, n2)
        rhs = tc.padd({k: v/8 for k, v in tc.pmul(e1, e1).items()}, {k: v/2 for k, v in e2.items()})
        for f in ({(0,)*6: F(1)}, e1, e2, tc.pmul(e1, e2)):
            self.assertEqual(tc.free_mean(tc.pmul(lhs, f)), tc.free_mean(tc.pmul(rhs, f)))

    def test_share_and_line(self):
        r = dict(eta=F(26, 10), overlap2=F(96, 100))
        share, p0 = pc.one_turn_share(r, F(32, 10), F(259, 100))
        self.assertTrue(0 < share < F(96, 100) and 0 < p0 < 1)
        # more share, a higher line ; a higher floor, a higher line
        self.assertLess(pc.cluster_line(3, F(32, 10), F(26, 10), F(9, 10)), pc.cluster_line(3, F(32, 10), F(26, 10), F(95, 100)))
        self.assertLess(pc.cluster_line(3, F(32, 10), F(26, 10), F(9, 10)), pc.cluster_line(3, F(33, 10), F(26, 10), F(9, 10)))
        # with share one the line is kappa (2 rho + E): three turns from two
        self.assertAlmostEqual(float(pc.cluster_line(3, F(3), F(2), F(1))), 2**(-2/3)*8, places=6)

    def test_numbers(self):
        two, three, four = self.num['two turns'], self.num['three turns'], self.num['four turns']
        self.assertEqual(two['lowest_rate'], [2.6592, 2.6594])
        self.assertEqual(two['gap_all_readings'], 0.5536)
        self.assertGreater(three['line'], 5.69)
        self.assertGreater(three['gap'][0], 0.5)
        self.assertGreater(four['gap_all_readings'], 0.46)
        with open(os.path.join(os.path.dirname(os.path.abspath(pc.__file__)), 'PC1_RESULT.json')) as f:
            self.assertTrue(json.load(f)['all_pass'])


if __name__ == '__main__':
    unittest.main()
