"""Tests for OL1.  python -m unittest test_ol1"""
from fractions import Fraction as F
from math import cos, exp, pi, sin, sqrt
import cmath
import json
import os
import random
import unittest

import ol1_one_site_lattice as ol


def lowest_two(diag, off):
    """Two lowest eigenvalues of a symmetric tridiagonal matrix, by counting sign changes (floats)."""
    def below(x):
        count, d = 0, 1.0
        for i in range(len(diag)):
            d = diag[i] - x - (off[i - 1]**2/d if i else 0.0)
            if d == 0.0:
                d = 1e-300
            if d < 0:
                count += 1
        return count
    out = []
    for want in (1, 2):
        lo, hi = min(diag) - 2*max(abs(t) for t in off) - 1, max(diag) + 1
        for _ in range(80):
            mid = (lo + hi)/2
            if below(mid) >= want:
                hi = mid
            else:
                lo = mid
        out.append(hi)
    return out


def class_even_levels(theta, n=4000):
    """-u''/2 - u/2 + V u on [0, pi/2], u(0) = 0, free end, by finite differences (floats)."""
    h = (pi/2)/n
    diag = [1/h**2 - 0.5 + sqrt(4 + 4*theta*sin(i*h)**2) - 2 for i in range(1, n + 1)]
    diag[-1] -= 0.5/h**2                              # free end
    return lowest_two(diag, [-0.5/h**2]*(n - 1))


class TestOL1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out, cls.num = ol.run()

    def test_all(self):
        bad = [k for k, v in self.out['detail'].items() if not v]
        self.assertEqual(bad, [])
        self.assertEqual(self.out['checks'], 19)

    def test_plaquette_with_matrices(self):
        rnd = random.Random(5)

        def link():
            psi = rnd.uniform(0, pi)
            m = [rnd.gauss(0, 1) for _ in range(3)]
            r = sqrt(sum(x*x for x in m))
            u = [sin(psi)*x/r for x in m]
            return cos(psi), u

        def matrix(c, u):
            return [[complex(c, u[2]), complex(u[1], u[0])], [complex(-u[1], u[0]), complex(c, -u[2])]]

        def mul(a, b):
            return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

        def dag(a):
            return [[a[j][i].conjugate() for j in range(2)] for i in range(2)]
        for _ in range(50):
            (c1, u), (c2, v) = link(), link()
            a, b = matrix(c1, u), matrix(c2, v)
            w = mul(mul(a, b), mul(dag(a), dag(b)))
            cross = [u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]]
            self.assertAlmostEqual((w[0][0] + w[1][1]).real/2, 1 - 2*sum(x*x for x in cross), places=12)

    def test_layer_on_the_sphere(self):
        # the Laplacian of the sphere as the flat Laplacian of the degree-zero extension, by differences
        kappa, h = 0.8, 1e-3

        def big(x):
            n2 = sum(t*t for t in x)
            return exp(-kappa*(x[1]**2 + x[2]**2)/n2)
        for x in ([0.5, 0.5, 0.5, 0.5], [0.1, 0.7, 0.1, sqrt(1 - 0.51)], [0.6, 0.0, 0.8, 0.0]):
            lap = 0.0
            for i in range(4):
                up, dn = list(x), list(x)
                up[i] += h
                dn[i] -= h
                lap += (big(up) - 2*big(x) + big(dn))/(h*h)
            t = x[1]**2 + x[2]**2
            self.assertAlmostEqual(lap/big(x), 4*(kappa*kappa*t*(1 - t) - kappa*(1 - 2*t)), places=4)
        # the layer's true lowest level (Legendre form: 4 n (n + 1) + mu/2 on the diagonal) against the floor
        for mu in (0.5, 2.0, 10.0, 50.0, 400.0):
            n = 120
            diag = [4*i*(i + 1) + mu/2 for i in range(n)]
            off = [-(mu/2)*(i + 1)/sqrt((2*i + 1)*(2*i + 3)) for i in range(n - 1)]
            level = lowest_two(diag, off)[0]
            floor = float(ol.layer_floor(F(mu)))
            self.assertLessEqual(floor, level + 1e-9)
            self.assertLess(level - floor, 2.3)

    def test_local_rates_on_one_link(self):
        kappa, h = 0.45, 1e-4

        def rate(f, psi, turning):
            d1 = (f(psi + h) - f(psi - h))/(2*h)
            d2 = (f(psi + h) - 2*f(psi) + f(psi - h))/(h*h)
            return -0.5*(d2 + 2*cos(psi)/sin(psi)*d1 - turning*f(psi)/sin(psi)**2)/f(psi)
        for psi in (0.2, 0.7, 1.3, pi/2, 2.4):
            y = sin(psi)**2
            self.assertAlmostEqual(rate(lambda p: exp(-kappa*sin(p)**2), psi, 0),
                                   3*kappa - 4*kappa*y - 2*kappa**2*y*(1 - y), places=5)
            self.assertAlmostEqual(rate(lambda p: sin(p)*exp(-kappa*sin(p)**2), psi, 2),
                                   1.5 + 5*kappa - 6*kappa*y - 2*kappa**2*y*(1 - y), places=5)

    def test_readings_by_quadrature(self):
        n = 20000
        grid = [(i + 0.5)*pi/n for i in range(n)]

        def means(f, df):
            norm = sum(f(p)**2*sin(p)**2 for p in grid)
            return (sum(df(p)**2*sin(p)**2 for p in grid)/norm, sum(f(p)**2*sin(p)**4 for p in grid)/norm)
        b = 0.3
        k, s = means(lambda p: 1 + b*cos(2*p), lambda p: -2*b*sin(2*p))
        self.assertAlmostEqual(3*k + 4*2.0*s*s, float(ol.two_term(F(2), F(3, 10))), places=8)
        for m in (1, 3, 7):
            k, s = means(lambda p: cos(p)**(2*m), lambda p: -2*m*cos(p)**(2*m - 1)*sin(p))
            self.assertAlmostEqual(3*k + 4*5.0*s*s, float(ol.power_reading(F(5), m)), places=7)

    def test_floors_against_a_direct_solve(self):
        # strong side: the floor of the first level is under the true level, and close to it
        for theta in (0.5, 2.0, 3.5):
            e0, e1 = class_even_levels(theta, 2000)
            th = F(theta)
            floor = float(ol.class_floor(th, ol.best_kappa(th, lambda x: 3*x, lambda x: 4*x)))
            self.assertLess(floor, e0)
            self.assertLess(e0 - floor, 0.05)
            self.assertGreater(e1, 4.0)
        # weak side: the two levels of the comparison at theta = 10^7 against the floors by the half line
        theta = 1e7
        e0, e1 = class_even_levels(theta, 6000)
        a = float(ol.CUT)
        unit = (2*theta*(sin(a)/a)**2)**(1/3)
        self.assertGreater(e0, unit*float(ol.E0) - 2.5)
        self.assertGreater(e1, unit*float(ol.E1) - 2.5)
        self.assertLess(e0, (2*theta)**(1/3)*2.3382)          # and the half line is the right scale
        weak = self.num['weak side']
        self.assertGreater(2*e0 + e1, weak['spots_theta_n_line_reading_gap_over_cube_root'][0][2])

    def test_numbers(self):
        self.assertEqual(self.num['strong side']['window'], [0.0, 3.75])
        self.assertGreater(self.num['strong side']['least_margin'], 0.3)
        weak = self.num['weak side']
        self.assertEqual((weak['start'], weak['gap_coefficient'], weak['limit_coefficient']), (10**7, 0.1314, 0.3271))
        self.assertEqual(self.num['in the units of the YM line']['weak_start_theta'], 2500000.0)
        with open(os.path.join(os.path.dirname(os.path.abspath(ol.__file__)), 'OL1_RESULT.json')) as f:
            self.assertTrue(json.load(f)['all_pass'])


if __name__ == '__main__':
    unittest.main()
