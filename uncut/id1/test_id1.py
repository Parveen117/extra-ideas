"""ID1: sympy, exact."""
import unittest
from unittest.mock import patch

import sympy as sp

import id1_information_dimension_curvature as y


class InformationDimensionCurvatureTests(unittest.TestCase):
    def test_run(self):
        res = y.run()
        self.assertEqual(res['sharing_among_directions'][3], 'radial : each transverse = 1 : -1/2')
        self.assertEqual(res['spreading'][3]['far_value'], '0')
        self.assertEqual(res['spreading'][1]['far_value'], 'oo')

    def test_a_null_reading_shared_by_three_cuts(self):
        n, rs = 3, (1, 2, 2)                                   # 9 = 1 + 4 + 4
        self.assertEqual(sum(q*q for q in rs), n*n)
        self.assertEqual(sum(sp.Rational(q*q, n*n) for q in rs), 1)
        self.assertEqual(1 - sp.Rational(rs[0]**2, n*n), sp.Rational(rs[1]**2 + rs[2]**2, n*n))

    def test_two_sources_add_and_stay_free(self):
        xs = y.coords(3)
        a = sp.Matrix([1, 0, 0])
        m = 1/y.radius(xs) + 2/sp.sqrt((xs[0] - a[0])**2 + xs[1]**2 + xs[2]**2)
        self.assertEqual(sp.simplify(sum(sp.diff(m, q, 2) for q in xs)), 0)

    def test_I_itself_is_not_free_to_add(self):
        xs = y.coords(3)
        rad = y.radius(xs)
        I1 = -sp.log(1 - 1/rad)/2
        lap = sp.simplify(sum(sp.diff(I1, q, 2) for q in xs))
        self.assertNotEqual(lap, 0)

    def test_a_law_that_ignores_the_dimension_is_rejected(self):
        with patch.object(y, 'spread', lambda d, rad: 1/rad):
            with self.assertRaises(ValueError):
                y.run()


if __name__ == '__main__':
    unittest.main()
