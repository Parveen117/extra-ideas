"""Independent local SU(2) derivative and geometry controls for YC8."""
from fractions import Fraction as F
from itertools import combinations
import unittest
import yc8_local_vacuum_dressing as y


def qm(p, q):
    a,b,c,d = p
    w,x,z,t = q
    return (a*w-b*x-c*z-d*t, a*x+b*w+c*t-d*z,
            a*z-b*t+c*w+d*x, a*t+b*z-c*x+d*w)


def qc(q):
    return (q[0],-q[1],-q[2],-q[3])


def trace_word(word, qs):
    out = (F(1),F(0),F(0),F(0))
    for edge, sign in word:
        out = qm(out, qs[edge] if sign == 1 else qc(qs[edge]))
    return out[0]


def first(word, qs, edge, axis):
    if edge not in {e for e,s in word}:
        return F(0)
    unit = tuple(F(int(i == axis)) for i in range(4))
    moved = list(qs)
    moved[edge] = qm(unit, qs[edge])
    return trace_word(word, moved)


def coefficients(word, qs, edge):
    out = []
    for axis in range(4):
        moved = list(qs)
        moved[edge] = tuple(F(int(i == axis)) for i in range(4))
        out.append(trace_word(word, moved))
    return out


class DressingTests(unittest.TestCase):
    def test_shared_link_harmonics_with_both_orientations(self):
        qs = [(F(1,2),)*4, (F(3,5),F(4,5),0,0), (0,F(3,5),F(4,5),0),
              (F(3,5),0,0,F(4,5)), (F(1,2),F(-1,2),F(1,2),F(-1,2)),
              (0,0,F(3,5),F(4,5)), (F(4,5),0,F(3,5),0)]
        p = ((0,1),(1,1),(2,-1),(3,-1))
        for sign in (-1,1):
            q = ((0,sign),(4,1),(5,-1),(6,-1))
            wp,wq = trace_word(p,qs),trace_word(q,qs)
            a,b = coefficients(p,qs,0),coefficients(q,qs,0)
            A,B = sum(x*z for x,z in zip(a,b)),wp*wq
            self.assertEqual(sum(x*x for x in a),1)
            self.assertEqual(sum(x*x for x in b),1)
            gamma = sum(first(p,qs,0,k)*first(q,qs,0,k) for k in (1,2,3))
            self.assertEqual(gamma,A-B)
            hB = F(0)
            for e in range(7):
                for k in (1,2,3):
                    ddp = -wp if e in {x for x,s in p} else F(0)
                    ddq = -wq if e in {x for x,s in q} else F(0)
                    hB -= ddp*wq+2*first(p,qs,e,k)*first(q,qs,e,k)+wp*ddq
            self.assertEqual(hB,26*B-2*A)
            self.assertEqual(F(2,39)*18*A-hB/26,gamma)

    def test_single_plaquette_casimir_and_gradient_norm(self):
        qs = [(F(1,2),)*4, (F(3,5),F(4,5),0,0),
              (F(3,5),0,F(4,5),0),(F(3,5),0,0,F(4,5))]
        word = ((0,1),(1,1),(2,-1),(3,-1))
        w = trace_word(word,qs)
        gamma = F(0)
        hchi = F(0)
        for e in range(4):
            local = sum(first(word,qs,e,k)**2 for k in (1,2,3))
            self.assertEqual(local,1-w*w)
            gamma += local
            for k in (1,2,3):
                hchi -= 8*(first(word,qs,e,k)**2-w*w)
        self.assertEqual(gamma,3-(4*w*w-1))
        self.assertEqual(hchi,32*(4*w*w-1))

    def test_closed_plaquette_words(self):
        g = y.geometry((3,4,5))
        for word in g['words']:
            balance = {}
            for e, sign in word:
                x, direction = g['edges'][e]
                end = list(x)
                end[direction] = (end[direction]+1) % (3,4,5)[direction]
                end = tuple(end)
                balance[x] = balance.get(x,0)+sign
                balance[end] = balance.get(end,0)-sign
            self.assertTrue(all(v == 0 for v in balance.values()))

    def test_unique_shared_links(self):
        g = y.geometry((3,3,4))
        for p,q in combinations(range(len(g['faces'])),2):
            shared = g['faces'][p] & g['faces'][q]
            self.assertLessEqual(len(shared),1)
            self.assertEqual(bool(shared),(p,q) in g['pairs'])

    def test_stencil_size_and_transpose_incidence(self):
        for shape in ((3,3,3),(4,4,4),(5,5,5)):
            g = y.geometry(shape)
            columns = [0]*len(g['edges'])
            for e, stencil in enumerate(g['stencils']):
                self.assertIn(e,stencil)
                self.assertLessEqual(len(stencil),121)
                for f in stencil:
                    columns[f] += 1
                    self.assertIn(e,g['stencils'][f])
            self.assertLessEqual(max(columns),121)

    def test_local_gradient_and_hessian_constants(self):
        g,h = y.local_constants()
        self.assertEqual(g,F(1,288)+F(1,39)+F(1,39))
        self.assertEqual(h,F(7,144)+F(2,13)+F(49,78))
        self.assertEqual(g,y.GRAD_F2)
        self.assertEqual(h,y.HESS_F2)

    def test_exact_cancellation(self):
        checks = y.symbolic_checks()
        self.assertTrue(all(checks.values()),checks)

    def test_uniform_density_bounds_have_correct_second_order(self):
        for theta in (F(1,1000),F(1,100),F(1,10),F(1,4),F(1,2)):
            lo,hi = y.energy_density(theta)
            centre = theta-theta*theta/48
            self.assertLessEqual(lo,centre)
            self.assertGreaterEqual(hi,centre)
            self.assertEqual(centre-lo,y.residual_density(theta))
            self.assertLessEqual(hi,theta)

    def test_parent_gap_is_decreasing_and_positive_on_its_window(self):
        values = [y.comparison_gap(F(i,20)) for i in range(11)]
        self.assertTrue(all(a>b for a,b in zip(values,values[1:])))
        self.assertEqual(values[-1],F(941,3744))
        self.assertGreater(values[-1],F(1,4))

    def test_density_bound_does_not_imply_volume_uniform_actual_gap(self):
        theta = F(1,10)
        eps = y.residual_density(theta)
        # Ordinary norm perturbation loses 2 P epsilon, despite fixed density.
        p = 10**6
        self.assertLess(y.comparison_gap(theta)-2*p*eps,0)
        self.assertGreater(y.comparison_gap(theta),1)

    def test_collapsed_geometries_and_bad_parameters_are_rejected(self):
        for shape in ((2,3,3),(1,1,1),(3,3),(True,3,3),(3.0,3,3)):
            with self.assertRaises(ValueError):
                y.geometry(shape)
        for t in (-1,F(3,4)):
            with self.assertRaises(ValueError):
                y.comparison_gap(t)
        with self.assertRaises(ValueError):
            y.energy_density(-1)


if __name__ == '__main__':
    unittest.main()
