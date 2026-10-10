"""Independent instruments, state sensitivity and field normalization checks."""
import hashlib
import json
import unittest
import sympy as sp
import yc23_cut_agitation as y


class YC23Tests(unittest.TestCase):
    def test_complex_sharp_cut_and_balanced_minimum(self):
        h, om = sp.diag(0, 7), sp.Matrix([1, 0])
        k = sp.Matrix([[0, -sp.I], [sp.I, 0]])
        e, q, direct, double, second = y.sharp_data(h, om, k)
        self.assertEqual((e, q), (sp.Rational(7, 2), sp.Rational(1, 2)))
        self.assertEqual(e, direct)
        self.assertEqual(e, double)
        self.assertEqual(second, 49)

    def test_noninvolutive_compressed_reading_is_not_a_sharp_cut(self):
        with self.assertRaises(ValueError):
            y.sharp_data(sp.diag(0, 2), sp.Matrix([1, 0]), sp.diag(1, 0))

    def test_global_noncommutation_can_be_silent_on_vacuum(self):
        h = sp.diag(0, 2, 9)
        k = sp.diag(1, sp.Matrix([[0, 1], [1, 0]]))
        self.assertNotEqual(y.comm(h, k), sp.zeros(3))
        self.assertEqual(y.sharp_data(h, sp.Matrix([1, 0, 0]), k)[0], 0)

    def test_smooth_balanced_instrument_with_nonconstant_true_ground(self):
        omega, h, edges = y.graph_ground()
        f = (sp.Rational(1), sp.Rational(0), -sp.Rational(1, 4))
        eta = sp.Rational(3, 5)
        ms = [sp.diag(*(sp.sqrt((1+sign*eta*x)/2) for x in f)) for sign in (1, -1)]
        rho = omega*omega.T
        post = sum((m*rho*m for m in ms), sp.zeros(3))
        self.assertTrue(y.zero(sp.trace(post)-1))
        for m in ms:
            self.assertTrue(y.zero((omega.T*m*m*omega)[0]-sp.Rational(1, 2)))
        energy = sp.trace(h*post)
        edge_energy = sum(c*sum((m[i, i]-m[j, j])**2 for m in ms)
                          for (i, j), c in edges.items())
        self.assertTrue(y.zero(energy-edge_energy))
        # The retained measure is omega², not uniform counting measure.
        self.assertNotEqual(sum(x*x for x in f)/3, sum(omega[i]**2*f[i]**2 for i in range(3)))

    def test_metric_can_be_nonzero_with_both_transport_curvatures_zero(self):
        x, z = sp.symbols('x z', real=True)
        a = x*z+x*x
        w = sp.Matrix([sp.cos(a/2), sp.sin(a/2)])
        g = sp.Matrix([[(w.diff(i).T*w.diff(j))[0] for j in (x, z)] for i in (x, z)])
        self.assertTrue(y.zero(g-g.T))
        self.assertTrue(y.zero(w.T*w.diff(x)))
        self.assertTrue(y.zero(w.T*w.diff(z)))
        self.assertEqual(sp.simplify(g[0, 0].subs({x: 1, z: 1})), sp.Rational(9, 4))
        k = 2*w*w.T-sp.eye(2)
        r = sp.Matrix([[0, -1], [1, 0]])
        self.assertTrue(y.zero((r*k)**2-sp.eye(2)))

    def test_haar_moments_from_independent_sphere_moment_recurrence(self):
        # Normalized S^3 coordinate recurrence: E[c^(2n)]/E[c^(2n-2)] = (2n-1)/(2n+2).
        moment = sp.S.One
        for n in range(1, 14):
            moment *= sp.Rational(2*n-1, 2*n+2)
            self.assertEqual(moment, y.haar_moment(2*n))
            self.assertEqual(y.haar_moment(2*n-1), 0)

    def test_single_link_quaternion_with_unit_spectator(self):
        q = sp.symbols('q0:4', real=True)
        spectator = (sp.Rational(3, 5), sp.Rational(4, 5), 0, 0)
        f = y.quat(q, spectator)[0]
        self.assertEqual(y.sphere_laplacian(f, q), -3*f)
        self.assertEqual(sp.expand(sum(sp.diff(f, x)**2 for x in q)), 1)

    def test_free_energy_closed_form_via_outward_geometric_tail(self):
        # Positive series tail <= t^(N+1)/(1-t); independent of closed summation.
        for s in (sp.Rational(3, 5), sp.Rational(4, 5), sp.Rational(12, 13)):
            t, n = 1-s*s, 24
            lower = sum((y.haar_moment(2*j)-y.haar_moment(2*j+2))*t**(j+1) for j in range(n))
            upper = lower+t**(n+1)/(1-t)
            energy = y.free_energy_from_s(s)
            self.assertLessEqual(lower, energy)
            self.assertLessEqual(energy, upper)

    def test_escape_normalization_not_outcome_entropy(self):
        # Equal probabilities at eta=0 have entropy log2 and zero excitation.
        omega = sp.Matrix([1, 0])
        m = sp.eye(2)/sp.sqrt(2)
        post = 2*m*omega*omega.T*m
        self.assertEqual(post, omega*omega.T)
        self.assertEqual((omega.T*m*m*omega)[0], sp.Rational(1, 2))
        self.assertEqual(1-(omega.T*post*omega)[0], 0)

    def test_one_probe_is_only_an_upper_gap_bound(self):
        h = sp.diag(0, sp.Rational(1, 100), 12)
        k = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        e, q, *_ = y.sharp_data(h, sp.Matrix([1, 0, 0]), k)
        self.assertEqual(e/q, 12)
        self.assertLess(h[1, 1], e/q)

    def test_source_pins_and_no_continuum_overclaim(self):
        record = json.loads((y.HERE/'YC23_RESULT.json').read_text())
        for path, digest in record['source_sha256'].items():
            self.assertEqual(hashlib.sha256((y.ROOT/path).read_bytes()).hexdigest(), digest)
        self.assertTrue(all(record['checks'].values()))
        self.assertFalse(record['claims']['continuum_mass_gap'])
        self.assertFalse(record['claims']['primitive_lambda_space_uniquely_derived'])


if __name__ == '__main__':
    unittest.main()
