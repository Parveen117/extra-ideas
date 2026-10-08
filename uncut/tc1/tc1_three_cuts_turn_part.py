"""TC1: with three cuts the part that a self-dagger first cut leaves out is a second four-reading.

Carrier and cuts as in FR1/IN1: C1 = K, C2 = RK, C3 = iota R, iota central, iota^2 = -1.
An element is rho = h + iota k with h = n + r.C and k = d + e.C both self-dagger.  A change of frame is rho -> G rho G-dagger, det G = 1.
Sources read (unchanged): extra-ideas physics fr1 (three cuts, no fourth), in1 (reading tensor, det = n^2 - r.r),
  em1 (F = E + iota B; F.F = (E.E - B.B) + 2 iota E.B), dm1 (the third cut is the unseen leg), gb1 B1 (a reading is blind to the phase turn),
  ob1 (content q counted along the phase); uncut fc1 (two cuts: m^2 and W), up7.
No measurement.  sympy, exact.  Python 3.12.
"""
import json

import sympy as sp

i = sp.I
I2 = sp.eye(2)
K = sp.Matrix([[1, 0], [0, -1]])
R = sp.Matrix([[0, -1], [1, 0]])
C = [K, R*K, i*R]


def dag(M):
    return M.H


def vec(n, r):
    return n*I2 + sum((r[j]*C[j] for j in range(3)), sp.zeros(2, 2))


def split(M):
    """rho = h + iota k with h, k self-dagger; returns ((n, r), (d, e))"""
    h, k = (M + dag(M))/2, (M - dag(M))/(2*i)
    def comps(X):
        return sp.simplify(X.trace()/2), [sp.simplify((X*C[j]).trace()/2) for j in range(3)]
    return comps(h), comps(k)


def dot(a, b):
    return a[0]*b[0] - sum(x*y for x, y in zip(a[1], b[1]))


def zero(e):
    return sp.simplify(sp.expand(sp.sympify(e).rewrite(sp.exp))) == 0


def boost(eta, j):
    return sp.cosh(eta/2)*I2 + sp.sinh(eta/2)*C[j]


def turn(theta, j):
    return sp.cos(theta/2)*I2 - i*sp.sin(theta/2)*C[j]


def run():
    res = {}
    for a_, b_ in ((0, 1), (0, 2), (1, 2)):
        if C[a_]*C[b_] + C[b_]*C[a_] != sp.zeros(2, 2):
            raise ValueError('cuts')
    if C[0]*C[1]*C[2] != -i*I2 and C[0]*C[1]*C[2] != i*I2:
        raise ValueError('iota is not the product of the three cuts')
    if R != -i*C[2]:
        raise ValueError('the two-cut turn')
    res['turn_of_two_cuts'] = 'R = -iota C3 : the turn of the two-cut carrier is iota times the third cut'

    n, d = sp.symbols('n d', real=True)
    r = sp.symbols('r1:4', real=True)
    e = sp.symbols('e1:4', real=True)
    rho = vec(n, r) + i*vec(d, e)
    h, k = (n, list(r)), (d, list(e))

    # T1 the split is unique and each half is a four-reading under every change of frame
    (hn, hr), (kd, ke) = split(rho)
    if not (zero(hn - n) and zero(kd - d) and all(zero(x - y) for x, y in zip(hr, r)) and all(zero(x - y) for x, y in zip(ke, e))):
        raise ValueError('split')
    eta, th = sp.Rational(0), sp.Rational(0)
    frames = [sp.Matrix([[2, 0], [0, sp.Rational(1, 2)]]),                              # boost along C1, exp(eta) = 4
              sp.Matrix([[sp.Rational(5, 4), sp.Rational(3, 4)], [sp.Rational(3, 4), sp.Rational(5, 4)]]),   # boost along C2
              sp.Matrix([[sp.Rational(3, 5), -sp.Rational(4, 5)], [sp.Rational(4, 5), sp.Rational(3, 5)]]),  # turn about C3
              sp.Matrix([[sp.Rational(3, 5) - sp.Rational(4, 5)*i, 0], [0, sp.Rational(3, 5) + sp.Rational(4, 5)*i]])]  # turn about C1
    frames.append(frames[1]*frames[3]*frames[0])
    for G in frames:
        if sp.simplify(G.det()) != 1:
            raise ValueError('frame')
        h2, k2 = split(sp.simplify(G*rho*dag(G)))
        if not (zero(dot(h2, h2) - dot(h, h)) and zero(dot(k2, k2) - dot(k, k)) and zero(dot(h2, k2) - dot(h, k))):
            raise ValueError('invariants under a frame change')
        hh, _ = split(sp.simplify(G*vec(n, r)*dag(G)))
        if not all(zero(x - y) for x, y in zip([h2[0]] + h2[1], [hh[0]] + hh[1])):
            raise ValueError('halves do not transform separately')
    res['invariants'] = 'three: h.h , k.k , h.k   (two cuts had two: m^2 and W)'

    # T2 the determinant: det rho = (h.h - k.k) + 2 iota h.k   - the same form as F.F for the light field
    if not zero(rho.det() - (dot(h, h) - dot(k, k)) - 2*i*dot(h, k)):
        raise ValueError('determinant')
    res['determinant'] = 'det rho = (h.h - k.k) + 2 iota h.k'

    # T3 the two-cut case inside it: only C1, C2 and the turn R = -iota C3
    W = sp.Symbol('W', real=True)
    two = n*I2 + r[0]*C[0] + r[1]*C[1] + W*R
    (hn2, hr2), (kd2, ke2) = split(two)
    if not (zero(kd2) and zero(ke2[0]) and zero(ke2[1]) and zero(ke2[2] + W) and zero(hr2[2])):
        raise ValueError('two-cut embedding')
    if not zero(two.det() - (n*n - r[0]**2 - r[1]**2) - W*W):
        raise ValueError('two-cut determinant')
    res['two_cut_case'] = 'W = minus the third component of k : the missed part points along the cut the two-cut frame does not have'

    # T4 the time part of k is a phase and no reading sees it; multiplying by a phase turns h and k into each other
    chi = sp.Symbol('chi', real=True)
    Ph = sp.exp(i*chi/2)*I2
    if sp.simplify(Ph*rho*dag(Ph) - rho) != sp.zeros(2, 2):
        raise ValueError('phase as a frame change')
    h3, k3 = split(sp.simplify((sp.cos(chi) + i*sp.sin(chi))*rho))
    if not (zero(h3[0] - (sp.cos(chi)*n - sp.sin(chi)*d)) and zero(k3[0] - (sp.sin(chi)*n + sp.cos(chi)*d))):
        raise ValueError('phase multiplication')
    res['phase'] = 'a frame change by a phase does nothing ; multiplying the element by Exp(iota chi) turns (h, k) by chi'

    # T5 at rest: h = (m; 0).  h.k = m d.  With d = 0 the determinant is real: M^2 = m^2 + |e|^2, and e turns as a direction
    m = sp.Symbol('m', positive=True)
    rest = vec(m, (0, 0, 0)) + i*vec(0, e)
    if not zero(rest.det() - m*m - (e[0]**2 + e[1]**2 + e[2]**2)):
        raise ValueError('rest determinant')
    G = frames[2]                                           # turn about C3
    _, k4 = split(sp.simplify(G*rest*dag(G)))
    if not (zero(k4[0]) and zero(k4[1][2] - e[2]) and zero(k4[1][0]**2 + k4[1][1]**2 - e[0]**2 - e[1]**2)):
        raise ValueError('turning e')
    G = frames[0]                                           # boost along C1
    h5, k5 = split(sp.simplify(G*rest*dag(G)))
    if not (zero(dot(h5, k5)) and zero(k5[0] - sp.Rational(15, 8)*e[0]) and zero(k5[1][1] - e[1])):
        raise ValueError('boosting e')
    res['at_rest'] = 'M^2 = m^2 + |e|^2 when h.k = 0 ; e is a direction in space ; moving, k gains a time part with h.k still 0'

    # exact witness
    w = vec(4, (0, 0, 0)) + i*vec(0, (0, 0, -3))
    if w != 4*I2 + 3*R or w.det() != 25:
        raise ValueError('witness')
    res['witness'] = 'FC1: 4 + 3R is h = (4; 0,0,0), k = (0; 0,0,-3), det 25'
    return res


if __name__ == '__main__':
    out = run()
    with open('TC1_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
