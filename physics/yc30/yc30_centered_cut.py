#!/usr/bin/env python3
"""Focused exact checks for YC30's new normalization and cut bounds.

No harmonic cutoff, predecessor regression, or new gap computation.
"""

from fractions import Fraction as F
from itertools import combinations

import sympy as s


def equal(a, b):
    assert all(s.simplify(x) == 0 for x in a - b)


def psd(a):
    equal(a, a.T)
    for size in range(1, a.rows + 1):
        for idx in combinations(range(a.rows), size):
            assert a.extract(idx, idx).det() >= 0


def cut(a, retained):
    p = tuple(retained)
    q = tuple(i for i in range(a.rows) if i not in p)
    ap = a.extract(p, p)
    aq = a.extract(q, q)
    cross = a.extract(q, p)
    hidden = -aq.inv() * cross
    w = s.zeros(a.rows, len(p))
    for column, row in enumerate(p):
        w[row, column] = 1
    for i, row in enumerate(q):
        for j in range(len(p)):
            w[row, j] = hidden[i, j]
    reduced = ap - a.extract(p, q) * aq.inv() * cross
    equal(w.T * a * w, reduced)
    return reduced, w, q


def sign_support():
    # A rectangular source has an unpaired retained direction.
    u = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5)]])
    x = s.zeros(3)
    x[:2, 2:3] = u.T
    x[2:3, :2] = u
    j = s.diag(1, 1, -1)
    support = s.diag(1, 1, 1)
    support[:2, :2] = u.T * u
    equal(x*x, support)
    equal(j*x, -x*j)
    equal((j*x)**2, -support)
    assert support != s.eye(3)
    print("Sign exchange: involution on source support; zero mode retained")


def exact_cuts():
    a = s.Matrix([[5, 1, 0, 1], [1, 6, 1, 0],
                  [0, 1, 7, 1], [1, 0, 1, 8]])
    tau, m = s.Integer(4), s.Integer(3)
    psd(a - m*s.eye(4))
    k = tau * (a + tau*s.eye(4)).inv()
    f = s.eye(4) - k
    q = tau/(m+tau)
    delta = 1-q
    psd(k)
    psd(q*s.eye(4)-k)

    sa, wa, _ = cut(a, (0, 1))
    sf, wf, _ = cut(f, (0, 1))
    ke = s.eye(2)-sf
    equal(sf, sa * (sa+tau*s.eye(2)).inv())
    equal(wf, (a+tau*s.eye(4))*wa*(sa+tau*s.eye(2)).inv())
    psd(ke)
    psd(q*s.eye(2)-ke)
    psd(s.eye(2)/delta-wf.T*wf)
    psd(sa/m-wa.T*wa)

    first, w1, _ = cut(f, (0, 1, 2))
    final, w2, _ = cut(first, (0,))
    direct, wd, hidden = cut(f, (0,))
    equal(final, direct)
    equal(w1*w2, wd)
    source = s.Matrix([1, -2, 3, -1])
    equal(w2.T*w1.T*source, wd.T*source)
    particular = s.zeros(4, 1)
    hidden_solution = f.extract(hidden, hidden).inv()*source.extract(hidden, [0])
    for i, row in enumerate(hidden):
        particular[row, 0] = hidden_solution[i, 0]
    equal(wd*direct.inv()*wd.T*source+particular, f.inv()*source)
    print("Exact energy normalization, nested lifts, full sources and both metric bounds")

    epsilon = s.Rational(1, 64)
    v = s.Matrix([0, 0, 0, 1])
    kt = k + epsilon*v*v.T
    ft = s.eye(4)-kt
    psd(kt)
    psd(q*s.eye(4)-kt)
    st, wt, _ = cut(ft, (0, 1))
    equal(sf-st, wt.T*(f-ft)*wf)
    hidden_inverse = s.zeros(4)
    hidden_inverse[2:4, 2:4] = f[2:4, 2:4].inv()
    equal(wf-wt, -hidden_inverse*(f-ft)*wt)
    diff = sf-st
    psd((epsilon/delta)**2*s.eye(2)-diff*diff)
    metric_diff = wf.T*wf-wt.T*wt
    psd((2*epsilon/delta**2)**2*s.eye(2)-metric_diff*metric_diff)
    print("Full-input return error and bounded-lift metric error: exact identities and bounds")


def derivative_metric():
    z = s.symbols('z', real=True)
    a = s.Matrix([[3, 1], [1, 5]])-z*s.eye(2)
    tau = s.Integer(4)
    f = a*(a+tau*s.eye(2)).inv()
    sa, wa, _ = cut(a, (0,))
    sf, wf, _ = cut(f, (0,))
    ma = wa.T*wa
    equal(-sa.diff(z), ma)
    equal(-sf.diff(z), tau*(sa+tau*s.eye(1)).inv()*ma*(sa+tau*s.eye(1)).inv())
    equal(-sf.diff(z), wf.T*(tau*(a+tau*s.eye(2)).inv()**2)*wf)
    print("Nonlinear normalized pencil: derivative metric equals the converted physical metric")


def ym_arithmetic():
    assert F(1, 4*12*24) == F(1, 1152)
    assert F(1, 4*13*25) == F(1, 1300)
    for name, gap, expected in (
        ('28-link', F(592, 175), (F(175, 323), F(148, 323), F(323, 148), F(175, 148))),
        ('64-link', F(11188, 3325), (F(3325, 6122), F(2797, 6122), F(6122, 2797), F(3325, 2797))),
    ):
        q = 4/(gap+4)
        reserve = 1-q
        actual = (q, reserve, 1/reserve, 4/gap)
        assert actual == expected
        print(f'{name}: q={q}, reserve={reserve}, bounded metric<={1/reserve}, physical E<=4 metric<={4/gap}')
    print('Raw reference return lower bounds: lambda^2 M2/1300; free point /1152')


if __name__ == '__main__':
    sign_support()
    exact_cuts()
    derivative_metric()
    ym_arithmetic()
