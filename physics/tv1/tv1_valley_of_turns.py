"""TV1 -- the valley of d turns: the quadratic layer closes for three turns and is logarithmic for four.

Record: d unit blocks U_i = (cos a_i ; u_i), |u_i| = s_i = sin a_i, with the commutator weight of every pair,
1 - 2 X_ij, X_ij = |u_i x u_j|^2 (HL1-H7).  The valley is X = 0: all vector parts on one axis n.

V1  exact split about any axis n, u_i = p_i n + t_i (t_i across n):
        X_ij = |p_i t_j - p_j t_i|^2 + |t_i x t_j|^2 :   a quadratic layer and a quartic core.
    At a centre point (all p = 0) the quadratic layer is zero.
V2  on the valley (p_i = s_i, t_i = s_i delta_i) the quadratic layer is the graph form
        sum over pairs of s_i^2 s_j^2 |delta_i - delta_j|^2 ;
    its reduced determinant is (prod s_i^2)(sum s_i^2)^(d-2)  (one common shift of the axis is free).
V3  so the weight the quadratic layer leaves on the valley is  (sum sin^2 a_i)^-(d-2)  (the blocks' own measure
    cancels the product).  Its mean over the angles is a sum over closed paths on the lattice of d directions:
        d = 2 :  1
        d = 3 :  (2/3) sum c_n(3)/6^n                 finite
        d = 4 :  (1/4) sum (n + 1) c_n(4)/8^n         grows without bound, by Log 2/pi^2 per doubling
    c_n(d) = closed paths of n steps = (n)! [x^n] (RW1's returned count at nu = 0)^d.

Exact integers and rationals.  Stdlib only.
"""
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))

M_TOP = 512                       # closed paths up to 2 M_TOP steps


# ---------------------------------------------------------------- vectors
def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c*x for x in a)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def across(w, n):
    """the part of w across the unit axis n"""
    return sub(w, scale(dot(w, n), n))


AXES = [(F(2, 3), F(2, 3), F(1, 3)), (F(3, 5), F(0), F(4, 5)), (F(6, 7), F(2, 7), F(3, 7)), (F(0), F(0), F(1))]


# ---------------------------------------------------------------- polynomials in x_1 .. x_d (dict exponent -> integer)
def pmul(p, q):
    out = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            out[e] = out.get(e, 0) + c1*c2
    return {e: c for e, c in out.items() if c}


def padd(p, q, c=1):
    out = dict(p)
    for e, x in q.items():
        out[e] = out.get(e, 0) + c*x
    return {e: x for e, x in out.items() if x}


def var(i, d):
    e = [0]*d
    e[i] = 1
    return {tuple(e): 1}


def perm_sign(p):
    sign, seen = 1, [False]*len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, length = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                length += 1
            if length % 2 == 0:
                sign = -sign
    return sign


def pdet(mat, d):
    n = len(mat)
    total = {}
    for p in permutations(range(n)):
        term = {tuple([0]*d): perm_sign(p)}
        for i in range(n):
            term = pmul(term, mat[i][p[i]])
            if not term:
                break
        total = padd(total, term)
    return total


def graph_form(d):
    """the reduced matrix of sum_(i<j) x_i x_j (y_i - y_j)^2 with y_d = 0, entries polynomials in x"""
    xs = [var(i, d) for i in range(d)]
    mat = [[{} for _ in range(d - 1)] for _ in range(d - 1)]
    for i in range(d - 1):
        for j in range(d - 1):
            if i == j:
                for k in range(d):
                    if k != i:
                        mat[i][j] = padd(mat[i][j], pmul(xs[i], xs[k]))
            else:
                mat[i][j] = padd({}, pmul(xs[i], xs[j]), -1)
    return mat


def tree_side(d):
    """(prod x_i)(sum x_i)^(d-2)"""
    xs = [var(i, d) for i in range(d)]
    total = {}
    for x in xs:
        total = padd(total, x)
    out = {tuple([1]*d): 1}
    for _ in range(d - 2):
        out = pmul(out, total)
    return out


# ---------------------------------------------------------------- closed paths
_WEIGHTS = {}


def closed_paths(d, top=M_TOP):
    """c_(2m)(d), m = 0 .. top: closed paths of 2m steps on the lattice of d directions.
    One more direction: choose which 2k of the 2m steps it takes, and its own closed path (C(2k, k) ways)."""
    if top not in _WEIGHTS:
        _WEIGHTS[top] = [[comb(2*m, 2*k)*comb(2*k, k) for k in range(m + 1)] for m in range(top + 1)]
    w = _WEIGHTS[top]
    key = ('paths', top)
    table = _WEIGHTS.setdefault(key, {0: [1] + [0]*top})
    have = max(table)
    while have < d:
        prev = table[have]
        table[have + 1] = [sum(w[m][k]*prev[m - k] for k in range(m + 1)) for m in range(top + 1)]
        have += 1
    return table[d]


def closed_paths_by_walking(d, steps):
    moves = [tuple((s if i == k else 0) for i in range(d)) for k in range(d) for s in (1, -1)]
    count = 0
    for path in product(moves, repeat=steps):
        if all(sum(m[i] for m in path) == 0 for i in range(d)):
            count += 1
    return count


def returned_count_power(d, top):
    """(2m)! [x^(2m)] (sum_j x^(2j)/(j!)^2)^d : d returned counts of two strands (RW1, nu = 0) sharing the steps"""
    base = [F(1, factorial(j)**2) for j in range(top + 1)]
    cur = [F(1)] + [F(0)]*top
    for _ in range(d):
        nxt = [F(0)]*(top + 1)
        for i, x in enumerate(cur):
            for j in range(top + 1 - i):
                nxt[i + j] += x*base[j]
        cur = nxt
    return [int(cur[m]*factorial(2*m)) for m in range(top + 1)]


def run():
    out, num = {}, {}
    rnd = random.Random(20261009)

    def rv():
        return tuple(F(rnd.randint(-9, 9), rnd.randint(1, 5)) for _ in range(3))

    # ---- V1: the exact split
    ok, core_only, layer_only = True, True, True
    for trial in range(200):
        n = AXES[trial % 4]
        u, v = rv(), rv()
        pu, pv, tu, tv = dot(u, n), dot(v, n), across(u, n), across(v, n)
        x = dot(cross(u, v), cross(u, v))
        layer = sub(scale(pu, tv), scale(pv, tu))
        core = cross(tu, tv)
        ok &= dot(n, n) == 1 and dot(tu, n) == 0 and x == dot(layer, layer) + dot(core, core)
        # all p = 0 (the vectors t_u, t_v themselves): the layer is zero, only the core is left
        zero = sub(scale(dot(tu, n), tv), scale(dot(tv, n), tu))
        core_only &= zero == (0, 0, 0) and dot(cross(tu, tv), cross(tu, tv)) == dot(core, core) and \
            (dot(core, core) == 0) == (cross(tu, tv) == (0, 0, 0))
        # one vector on the axis: only the layer is left
        w = scale(pu, n)
        layer_only &= dot(cross(w, v), cross(w, v)) == pu*pu*dot(tv, tv)
    out['V1 |u x v|^2 = |p_u t_v - p_v t_u|^2 + |t_u x t_v|^2 about any axis (200 exact pairs, four axes)'] = ok
    out['V1 with no part along the axis only the quartic core is left; with one vector on the axis only the layer'] = \
        core_only and layer_only
    ok = True
    for trial in range(100):
        n = AXES[trial % 4]
        s = [F(rnd.randint(1, 9), 10) for _ in range(3)]
        dl = [across(rv(), n) for _ in range(3)]
        us = [scale(s[i], add(n, dl[i])) for i in range(3)]
        for i, j in ((0, 1), (1, 2), (2, 0)):
            x = dot(cross(us[i], us[j]), cross(us[i], us[j]))
            diff = sub(dl[i], dl[j])
            ok &= x == s[i]**2*s[j]**2*(dot(diff, diff) + dot(cross(dl[i], dl[j]), cross(dl[i], dl[j])))
    out['V1 near the valley, u_i = s_i (n + delta_i): X_ij = s_i^2 s_j^2 (|delta_i - delta_j|^2 + |delta_i x delta_j|^2) exactly'] = ok

    # ---- V2: the graph form and its determinant
    ok, trees = True, []
    for d in range(2, 7):
        det = pdet(graph_form(d), d)
        ok &= det == tree_side(d)
        trees.append(sum(det.values()))
    out['V2 reduced determinant of sum x_i x_j (y_i - y_j)^2 is (prod x_i)(sum x_i)^(d-2), d = 2..6 (polynomial identity)'] = ok
    out['V2 at x = 1 it counts the trees on d points: 1, 3, 16, 125, 1296'] = trees == [1, 3, 16, 125, 1296]
    ok = True
    for d in (3, 4):
        for _ in range(50):
            x = [F(rnd.randint(1, 9), rnd.randint(1, 9)) for _ in range(d)]
            # the full form kills the common shift and nothing else
            full = [[(x[i]*(sum(x) - x[i]) if i == j else -x[i]*x[j]) for j in range(d)] for i in range(d)]
            ok &= all(sum(row) == 0 for row in full)
            red = sum(c*_prod(x, e) for e, c in pdet(graph_form(d), d).items())
            prod_x = F(1)
            for t in x:
                prod_x *= t
            ok &= red == prod_x*sum(x)**(d - 2) and red > 0
        # one turn at a centre point: the determinant is zero
        x = [F(0)] + [F(1)]*(d - 1)
        ok &= sum(c*_prod(x, e) for e, c in pdet(graph_form(d), d).items()) == 0
    out['V2 the common shift of the axis is the only free direction when no s_i is zero; a turn at a centre point frees more'] = ok
    out['V2 weight left on the valley: (prod s_i^2) from the blocks over (prod s_i^2)(sum s_i^2)^(d-2) = (sum s_i^2)^-(d-2)'] = \
        all(tree_side(d) == pmul({tuple([1]*d): 1}, _power_sum(d, d - 2)) for d in range(2, 7))

    # ---- V3: the mean of the valley weight as a sum over closed paths
    c3, c4, c2 = closed_paths(3), closed_paths(4), closed_paths(2, 40)
    out['V3 closed paths: 6, 90, 1860 in three directions; 8, 168, 5120 in four (three routes for the first)'] = \
        c3[1:4] == [6, 90, 1860] and c4[1:4] == [8, 168, 5120] and \
        [closed_paths_by_walking(3, k) for k in (2, 4, 6)] == [6, 90, 1860] and \
        [closed_paths_by_walking(4, k) for k in (2, 4)] == [8, 168] and \
        returned_count_power(3, 12) == c3[:13] and returned_count_power(4, 12) == c4[:13] and \
        c2[:6] == [comb(2*m, m)**2 for m in range(6)]
    # (sum sin^2 a_i) = (d/2)(1 - y), y = (1/d) sum cos 2a_i ; mean of y^n = c_n/(2d)^n ;
    # (1 - y)^-(d-2) = sum C(n + d - 3, d - 3) y^n
    sums3, sums4, s3, s4 = {}, {}, F(0), F(0)
    for m in range(M_TOP + 1):
        s3 += F(2, 3)*F(c3[m], 6**(2*m))
        s4 += F(1, 4)*(2*m + 1)*F(c4[m], 8**(2*m))
        if m in (8, 16, 32, 64, 128, 256, 512):
            sums3[m], sums4[m] = s3, s4
    marks = (8, 16, 32, 64, 128, 256, 512)
    out['V3 three turns: the partial sums rise and stay below 11/4 (the mean is finite)'] = \
        all(sums3[x] < sums3[y] for x, y in zip(marks, marks[1:])) and sums3[512] < F(11, 4) and sums3[512] > F(99, 100)
    steps3 = [sums3[y] - sums3[x] for x, y in zip(marks, marks[1:])]
    out['V3 three turns: each doubling adds less than the one before (by more than a fifth)'] = \
        all(5*b < 4*a for a, b in zip(steps3, steps3[1:]))
    steps4 = [sums4[y] - sums4[x] for x, y in zip(marks, marks[1:])]
    out['V3 four turns: each doubling adds more than the one before, between 1/15 and 71/1000 from 2 x 16 steps on'] = \
        all(b > a for a, b in zip(steps4, steps4[1:])) and all(F(1, 15) < t < F(71, 1000) for t in steps4[1:])
    num['three turns: partial sums at 2m = 16 .. 1024'] = [round(float(sums3[x]), 6) for x in marks]
    num['four turns: partial sums at 2m = 16 .. 1024'] = [round(float(sums4[x]), 6) for x in marks]
    num['four turns: added per doubling'] = [round(float(t), 6) for t in steps4]
    c5 = closed_paths(5)
    s5, sums5 = F(0), {}
    for m in range(M_TOP + 1):
        s5 += comb(2*m + 2, 2)*F(c5[m], 10**(2*m))                          # C(n + 2, 2) c_n / 10^n
        if m in marks:
            sums5[m] = s5
    steps5 = [sums5[y] - sums5[x] for x, y in zip(marks, marks[1:])]
    out['V3 two turns: the weight is 1 (the power is zero); five turns: each doubling adds about sqrt 2 times the one before'] = \
        tree_side(2) == {(1, 1): 1} and all(F(13, 10)*a < b < F(3, 2)*a for a, b in zip(steps5, steps5[1:]))
    # the walk inequality's exponents: near a centre point the weight is |a|^-2(d-2) against the volume |a|^(d-1) d|a|
    out['V3 the count of powers: radial exponent 3 - d ; integrable at a centre point iff d < 4 ; d = 4 is the logarithm'] = \
        [(d - 1) - 2*(d - 2) for d in (2, 3, 4, 5)] == [1, 0, -1, -2]
    return out, num


def _prod(x, e):
    out = F(1)
    for t, k in zip(x, e):
        out *= t**k
    return out


def _power_sum(d, k):
    total = {}
    for i in range(d):
        total = padd(total, var(i, d))
    out = {tuple([0]*d): 1}
    for _ in range(k):
        out = pmul(out, total)
    return out


if __name__ == '__main__':
    out, num = run()
    for k, v in out.items():
        print('PASS' if v else 'FAIL', k)
    for k, v in num.items():
        print(k, v)
    with open(os.path.join(HERE, 'TV1_RESULT.json'), 'w') as fh:
        json.dump({'checks': {k: bool(v) for k, v in out.items()}, 'numbers': num}, fh, indent=1)
        fh.write('\n')
