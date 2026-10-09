"""TC1 -- the core of d turns: three numbers, a difference form whose diagonal is the valley, a flat measure at
four turns, the constant of the logarithm, and the scale the core makes for itself.

The core: near a centre point the vector parts of d unit blocks are a free 3 x d matrix C (three cut components,
d directions), and the commutator weight of all pairs is  X = sum_(i<j) |c_i x c_j|^2  (HL1-H7, TV1-V1).

C1  X = (1/2)[(tr G)^2 - tr G^2],  G = C C^t (3 x 3):  the second symmetric function e2 of three numbers x_1, x_2, x_3
    (the squared singular values).  Unchanged by turning the three cuts and by turning the d directions.
    A fourth direction adds no fourth number.
C2  the one law:  X = R - D,  R = (1/2)(sum x)^2,  D = (1/2) sum x^2 ;  0 <= X ;  X = 0 exactly on the valley ;
    lost/seen = D/R in [1/3, 1].  As a form in x it has signature (1, 2); the valley is its null line.
C3  about a valley point C = n p^t + T:  X = |p|^2 |T|^2 - |T p|^2 + X(T).  The layer is |p|^2 times the projector
    across p in the space of turns: every transverse rate is the same, and the layer's ground rate is linear in |p|.
C4  the measure in x is  prod |x_a - x_b| prod x_a^((d-4)/2) d^3x :  flat exactly at four turns
    (checked on moments of the free weight, d = 4 and 6).
C5  four turns, flat core, largest number cut at Lambda:  K(2 Lambda) - K(Lambda) = (Log 2)/16 - delta,
    0 <= delta <= delta(Lambda) -> 0.  With a coupling kappa the record is kappa^-3 K(Lambda sqrt(kappa)).
C6  X has degree four in C: the Hamiltonian -(g^2/2) Laplacian + X/g^2 is g^(2/3) times the one at g = 1.
    Gaussian reading: ground rate <= c_d g^(2/3), c_d^3 = (729/256) d^3 (d - 1).

Exact rationals.  Stdlib only.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, factorial
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- small exact linear algebra
def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])


def mat_mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(r) for r in zip(*a)]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def det(a):
    a = [row[:] for row in a]
    n, sign, out = len(a), 1, F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            sign = -sign
        out *= a[c][c]
        for r in range(c + 1, n):
            f = a[r][c]/a[c][c]
            for k in range(c, n):
                a[r][k] -= f*a[c][k]
    return sign*out


def inverse(a):
    n = len(a)
    m = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(a)]
    for c in range(n):
        piv = next(r for r in range(c, n) if m[r][c] != 0)
        m[c], m[piv] = m[piv], m[c]
        m[c] = [x/m[c][c] for x in m[c]]
        for r in range(n):
            if r != c and m[r][c] != 0:
                f = m[r][c]
                m[r] = [x - f*y for x, y in zip(m[r], m[c])]
    return [row[n:] for row in m]


def rotation(n, rnd):
    """an exact rotation of n dimensions: (1 - A)(1 + A)^-1 for a random skew A"""
    a = [[F(0)]*n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = F(rnd.randint(-4, 4), rnd.randint(1, 4))
            a[j][i] = -a[i][j]
    one = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    plus = [[one[i][j] + a[i][j] for j in range(n)] for i in range(n)]
    minus = [[one[i][j] - a[i][j] for j in range(n)] for i in range(n)]
    return mat_mul(minus, inverse(plus))


def weight(c):
    """X = sum over pairs of directions of |c_i x c_j|^2 for a 3 x d matrix"""
    cols = list(zip(*c))
    return sum(dot(cross(cols[i], cols[j]), cross(cols[i], cols[j])) for i, j in combinations(range(len(cols)), 2))


def gram(c):
    return mat_mul(c, transpose(c))


def e2_of(g):
    return (trace(g)**2 - trace(mat_mul(g, g)))/2


# ---------------------------------------------------------------- polynomials (dict exponent tuple -> rational)
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


def odd_fact(k):
    out = 1
    for j in range(1, 2*k, 2):
        out *= j
    return out


def free_mean(poly, var=F(1, 2)):
    """mean under the free weight exp(-sum c^2) (every entry has variance 1/2), by pairing"""
    total = F(0)
    for e, c in poly.items():
        if any(k % 2 for k in e):
            continue
        term = F(c)
        for k in e:
            term *= odd_fact(k//2)*var**(k//2)
        total += term
    return total


def entry_polys(d):
    """e1, e2, e3 of G = C C^t as polynomials in the 3d entries of C"""
    n = 3*d

    def var(a, i):
        e = [0]*n
        e[a*d + i] = 1
        return {tuple(e): 1}
    g = [[{} for _ in range(3)] for _ in range(3)]
    for a in range(3):
        for b in range(3):
            for i in range(d):
                g[a][b] = padd(g[a][b], pmul(var(a, i), var(b, i)))
    e1 = padd(padd(g[0][0], g[1][1]), g[2][2])
    e2 = {}
    for a, b in ((0, 1), (0, 2), (1, 2)):
        e2 = padd(e2, padd(pmul(g[a][a], g[b][b]), pmul(g[a][b], g[b][a]), -1))
    e3 = {}
    for p in permutations(range(3)):
        sign = 1 if p in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
        e3 = padd(e3, pmul(pmul(g[0][p[0]], g[1][p[1]]), g[2][p[2]]), sign)
    return e1, e2, e3


def number_side(d):
    """means of e1, e2, e3, e1^2, e1 e2 under prod|x_a - x_b| prod x^((d-4)/2) exp(-sum x), d even,
    on x_1 > x_2 > x_3 > 0 with x_3 = a, x_2 = a + b, x_1 = a + b + c"""
    a, b, c = {(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}
    x3, x2 = a, padd(a, b)
    x1 = padd(x2, c)
    delta = pmul(pmul(c, padd(b, c)), b)
    e1 = padd(padd(x1, x2), x3)
    e2 = padd(padd(pmul(x1, x2), pmul(x1, x3)), pmul(x2, x3))
    e3 = pmul(pmul(x1, x2), x3)
    w = delta
    for _ in range((d - 4)//2):
        w = pmul(w, e3)

    def integral(p):
        return sum(co*F(factorial(i), 3**(i + 1))*F(factorial(j), 2**(j + 1))*factorial(k)
                   for (i, j, k), co in p.items())
    z = integral(w)
    return [integral(pmul(w, f))/z for f in (e1, e2, e3, pmul(e1, e1), pmul(e1, e2))]


# ---------------------------------------------------------------- C5: the logarithm of four turns
def ray_series():
    """r(x1) = integral over x2 > x3 > 0 of (x1-x2)(x1-x3)(x2-x3)(1 - 2 x2 x3) exp(-2 x1 (x2 + x3)),
    as {power of x1: coefficient}; x3 = c, x2 = c + t; monomial x1^m c^i t^j -> x1^(m-i-j-2) i! j!/(4^(i+1) 2^(j+1))"""
    x1, c, t = {(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}
    one = {(0, 0, 0): 1}
    poly = pmul(pmul(t, padd(padd(x1, c, -1), t, -1)), padd(x1, c, -1))
    poly = pmul(poly, padd(one, pmul(c, padd(c, t)), -2))
    out = {}
    for (m, i, j), co in poly.items():
        k = m - i - j - 2
        out[k] = out.get(k, 0) + co*F(factorial(i)*factorial(j), 4**(i + 1)*2**(j + 1))
    return {k: v for k, v in out.items() if v}


def outside_bound(lam):
    """upper bound of the integral over x1 in [lam, 2 lam] of the part of the comparison integral with x2 > x1:
    |integrand| <= x2^3 (1 + 2 x2^2) exp(-s x2) exp(-s x3), s = 2 x1"""
    def tail(n, x):          # integral over [x, inf) of y^n exp(-2 x y) dy, without the factor exp(-2 x^2)
        return sum(F(factorial(n), factorial(k))*x**k/(2*x)**(n - k + 1) for k in range(n + 1))
    # times exp(-2 x^2) every term of the bracket is x^k exp(-2 x^2) with k <= 3, which falls for x >= 1:
    # the value at lam bounds the product on [lam, 2 lam] (lam >= 2 is used)
    assert lam >= 2
    bracket = (tail(3, lam) + 2*tail(5, lam))/(2*lam)
    t = 2*lam*lam
    exp_bound = 1/sum(t**k/F(factorial(k)) for k in range(12))          # exp(-t) <= 1/(partial sum of exp(t))
    return lam*exp_bound*bracket


def delta_of(lam, series):
    """delta(Lambda): (Log 2)/16 - delta(Lambda) <= K(2 Lambda) - K(Lambda) <= (Log 2)/16"""
    total = F(0)
    for k, v in series.items():
        if k == -1:
            continue
        p = k + 1                                           # integral of x^k over [lam, 2 lam], k != -1
        total += v*((2*lam)**p - lam**p)/p
    return -total + outside_bound(lam)


def run():
    out, num = {}, {}
    rnd = random.Random(20261009)

    def rmat(rows, cols):
        return [[F(rnd.randint(-5, 5), rnd.randint(1, 4)) for _ in range(cols)] for _ in range(rows)]

    # ---- C1: three numbers
    ok, ok_rot, ok_rank = True, True, True
    for d in (2, 3, 4, 5, 6):
        for _ in range(20):
            c = rmat(3, d)
            g = gram(c)
            minors = sum(g[a][a]*g[b][b] - g[a][b]*g[b][a] for a, b in ((0, 1), (0, 2), (1, 2)))
            ok &= weight(c) == e2_of(g) == minors
            small = mat_mul(transpose(c), c)                                    # d x d
            ok &= e2_of(small) == e2_of(g)
            ok_rot &= weight(mat_mul(rotation(3, rnd), c)) == weight(c) == weight(mat_mul(c, rotation(d, rnd)))
            if d >= 4:
                ok_rank &= all(det([[small[i][j] for j in idx] for i in idx]) == 0
                               for idx in combinations(range(d), 4))
    out['C1 X = (1/2)[(tr G)^2 - tr G^2] = sum of the 2 x 2 minors of G = C C^t, d = 2..6 (100 exact matrices)'] = ok
    out['C1 X is unchanged by turning the three cuts and by turning the d directions (exact rotations)'] = ok_rot
    out['C1 a fourth direction adds no fourth number: every 4 x 4 minor of C^t C is zero (d = 4, 5, 6)'] = ok_rank

    # ---- C2: the one law
    ok = True
    for _ in range(100):
        c = rmat(3, rnd.choice((3, 4)))
        g = gram(c)
        r, dd = trace(g)**2/2, trace(mat_mul(g, g))/2
        ok &= weight(c) == r - dd and r - dd >= 0 and 3*dd >= r >= dd
    n, p = [F(2, 3), F(2, 3), F(1, 3)], [F(3), F(-1), F(2), F(5)]
    valley = [[n[a]*p[i] for i in range(4)] for a in range(3)]
    spread = [[F(int(a == i)) for i in range(4)] for a in range(3)]
    gs = gram(spread)
    out['C2 X = R - D with R = (tr G)^2/2, D = tr G^2/2 ; X >= 0 ; 1/3 <= D/R <= 1'] = ok
    out['C2 X = 0 on the valley (D = R, the diagonal) ; D/R = 1/3 when the three numbers are equal'] = \
        weight(valley) == 0 and trace(gram(valley))**2 == trace(mat_mul(gram(valley), gram(valley))) and \
        3*trace(mat_mul(gs, gs)) == trace(gs)**2 and weight(spread) == 3
    form = [[F(0), F(1), F(1)], [F(1), F(0), F(1)], [F(1), F(1), F(0)]]               # 2 e2(x) = x^t form x
    ok = True
    for lam, mult in ((F(2), 1), (F(-1), 2)):
        shifted = [[form[i][j] - (lam if i == j else 0) for j in range(3)] for i in range(3)]
        ok &= det(shifted) == 0
    out['C2 as a form in the three numbers e2 has values 2, -1, -1: signature (1, 2), null on the valley'] = \
        ok and trace(form) == 0 and det(form) == 2 and \
        all(sum(form[i][j]*v[i]*v[j] for i in range(3) for j in range(3)) == 0 for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))

    # ---- C3: the layer about a valley point
    ok, ok_proj = True, True
    axes = [(F(2, 3), F(2, 3), F(1, 3)), (F(3, 5), F(0), F(4, 5)), (F(6, 7), F(2, 7), F(3, 7))]
    for trial in range(60):
        d = rnd.choice((3, 4, 5))
        n = axes[trial % 3]
        p = [F(rnd.randint(-5, 5), rnd.randint(1, 3)) for _ in range(d)]
        raw = rmat(3, d)
        cols = []
        for i in range(d):
            w = tuple(raw[a][i] for a in range(3))
            k = dot(w, n)
            cols.append(tuple(w[a] - k*n[a] for a in range(3)))                 # across n
        t = [[cols[i][a] for i in range(d)] for a in range(3)]
        c = [[n[a]*p[i] + t[a][i] for i in range(d)] for a in range(3)]
        tp = [sum(t[a][i]*p[i] for i in range(d)) for a in range(3)]
        frob = sum(x*x for row in t for x in row)
        ok &= weight(c) == dot(p, p)*frob - dot(tp, tp) + weight(t)
        # the layer in the space of turns: q(v) = |p|^2 |v|^2 - (p.v)^2
        v = [F(rnd.randint(-5, 5)) for _ in range(d)]
        along = [x for x in p]
        across_v = [v[i] - (dot(v, p)/dot(p, p))*p[i] for i in range(d)] if dot(p, p) else v
        q = lambda y: dot(p, p)*dot(y, y) - dot(p, y)**2
        ok_proj &= q(along) == 0 and (dot(p, p) == 0 or q(across_v) == dot(p, p)*dot(across_v, across_v))
    out['C3 about a valley point C = n p^t + T: X = |p|^2 |T|^2 - |T p|^2 + X(T) (60 exact cases, d = 3, 4, 5)'] = ok
    out['C3 the layer is |p|^2 times the projector across p in the space of turns: one rate for all transverse directions'] = ok_proj
    ok = True
    for _ in range(50):
        a, t = F(rnd.randint(1, 9), rnd.randint(1, 9)), F(rnd.randint(-9, 9), rnd.randint(1, 9))
        # psi = exp(-a t^2): psi''/psi = 4 a^2 t^2 - 2a ; with |p|^2 = 2 a^2 the layer -(1/2) d^2 + |p|^2 t^2 gives rate a
        ok &= -F(1, 2)*(4*a*a*t*t - 2*a) + 2*a*a*t*t == a
    out['C3 ground rate of one transverse direction is |p|/sqrt 2, so the layer rises along the valley as sqrt 2 (d - 1)|p|'] = ok
    ok = True
    for d in (3, 4):
        for _ in range(20):
            p = [F(rnd.randint(1, 6), rnd.randint(1, 3)) for _ in range(d)]
            pp = dot(p, p)
            m = [[p[i]*((pp if i == j else 0) - p[i]*p[j])*p[j] for j in range(d - 1)] for i in range(d - 1)]
            prod = F(1)
            for x in p:
                prod *= x*x
            ok &= det(m) == prod*pp**(d - 2)
    out['C3 in the variables of TV1 the layer has the tree determinant (prod p_i^2)(sum p_i^2)^(d-2)'] = ok

    # ---- C4: the measure in the three numbers
    rows = {}
    for d in (4, 6):
        e1, e2, e3 = entry_polys(d)
        entries = [free_mean(f) for f in (e1, e2, e3, pmul(e1, e1), pmul(e1, e2))]
        rows[d] = entries
        out['C4 d = %d: means of e1, e2, e3, e1^2, e1 e2 from the entries equal those from '
            'prod|x_a - x_b| prod x^(%d) exp(-sum x)' % (d, (d - 4)//2)] = \
            entries == number_side(d) and entries[:3] == [F(3*d, 2), F(3, 2)*comb(d, 2), F(3, 4)*comb(d, 3)]
    out['C4 control: the flat measure of four turns does not give the means of six'] = rows[6] != number_side(4)
    num['four turns: means of e1, e2, e3, e1^2, e1 e2'] = [str(x) for x in rows[4]]

    # ---- C5: the constant of the logarithm at four turns
    series = ray_series()
    out['C5 along the valley the comparison integral is 1/(16 x1) + lower powers, and the inner integral without '
        'corrections is exactly 1/(16 x1)'] = \
        series[-1] == F(1, 16) and max(series) == -1 and F(1, 2)/F(2)**3 == F(1, 16)
    deltas = {lam: delta_of(F(lam), series) for lam in (2, 4, 8, 16, 32)}
    out['C5 (Log 2)/16 - delta <= K(2 Lambda) - K(Lambda) <= (Log 2)/16 with 0 < delta(Lambda) < 1/(8 Lambda^2), Lambda = 2 .. 32'] = \
        all(0 < deltas[lam] < F(1, 8*lam*lam) for lam in deltas) and \
        all(deltas[2*lam] < deltas[lam] for lam in (2, 4, 8, 16))
    num['delta(Lambda), Lambda = 2 .. 32'] = [float(deltas[lam]) for lam in (2, 4, 8, 16, 32)]
    num['comparison series: powers of x1 and coefficients'] = {str(k): str(v) for k, v in sorted(series.items())}
    x = [F(5), F(3), F(2)]
    cx = [F(7, 3)*t for t in x]
    vand = lambda y: (y[0] - y[1])*(y[0] - y[2])*(y[1] - y[2])
    sym2 = lambda y: y[0]*y[1] + y[0]*y[2] + y[1]*y[2]
    out['C5 scaling: the measure has degree 6 and e2 degree 2 in the numbers, so a coupling kappa gives kappa^-3 K(Lambda sqrt kappa)'] = \
        vand(cx) == F(7, 3)**3*vand(x) and sym2(cx) == F(7, 3)**2*sym2(x)

    # ---- C6: the scale the core makes for itself
    ok, cubes = True, {}
    for d in (3, 4):
        e1, e2, _ = entry_polys(d)
        kin = free_mean(e1)/2                    # <-(1/2) Laplacian> for psi = exp(-|C|^2/2): <|C|^2>/2
        pot = free_mean(e2)
        ok &= kin == F(3*d, 4) and pot == F(3, 2)*comb(d, 2)
        for lam in (F(2), F(3)):                 # psi(C/lam): entries have variance lam^2/2
            ok &= free_mean(e2, lam*lam/2) == lam**4*pot and free_mean(e1, lam*lam/2) == lam**2*free_mean(e1)
        # minimum over lam of g^2 kin/lam^2 + pot lam^4/g^2 is (3/2)(kin^2 2 pot)^(1/3) g^(2/3)
        cubes[d] = F(27, 8)*kin**2*2*pot
        ok &= cubes[d] == F(729, 256)*d**3*(d - 1)
    out['C6 X has degree four: under C -> lam C the kinetic part goes as lam^-2 and X as lam^4 (free readings, d = 3, 4)'] = ok
    out['C6 Gaussian reading: ground rate <= c_d g^(2/3), c_d^3 = (729/256) d^3 (d - 1) ; 5.35 < c_3 < 5.36 ; 8.17 < c_4 < 8.18'] = \
        F(535, 100)**3 < cubes[3] < F(536, 100)**3 and F(817, 100)**3 < cubes[4] < F(818, 100)**3
    ok = True
    for _ in range(50):
        g2, a, t = (F(rnd.randint(1, 9), rnd.randint(1, 9)) for _ in range(3))
        pp = 2*a*a*g2*g2                         # |p|^2 = 2 a^2 g^4
        # -(g^2/2) psi''/psi + (|p|^2/g^2) t^2 for psi = exp(-a t^2)
        ok &= -(g2/2)*(4*a*a*t*t - 2*a) + (pp/g2)*t*t == g2*a and (g2*a)**2*2 == pp
    out['C6 the layer\'s ground rate |p|/sqrt 2 does not depend on g: the valley is closed by a rise that carries no coupling'] = ok
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out.items():
        print('PASS' if v else 'FAIL', k)
    for k, v in num.items():
        print(k, v)
    with open(os.path.join(HERE, 'TC1_RESULT.json'), 'w') as fh:
        json.dump({'checks': {k: bool(v) for k, v in out.items()}, 'numbers': num}, fh, indent=1)
        fh.write('\n')
