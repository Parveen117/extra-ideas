"""CR1 -- the pure numbers of the core: its two lowest rates in the sector unchanged by both turnings.

The core (TC1): C a free 3 x d matrix, X = e2 of the three numbers x_a (squared singular values),
    h = -(1/2) Laplacian + X ,     rates of -(g^2/2) Laplacian + X/g^2 = g^(2/3) x rates of h.
Readings unchanged by turning the cuts and the directions are functions of e1, e2, e3 (the symmetric functions
of the three numbers).

R1  on such functions (1/2) Laplacian is
        L = 3d D1 + 2(d-1) e1 D2 + (d-2) e2 D3
            + 2 e1 D11 + 8 e2 D12 + 12 e3 D13 + (2 e1 e2 + 6 e3) D22 + 8 e1 e3 D23 + 2 e2 e3 D33 :
    the number of directions enters in three coefficients only.
R2  free means (weight exp(-w e1)) by one rule:  mean(P) = mean(L P)/(2 w m)  for P of degree m (e1: 1, e2: 2, e3: 3).
R3  the ladder: readings (monomial of degree <= D) x exp(-w e1/2); exact matrices S, H; the number of rates of the
    ladder below mu is the count of negative directions of H - mu S, and it bounds the rates of h from above:
        three turns   lowest rate < 5.1868 ,  second < 8.0480
        four turns    lowest rate < 8.0030 ,  second < 11.5660
R4  from below: for every pair -(a/2) Laplacian_y + b |x cross y|^2 >= |x| sqrt(2ab) (the layer's ground rate), and
    -(alpha/2) Laplacian + beta |c| >= 3 (25 alpha beta^2/128)^(1/3) (local rate of exp(-k |c|^(3/2))); so
        lowest rate^3 >= 675 d^3 (d - 1)/512 :   4.14 (three turns),  6.32 (four turns).

Exact rationals.  Stdlib only.
"""
from fractions import Fraction as F
from math import comb, lcm
import importlib.util
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def load(folder, name):
    path = os.path.join(HERE, '..', folder, name + '.py')
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


tc1 = load('tc1', 'tc1_core_of_turns')

RATE = {3: F(8, 5), 4: F(9, 5)}          # w of the ladder's weight exp(-w e1/2)
TOP = 10                                  # degree of the certified ladder


# ---------------------------------------------------------------- polynomials in e1, e2, e3
def padd(p, q, c=1):
    out = dict(p)
    for e, x in q.items():
        out[e] = out.get(e, 0) + c*x
    return {e: x for e, x in out.items() if x}


def pmul(p, q):
    out = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
            out[e] = out.get(e, 0) + c1*c2
    return {e: x for e, x in out.items() if x}


def degree(e):
    return e[0] + 2*e[1] + 3*e[2]


def gen_mono(e, d):
    """L on the monomial e1^a e2^b e3^c"""
    a, b, c = e
    out = {}

    def put(k, v):
        if v:
            out[k] = out.get(k, 0) + v
    put((a - 1, b, c), 3*d*a + 2*a*(a - 1) + 8*a*b + 12*a*c)
    put((a + 1, b - 1, c), 2*(d - 1)*b + 2*b*(b - 1) + 8*b*c)
    put((a, b + 1, c - 1), (d - 2)*c + 2*c*(c - 1))
    put((a, b - 2, c + 1), 6*b*(b - 1))
    return out


def gen(p, d):
    out = {}
    for e, x in p.items():
        for k, v in gen_mono(e, d).items():
            out[k] = out.get(k, 0) + x*v
    return {e: x for e, x in out.items() if x}


def free_means(d, w):
    memo = {(0, 0, 0): F(1)}

    def mono(e):
        if e not in memo:
            memo[e] = sum(v*mono(k) for k, v in gen_mono(e, d).items())/(2*w*degree(e))
        return memo[e]

    def mean(p):
        return sum(x*mono(e) for e, x in p.items())
    return mean


def ladder(d, w, top):
    """exact S and H of the readings m exp(-w e1/2), m a monomial of degree <= top"""
    mean = free_means(d, w)
    basis = [(a, b, c) for a in range(top + 1) for b in range(top//2 + 1) for c in range(top//3 + 1)
             if a + 2*b + 3*c <= top]
    basis.sort(key=lambda e: (degree(e), e))
    n = len(basis)
    s = [[None]*n for _ in range(n)]
    h = [[None]*n for _ in range(n)]
    for j, ej in enumerate(basis):
        m = {ej: F(1)}
        # h (m g)/g = -L m + 2 w deg m + (3 d w/2 - w^2 e1/2) m + e2 m ,  g = exp(-w e1/2)
        hp = padd({}, gen(m, d), -1)
        hp = padd(hp, m, 2*w*degree(ej) + F(3*d)*w/2)
        hp = padd(hp, pmul(m, {(1, 0, 0): F(1)}), -w*w/2)
        hp = padd(hp, pmul(m, {(0, 1, 0): F(1)}))
        for i, ei in enumerate(basis):
            s[i][j] = mean(pmul({ei: F(1)}, m))
            h[i][j] = mean(pmul({ei: F(1)}, hp))
    return basis, s, h


def below(s, h, mu, size=None):
    """number of negative directions of H - mu S, by exact elimination: the count of negative pivots"""
    n = len(s) if size is None else size
    a = [[h[i][j] - mu*s[i][j] for j in range(n)] for i in range(n)]
    neg = 0
    for k in range(n):
        p = a[k][k]
        assert p != 0
        if p < 0:
            neg += 1
        for i in range(k + 1, n):
            if a[i][k]:
                f = a[i][k]/p
                for j in range(i, n):
                    a[i][j] -= f*a[k][j]
        for i in range(k + 1, n):
            for j in range(k + 1, i):
                a[i][j] = a[j][i]
    return neg


def below_by_minors(s, h, mu, size=None):
    """second route (used in the tests): the matrix is made whole by one positive factor and the signs of its
    leading minors are read by fraction-free elimination; the count is the number of sign changes"""
    n = len(s) if size is None else size
    rows = [[h[i][j] - mu*s[i][j] for j in range(n)] for i in range(n)]
    den = 1
    for row in rows:
        for x in row:
            den = lcm(den, x.denominator)
    a = [[int(x*den) for x in row] for row in rows]
    prev, changes, last_sign = 1, 0, 1
    for k in range(n):
        p = a[k][k]                                              # the k-th leading minor (Bareiss)
        assert p != 0
        sign = 1 if p > 0 else -1
        if sign != last_sign:
            changes += 1
        last_sign = sign
        for i in range(k + 1, n):
            for j in range(i, n):
                a[i][j] = (a[i][j]*p - a[i][k]*a[k][j])//prev
            for j in range(k + 1, i):
                a[i][j] = a[j][i]
        prev = p
    return changes


def bracket(s, h, index, lo, hi, size, grid=10**5):
    """the index-th rate of a ladder between two neighbours of the grid 1/grid, by counting"""
    lo, hi = int(F(lo)*grid), int(F(hi)*grid) + 1
    assert below(s, h, F(lo, grid), size) <= index < below(s, h, F(hi, grid), size)
    while hi - lo > 1:
        mid = (lo + hi)//2
        if below(s, h, F(mid, grid), size) > index:
            hi = mid
        else:
            lo = mid
    return F(lo, grid), F(hi, grid)


# ---------------------------------------------------------------- the entries of C (for R1, R2)
def entry_laplacian_half(poly):
    out = {}
    for e, c in poly.items():
        for i, k in enumerate(e):
            if k >= 2:
                f = list(e)
                f[i] -= 2
                f = tuple(f)
                out[f] = out.get(f, 0) + F(c*k*(k - 1), 2)
    return {e: c for e, c in out.items() if c}


def to_entries(p, gens):
    """a polynomial in e1, e2, e3 written in the entries of C"""
    out = {}
    cache = {}

    def power(i, k):
        if (i, k) not in cache:
            cache[(i, k)] = {tuple([0]*len(next(iter(gens[0])))): 1} if k == 0 else tc1.pmul(power(i, k - 1), gens[i])
        return cache[(i, k)]
    for (a, b, c), x in p.items():
        term = tc1.pmul(tc1.pmul(power(0, a), power(1, b)), power(2, c))
        out = tc1.padd(out, term, x)
    return out


# ---------------------------------------------------------------- Laurent polynomials in t = sqrt(r) (for R4)
def ldiff_r(p):
    """d/dr of a Laurent polynomial in t = sqrt(r): (1/(2t)) d/dt"""
    return {k - 2: F(c*k, 2) for k, c in p.items() if k}


def lmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            out[k1 + k2] = out.get(k1 + k2, 0) + c1*c2
    return {k: c for k, c in out.items() if c}


def run():
    out, num = {}, {}
    rnd = random.Random(20261009)
    e1, e2, e3 = {(1, 0, 0): F(1)}, {(0, 1, 0): F(1)}, {(0, 0, 1): F(1)}
    gens_e = [e1, e2, e3]

    # ---- R1: the generator in the three symmetric functions
    for d in (3, 4, 5):
        gens = tc1.entry_polys(d)
        ok = True
        tests = gens_e + [pmul(gens_e[i], gens_e[j]) for i in range(3) for j in range(i, 3)]
        if d == 5:
            tests = tests[:8]                                   # without e3^2 (size); the eight fix every coefficient but D33's
        for p in tests:
            ok &= entry_laplacian_half(to_entries(p, gens)) == to_entries(gen(p, d), gens)
        out['R1 d = %d: (1/2) Laplacian in the entries of C equals L on e_i and e_i e_j (%d polynomial identities)'
            % (d, len(tests))] = ok
    out['R1 the number of directions enters three coefficients only: L e1 = 3d, L e2 = 2(d-1) e1, L e3 = (d-2) e2'] = \
        all(gen(e1, d) == {(0, 0, 0): 3*d} and gen(e2, d) == {(1, 0, 0): 2*(d - 1)} and
            gen(e3, d) == ({(0, 1, 0): d - 2} if d != 2 else {}) for d in (2, 3, 4, 5, 6))

    # ---- R2: free means by one rule
    ok = True
    for d in (3, 4):
        gens = tc1.entry_polys(d)
        for w in (F(1), F(8, 5)):
            mean = free_means(d, w)
            for p in (e1, e2, e3, pmul(e1, e1), pmul(e1, e2), pmul(e2, e2), pmul(e1, e3)):
                ok &= mean(p) == tc1.free_mean(to_entries(p, gens), 1/(2*w))
    mean4 = free_means(4, F(1))
    out['R2 mean(P) = mean(L P)/(2 w m): equal to the means by pairing of the entries (d = 3, 4 ; w = 1, 8/5 ; seven readings)'] = ok
    out['R2 four turns, w = 1: 6, 9, 3, 42, 72 (TC1-C4)'] = \
        [mean4(p) for p in (e1, e2, e3, pmul(e1, e1), pmul(e1, e2))] == [6, 9, 3, 42, 72]

    # ---- R3: the ladder and its counts
    shifts = {3: (F(51868, 10000), F(8048, 1000)), 4: (F(8003, 1000), F(11566, 1000))}
    ladders = {}
    for d in (3, 4):
        w = RATE[d]
        basis, s, h = ladder(d, w, TOP)
        ladders[d] = (basis, s, h)
        n = len(basis)
        out['R3 d = %d: the %d x %d matrices S and H are symmetric (exact)' % (d, n, n)] = \
            all(s[i][j] == s[j][i] and h[i][j] == h[j][i] for i in range(n) for j in range(i))
        out['R3 d = %d: the first reading alone gives 3dw/4 + (3/2)C(d,2)/w^2 (TC1-C6\'s Gaussian reading)' % d] = \
            h[0][0]/s[0][0] == F(3*d)*w/4 + F(3, 2)*comb(d, 2)/(w*w)
        lo, hi = shifts[d]
        out['R3 d = %d: one rate of the ladder below %s and two below %s, so the two lowest rates of h are below them'
            % (d, float(lo), float(hi))] = below(s, h, lo) >= 1 and below(s, h, hi) >= 2 and below(s, h, F(0)) == 0
    rows, last, lows = [], {}, {}
    for d in (3, 4):
        basis, s, h = ladders[d]
        for top in (0, 2, 4, 6, 8, TOP):
            size = sum(1 for e in basis if degree(e) <= top)
            grid = 10**5 if top < TOP else 10**4
            # a larger ladder contains the smaller: its rates are not above the previous ones
            step0, step1 = (F(1, 2), F(1)) if top < TOP else (F(1, 50), F(1, 50))
            hi0 = last.get((d, 0), F(30))
            b0 = bracket(s, h, 0, max(hi0 - step0, 0) if top else 0, hi0, size, grid)
            row = [d, top, size, float(b0[1])]
            last[(d, 0)], lows[(d, top, 0)] = b0[1], b0[0]
            if size > 1:
                hi1 = last.get((d, 1), F(60))
                b1 = bracket(s, h, 1, max(hi1 - step1, 0) if (d, 1) in last else 0, hi1, size, grid)
                row.append(float(b1[1]))
                last[(d, 1)], lows[(d, top, 1)] = b1[1], b1[0]
            rows.append(row)
    num['ladder: turns, degree, size, lowest rate, second rate (upper ends of exact brackets)'] = rows
    out['R3 the ladder only falls as it grows (degrees 0 to 8 on one grid; degree 10 starts below degree 8), d = 3, 4'] = \
        all(rows[i + 1][3] <= rows[i][3] for i in (0, 1, 2, 3, 6, 7, 8, 9)) and \
        all(rows[i + 1][4] <= rows[i][4] for i in (1, 2, 3, 7, 8, 9)) and \
        all(lows[(d, TOP, i)] < F(str(rows[4 + 6*(d - 3)][3 + i])) for d in (3, 4) for i in (0, 1))
    for d, i in ((3, 5), (4, 11)):
        num['d = %d, degree %d: lowest and second rate of the ladder, and their difference' % (d, TOP)] = \
            [rows[i][3], rows[i][4], round(rows[i][4] - rows[i][3], 4)]

    # ---- R4: the lowest rate from below
    ok = True
    for _ in range(50):
        a, b, x, y = (F(rnd.randint(1, 9), rnd.randint(1, 9)) for _ in range(4))
        # -(a/2) psi''/psi + b x^2 y^2 for psi = exp(-m y^2), with 2 m^2 a = b x^2 (m = x sqrt(b/(2a))): rate a m per direction
        m_sq = b*x*x/(2*a)
        # rate^2 = a^2 m^2 = a b x^2/2 ; two transverse directions: (2 a m)^2 = 2 a b x^2
        ok &= -(a/2)*(4*m_sq*y*y) + b*x*x*y*y == 0 and (2*a)**2*m_sq == 2*a*b*x*x
    out['R4 one pair: -(a/2) Laplacian_y + b |x cross y|^2 has ground rate |x| sqrt(2ab) (two transverse directions)'] = ok
    ok = True
    for d in (3, 4, 5):
        alpha = F(1, 2)
        a, b = (1 - alpha)/(d - 1), F(1, 2)
        # every turn is the x of d - 1 pairs: (d - 1) sqrt(2ab) = beta, and gives away (1 - alpha) of its kinetic part
        ok &= (d - 1)**2*2*a*b == (1 - alpha)*(d - 1) and (d - 1)*a == 1 - alpha
    out['R4 summing the pairs: h >= sum over turns of -(alpha/2) Laplacian + beta |c|, beta^2 = (1 - alpha)(d - 1)'] = ok
    k = F(3, 7)
    logd = {-2: F(1), 1: -F(3, 2)*k}                              # psi'/psi for psi = r exp(-k r^(3/2)), in t = sqrt r
    second = padd_l(ldiff_r(logd), lmul(logd, logd))             # psi''/psi
    out['R4 for u = r Exp(-k r^(3/2)): u\'\'/u = (9/4) k^2 r - (15/4) k r^(-1/2) (the same as Laplacian psi/psi in space)'] = \
        second == {2: F(9, 4)*k*k, -1: -F(15, 4)*k}
    ok = True
    for _ in range(50):
        alpha, beta = F(rnd.randint(1, 9), 10), F(rnd.randint(1, 9), rnd.randint(1, 5))
        # local rate (beta - 9 alpha k^2/8) r + (15 alpha k/8) r^(-1/2) = A r + B r^(-1/2) >= 3 (A B^2/4)^(1/3)
        # with k^2 = 4 beta/(9 alpha): A = beta/2 and (bound)^3 = 27 A B^2/4 = 27 . 25 alpha beta^2/128
        k_sq = 4*beta/(9*alpha)
        big_a = beta - 9*alpha*k_sq/8
        b_sq = (15*alpha/8)**2*k_sq
        ok &= big_a == beta/2 and 27*big_a*b_sq/4 == 27*25*alpha*beta*beta/128
        r = F(rnd.randint(1, 30), 7)                               # three-term mean: (A r)(B/2 r^-1/2)^2 = A B^2/4
        ok &= big_a*r*(b_sq/4)/r == big_a*b_sq/4
    out['R4 its least local rate: (rate)^3 >= 27 . 25 alpha beta^2/128 (three-term mean), so with alpha = 1/2: '
        'lowest rate of h cubed >= 675 d^3 (d - 1)/512'] = \
        ok and all(d**3*27*25*F(1, 2)*F(d - 1, 2)/128 == F(675*d**3*(d - 1), 512) for d in (3, 4, 5))
    low = {d: F(675*d**3*(d - 1), 512) for d in (3, 4)}
    out['R4 three turns: 4.14 < lowest rate < 5.1868 ; four turns: 6.32 < lowest rate < 8.0030'] = \
        F(414, 100)**3 < low[3] and F(632, 100)**3 < low[4] and \
        below(ladders[3][1], ladders[3][2], shifts[3][0]) >= 1 and below(ladders[4][1], ladders[4][2], shifts[4][0]) >= 1
    out['R4 control: the bound from below is under the ladder at every degree (it could not be otherwise)'] = \
        all(F(str(row[3]))**3 > low[row[0]] for row in rows)
    num['lowest rate cubed, from below'] = {str(d): str(low[d]) for d in low}
    return out, num


def padd_l(p, q):
    out = dict(p)
    for k, c in q.items():
        out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}


if __name__ == '__main__':
    out, num = run()
    for k, v in out.items():
        print('PASS' if v else 'FAIL', k)
    for k, v in num.items():
        print(k, v)
    with open(os.path.join(HERE, 'CR1_RESULT.json'), 'w') as fh:
        json.dump({'checks': {k: bool(v) for k, v in out.items()}, 'numbers': num}, fh, indent=1)
        fh.write('\n')
