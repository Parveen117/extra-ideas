"""DS1: one diagonal.

The owner's statement: in the drafts, and in the whole framework, the seam is the diagonal -- the 45-degree
line, a local reference.  On it one knows what remained and what is lost; no absolute reference is needed.
The Riemann line, the quantum and classical readings and the physics line have the same structure there.
This stage matches that structure against what is already certified, and proves the common part once.

Sources read before building: RKF theorum/24 Theorem 6.1 (S = R + D, F = Z* j Z = R - D, j = P - Q),
theorum/28 (finite signed cut form F_n = R_n - D_n; "decided by a matrix of size at most five";
dim ker F = dim ker(I - B)), theorum/02 (F = S - V V*, rank V <= 5, B = V* S^-1 V, k_Sigma = N_+(B - I)),
theorum/40 (W = B + A_Gamma - P_1/2; the zero-sum formula pinned, not rederived), physics ONE_LAW, DO1, QD1,
DU1, DG1, FD1, PT1, tools RW1, PM1, CY1, Publications LAM-2 (tau(chi_4)^2 = -4).

Carrier: exact cut-complex rationals (rad + iota turn), dagger = (rad, -turn).  No root is taken.

D1  THE FLIP IS THE MIRROR PAIRING.  For a cut P + Q = 1 put j = P - Q (j^2 = 1, j* = j).  Then
        S = Z* Z = R + D ,     F = Z* j Z = R - D        (theorum/24).
    For a mirror of mirrors (j x j' on a product reading) the flip is the product of the flips.
D2  THE DIAGONAL IS THE 45-DEGREE LINE.  For one reading z = a + b iota read by the cut (a seen, b lost):
        F/S = rad( z / z-dagger ) ,    (F/S)^2 + (turn(z / z-dagger))^2 = 1 ,    4 R D = S^2 - F^2 .
    Seen = lost exactly when z / z-dagger is a quarter turn: z on the 45-degree lines.  The diagonal
    direction is 1 + iota; its norm is the prime 2, its square is 2 iota = the Gauss sum of the character
    mod 4, its fourth power is -4.
D3  POSITIVITY IS A COUNT, AND THE COUNT IS FINITE.  If R is positive and D = V V* has rank r, then
        n_-(R - D) = n_+(B - 1) = number of zeros of det(1 - z B) in 0 < z < 1 ,    B = V* R^-1 V  (r x r),
    and dim ker(R - D) = dim ker(1 - B).  The count is at most r, and does not depend on the size of the
    space.  F >= 0 exactly when no direction has lost/seen above 1: nothing past the diagonal.
D4  THE MIRROR FORM.  For points x with a mirror x -> x', and weights m_x = m_x' > 0,
        W(f) = sum_x m_x f(x) dagger(f(x'))
    is a square at every point on the mirror's fixed line, and one plus and one minus for every pair off it:
        inertia(W) = ( fixed + pairs , pairs ) .       W >= 0  iff  every point is on the line.
    The same for the mirror s -> 1 - s-dagger (line rad = 1/2) and for J(z) = iota z-dagger (the 45-degree line);
    the chart z = (1 - iota)(s - 1/2) carries one into the other.
D5  THE SEAM MEASURE.  For two strands with amplitudes a, b (factorial series, F00-E):
        E(ab)^2 = Exp(-(a - b)^2) E(a^2) E(b^2) ,
        E(a^2) E(b^2) - E(ab)^2 = sum_(i<j) (a^i b^j - a^j b^i)^2 / (i! j!) ,
    so the squared mirror overlap of the strands is Exp(-w^2), w = a - b: 1 on the diagonal, below 1 off it.
    The share of returned histories lies between sinh(2ab)/(2ab) Exp(-a^2 - b^2) and
    (Exp(-(a-b)^2) + Exp(-(a+b)^2))/2.
D7  THE CELL BETWEEN THE DIAGONALS.  On a clock of 4M marks the window of M marks (one cell of the four
    units; a quarter turn wide) has harmonics W(k) = sum_(a<M) zeta^(ka).  W(k) = 0 exactly for the k that are
    multiples of 4 and not of 4M: the cell does not see a unit-free harmonic except the mean (the
    counterpart of PM1-P1).  sum_k |W(k)|^2 = 4 M^2: four cells.
D6  NO CROSSING.  For the returned count of RW1, x = w^2/y satisfies 0 < x < 1 and y dx/dy = w (2V - w)/y > 0:
    the reading x + iota turns from the cut toward the diagonal 1 + iota, always closer, never across.

Python 3.12, standard library only.
"""
from fractions import Fraction as F
from math import comb, factorial
import importlib.util
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def load(folder, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, '..', folder, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class G:
    """exact cut-complex rational  rad + iota*turn"""
    __slots__ = ('r', 't')

    def __init__(self, r=0, t=0):
        self.r, self.t = F(r), F(t)

    def __add__(self, o):
        o = g(o)
        return G(self.r + o.r, self.t + o.t)
    __radd__ = __add__

    def __neg__(self):
        return G(-self.r, -self.t)

    def __sub__(self, o):
        return self + (-g(o))

    def __rsub__(self, o):
        return g(o) - self

    def __mul__(self, o):
        o = g(o)
        return G(self.r*o.r - self.t*o.t, self.r*o.t + self.t*o.r)
    __rmul__ = __mul__

    def dag(self):
        return G(self.r, -self.t)

    def norm(self):
        return self.r*self.r + self.t*self.t

    def inv(self):
        n = self.norm()
        return G(self.r/n, -self.t/n)

    def __truediv__(self, o):
        return self*g(o).inv()

    def __pow__(self, n):
        out = G(1)
        for _ in range(n):
            out = out*self
        return out

    def __eq__(self, o):
        o = g(o)
        return self.r == o.r and self.t == o.t

    def __hash__(self):
        return hash((self.r, self.t))

    def __repr__(self):
        return f'({self.r} + {self.t} iota)'


def g(x):
    return x if isinstance(x, G) else G(x)


IOTA = G(0, 1)


# ---------------------------------------------------------------- matrices over the carrier
def mmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), G(0)) for j in range(len(b[0]))] for i in range(len(a))]


def madd(a, b, s=1):
    return [[a[i][j] + s*b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def dagger(a):
    return [[a[j][i].dag() for j in range(len(a))] for i in range(len(a[0]))]


def eye(n):
    return [[G(1 if i == j else 0) for j in range(n)] for i in range(n)]


def is_hermitian(a):
    return a == dagger(a)


def inverse(a):
    n = len(a)
    m = [row[:] + e for row, e in zip([r[:] for r in a], eye(n))]
    for c in range(n):
        p = next(r for r in range(c, n) if m[r][c] != 0)
        m[c], m[p] = m[p], m[c]
        iv = m[c][c].inv()
        m[c] = [x*iv for x in m[c]]
        for r in range(n):
            if r != c and m[r][c] != 0:
                f = m[r][c]
                m[r] = [x - f*y for x, y in zip(m[r], m[c])]
    return [row[n:] for row in m]


def inertia(a):
    """(n_plus, n_minus, n_zero) of a self-dagger matrix by exact congruences (no root, no float)"""
    a = [row[:] for row in a]
    assert is_hermitian(a), 'not self-dagger'
    n = len(a)
    plus = minus = 0
    live = list(range(n))
    while live:
        i = next((k for k in live if a[k][k] != 0), None)
        if i is None:
            pair = next(((p, q) for p in live for q in live if p < q and a[p][q] != 0), None)
            if pair is None:
                break
            p, q = pair
            c = a[p][q].dag()                               # e_p + c e_q has value 2 |a_pq|^2 > 0
            for k in range(n):
                a[k][p] = a[k][p] + a[k][q]*c
            for k in range(n):
                a[p][k] = a[p][k] + c.dag()*a[q][k]
            i = p
        d = a[i][i]
        assert d.t == 0
        if d.r > 0:
            plus += 1
        else:
            minus += 1
        col = [a[k][i] for k in range(n)]
        row = [a[i][k] for k in range(n)]
        for p in live:
            if p == i:
                continue
            for q in live:
                if q == i:
                    continue
                a[p][q] = a[p][q] - col[p]*row[q]/d
        live.remove(i)
        for k in range(n):
            a[k][i] = G(0)
            a[i][k] = G(0)
    return plus, minus, n - plus - minus


def charpoly_one_minus_zb(b):
    """coefficients c_0..c_r of det(1 - z B) = sum c_k z^k  (Faddeev-LeVerrier), all radial for self-dagger B"""
    r = len(b)
    coeffs = [F(1)]
    m = eye(r)
    c_prev = G(1)
    for k in range(1, r + 1):
        bm = mmul(b, m)
        tr = sum((bm[i][i] for i in range(r)), G(0))
        ck = -tr/k
        assert ck.t == 0
        coeffs.append(ck.r)
        m = madd(bm, [[ck*x for x in row] for row in eye(r)])
    return coeffs                                           # det(1 - zB) = sum coeffs[k] z^k


# ---------------------------------------------------------------- real polynomials: Sturm count with multiplicity
def ptrim(p):
    p = p[:]
    while p and p[-1] == 0:
        p.pop()
    return p


def pdiff(p):
    return ptrim([k*c for k, c in enumerate(p)][1:])


def prem(a, b):
    a = ptrim(a)
    while len(a) >= len(b) and a:
        f = a[-1]/b[-1]
        s = len(a) - len(b)
        a = ptrim([x - (f*b[i - s] if i >= s else 0) for i, x in enumerate(a)][:-1] + [0])
    return a


def pgcd(a, b):
    a, b = ptrim(a), ptrim(b)
    while b:
        a, b = b, prem(a, b)
    return [x/a[-1] for x in a]


def pval(p, x):
    out = F(0)
    for c in reversed(p):
        out = out*x + c
    return out


def distinct_roots(p, lo, hi):
    """number of distinct roots in the open interval (lo, hi); p(lo), p(hi) may vanish"""
    p = ptrim(p)
    if len(p) <= 1:
        return 0
    sq = ptrim(p)
    gcd = pgcd(p, pdiff(p))
    if len(gcd) > 1:                                        # square-free part
        q, rem = [], p[:]
        while len(rem) >= len(gcd):
            f = rem[-1]/gcd[-1]
            q.insert(0, f)
            s = len(rem) - len(gcd)
            rem = [x - (f*gcd[i - s] if i >= s else 0) for i, x in enumerate(rem)][:-1]
        sq = ptrim(q)
    while pval(sq, lo) == 0:                                # remove a root at an end point
        sq = ptrim([sum(sq[k]*comb(k, j)*lo**(k - j) for k in range(j, len(sq))) for j in range(len(sq))][1:])
        sq = ptrim([sum(sq[k]*comb(k, j)*(-lo)**(k - j) for k in range(j, len(sq))) for j in range(len(sq))])
    while pval(sq, hi) == 0:
        sq = ptrim([sum(sq[k]*comb(k, j)*hi**(k - j) for k in range(j, len(sq))) for j in range(len(sq))][1:])
        sq = ptrim([sum(sq[k]*comb(k, j)*(-hi)**(k - j) for k in range(j, len(sq))) for j in range(len(sq))])
    chain = [sq, pdiff(sq)]
    while chain[-1]:
        chain.append([-x for x in prem(chain[-2], chain[-1])])
    chain.pop()

    def changes(x):
        vals = [v for v in (pval(q, x) for q in chain) if v != 0]
        return sum(1 for u, v in zip(vals, vals[1:]) if (u > 0) != (v > 0))
    return changes(lo) - changes(hi)


def roots_with_multiplicity(p, lo, hi):
    p = ptrim(p)
    total = 0
    while len(p) > 1:
        total += distinct_roots(p, lo, hi)
        p = pgcd(p, pdiff(p))
    return total


# ---------------------------------------------------------------- bivariate series (dict (i, j) -> coefficient)
def bmul(p, q, deg):
    out = {}
    for (i, j), x in p.items():
        for (k, l), y in q.items():
            if i + j + k + l <= deg:
                out[(i + k, j + l)] = out.get((i + k, j + l), F(0)) + x*y
    return {k: v for k, v in out.items() if v != 0}


def bexp(p, deg):
    """Exp of a bivariate polynomial without constant term, truncated at total degree deg"""
    out, term, n = {(0, 0): F(1)}, {(0, 0): F(1)}, 0
    while term:
        n += 1
        term = {k: v/n for k, v in bmul(term, p, deg).items()}
        for k, v in term.items():
            out[k] = out.get(k, F(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def cyclotomic(n, cache={}):
    """integer coefficients of the n-th cyclotomic polynomial (low degree first)"""
    if n in cache:
        return cache[n]
    num = [-1] + [0]*(n - 1) + [1]                             # x^n - 1
    for d in range(1, n):
        if n % d == 0:
            den = cyclotomic(d)
            quo = [0]*(len(num) - len(den) + 1)
            rem = num[:]
            for k in range(len(quo) - 1, -1, -1):
                quo[k] = rem[k + len(den) - 1]//den[-1]
                for i, c in enumerate(den):
                    rem[k + i] -= quo[k]*c
            num = quo
    cache[n] = num
    return num


def reduce_sum_of_marks(exponents, n):
    """sum zeta_n^e  as a polynomial in zeta_n of degree below phi(n): exact reduction modulo the n-th
    cyclotomic polynomial"""
    poly = [0]*n
    for e in exponents:
        poly[e % n] += 1
    phi = cyclotomic(n)
    for k in range(n - 1, len(phi) - 2, -1):
        c = poly[k]
        if c:
            for i, d in enumerate(phi):
                poly[k - (len(phi) - 1) + i] -= c*d
    return poly[:len(phi) - 1]


def is_zero_sum_of_marks(exponents, n):
    return all(c == 0 for c in reduce_sum_of_marks(exponents, n))


def run():
    out, num = {}, {}
    rng = random.Random(45)
    gi = lambda lim=3: G(rng.randint(-lim, lim), rng.randint(-lim, lim))

    # D1 -------------------------------------------------------------------------------------------------
    ok = True
    for _ in range(20):
        n, k = 6, 3
        z = [[gi() for _ in range(k)] for _ in range(n)]
        cut = rng.sample(range(n), 3)                                  # P keeps these rows, Q the others
        p = [[G(1 if i == j and i in cut else 0) for j in range(n)] for i in range(n)]
        q = madd(eye(n), p, -1)
        j = madd(p, q, -1)
        s, r, d = mmul(dagger(z), z), mmul(dagger(z), mmul(p, z)), mmul(dagger(z), mmul(q, z))
        f = mmul(dagger(z), mmul(j, z))
        ok &= s == madd(r, d) and f == madd(r, d, -1) and mmul(j, j) == eye(n) and is_hermitian(j)
        ok &= inertia(r)[1] == 0 and inertia(d)[1] == 0
    out['D1_flip_is_the_mirror_pairing'] = ok
    # a mirror of mirrors: for a product reading the flip is the product of the flips
    ok = True
    for _ in range(20):
        z1, z2 = [gi(), gi()], [gi(), gi(), gi()]
        j1 = [G(1), G(-1)]                                              # diagonal mirrors on each factor
        j2 = [G(1), G(-1), G(-1)]
        f1 = sum((x.dag()*s*x for x, s in zip(z1, j1)), G(0))
        f2 = sum((x.dag()*s*x for x, s in zip(z2, j2)), G(0))
        f12 = sum(((x*y).dag()*(s*t)*(x*y) for x, s in zip(z1, j1) for y, t in zip(z2, j2)), G(0))
        ok &= f12 == f1*f2
    out['D1_nested_mirrors_multiply_the_flips'] = ok
    # a cut given by a mirror that is not a coordinate cut: K = swap of two labels, P = (1 + K)/2
    kk = [[G(0), G(1)], [G(1), G(0)]]
    pk = [[x/2 for x in row] for row in madd(eye(2), kk)]
    qk = [[x/2 for x in row] for row in madd(eye(2), kk, -1)]
    out['D1_swap_cut_projects_on_the_two_diagonals'] = (
        mmul(pk, pk) == pk and mmul(qk, qk) == qk and madd(pk, qk, -1) == kk
        and mmul(pk, [[G(1)], [G(1)]]) == [[G(1)], [G(1)]] and mmul(pk, [[G(1)], [G(-1)]]) == [[G(0)], [G(0)]])

    # D2 -------------------------------------------------------------------------------------------------
    ok_turn = ok_unit = ok_diag = True
    for a in range(-6, 7):
        for b in range(-6, 7):
            if (a, b) == (0, 0):
                continue
            z = G(a, b)
            u = z/z.dag()
            seen, lost = F(a*a), F(b*b)
            ok_turn &= u.r == (seen - lost)/(seen + lost)
            ok_unit &= u.norm() == 1 and 4*seen*lost == (seen + lost)**2 - (seen - lost)**2
            ok_diag &= (seen == lost) == (u == IOTA or u == -IOTA)
    out['D2_flip_over_source_is_the_rad_part_of_the_turn'] = ok_turn
    out['D2_unit_relation'] = ok_unit
    out['D2_seen_equals_lost_exactly_at_the_quarter_turn'] = ok_diag
    d = G(1, 1)
    gauss = sum((G(c)*IOTA**a for a, c in ((1, 1), (3, -1))), G(0))      # sum chi_4(a) iota^a
    out['D2_diagonal_vector_and_the_prime_two'] = (
        d.norm() == 2 and d/d.dag() == IOTA and d*d == 2*IOTA == gauss and d**4 == G(-4) == gauss*gauss)
    out['D2_no_rational_unit_point_on_the_diagonal'] = all(
        2*p*p != q*q for p in range(1, 200) for q in range(1, 300))     # 2 is not a rational square (to 300)

    # D3 -------------------------------------------------------------------------------------------------
    ok_three = ok_ker = ok_bound = ok_embed = ok_congr = True
    counts = []
    for trial in range(36):
        n = rng.choice((6, 7, 8))
        r = 1 + trial % 5
        gmat = [[gi(2) for _ in range(n)] for _ in range(n)]
        rr = madd(eye(n), mmul(dagger(gmat), gmat))                    # positive: 1 + G* G
        scale = F(rng.choice((1, 2, 3, 5)), rng.choice((1, 2)))
        v = [[gi(2)*scale for _ in range(r)] for _ in range(n)]
        dd = mmul(v, dagger(v))
        ff = madd(rr, dd, -1)
        bb = mmul(dagger(v), mmul(inverse(rr), v))
        in_f, in_b = inertia(ff), inertia(madd(bb, eye(r), -1))
        poly = charpoly_one_minus_zb(bb)
        zeros = roots_with_multiplicity(poly, F(0), F(1))
        ok_three &= in_f[1] == in_b[0] == zeros
        ok_ker &= in_f[2] == in_b[2]
        ok_bound &= in_f[1] <= r and inertia(rr) == (n, 0, 0)
        counts.append(in_f[1])
        # the two block congruences behind the count (Haynsworth):  [[R, V], [V*, 1]]
        big = [rr[i] + v[i] for i in range(n)] + [dagger(v)[i] + eye(r)[i] for i in range(r)]
        ok_congr &= inertia(big)[1] == in_f[1] == inertia(madd(eye(r), bb, -1))[1]
        # a larger space with nothing lost in the new directions: the count does not move
        extra = 3
        r2 = [rr[i] + [G(0)]*extra for i in range(n)] + [[G(0)]*n + [G(2 if i == k else 0) for k in range(extra)]
                                                          for i in range(extra)]
        v2 = v + [[G(0)]*r for _ in range(extra)]
        ok_embed &= inertia(madd(r2, mmul(v2, dagger(v2)), -1))[1] == in_f[1]
        ok_embed &= mmul(dagger(v2), mmul(inverse(r2), v2)) == bb
    out['D3_three_routes_to_the_count'] = ok_three
    out['D3_kernels_agree'] = ok_ker
    out['D3_count_is_at_most_the_rank_of_the_memory'] = ok_bound
    out['D3_bordered_form_carries_the_same_count'] = ok_congr
    out['D3_count_does_not_depend_on_the_size_of_the_space'] = ok_embed
    out['D3_instances_cover_several_counts'] = len(set(counts)) >= 4
    num['D3_counts_seen'] = sorted(set(counts))
    # on the diagonal exactly: lost = seen in one direction, a kernel and no count
    rdiag = [[G(4 if i == j == 0 else (1 if i == j else 0)) for j in range(4)] for i in range(4)]
    vdiag = [[G(2)], [G(0)], [G(0)], [G(0)]]
    fd = madd(rdiag, mmul(vdiag, dagger(vdiag)), -1)
    bd = mmul(dagger(vdiag), mmul(inverse(rdiag), vdiag))
    out['D3_on_the_diagonal_a_kernel_and_no_count'] = (
        bd == [[G(1)]] and inertia(fd) == (3, 0, 1) and roots_with_multiplicity(charpoly_one_minus_zb(bd), F(0), F(1)) == 0)
    out['D3_just_past_the_diagonal_one_count'] = inertia(
        madd(rdiag, mmul([[G(F(201, 100))], [G(0)], [G(0)], [G(0)]], [[G(F(201, 100)), G(0), G(0), G(0)]]), -1))[1] == 1

    # D4 -------------------------------------------------------------------------------------------------
    def mirror_form(points, mirror, weights, size):
        m = [[G(0)]*size for _ in range(size)]
        for x, w in zip(points, weights):
            xm = mirror(x).dag()
            px = [x**i for i in range(size)]
            pm = [xm**i for i in range(size)]
            for i in range(size):
                for jx in range(size):
                    m[i][jx] = m[i][jx] + w*px[i]*pm[jx]
        return m
    half = lambda s: G(1 - s.r, s.t)                                   # s -> 1 - s-dagger ; fixed line rad = 1/2
    swap = lambda s: G(s.t, s.r)                                       # J(z) = iota z-dagger ; fixed line rad = turn
    ok_half = ok_swap = True
    for fixed, pairs in ((3, 0), (2, 1), (0, 2), (1, 3), (4, 2)):
        pts, wts = [], []
        for k in range(fixed):
            pts.append(G(F(1, 2), k + 1))
            wts.append(F(k + 1))
        for k in range(pairs):
            a = G(F(1, 2) + F(k + 1, 5), 2*k + 1)
            pts += [a, half(a)]
            wts += [F(k + 2), F(k + 2)]
        size = len(pts) + 2
        m = mirror_form(pts, half, wts, size)
        ok_half &= is_hermitian(m) and inertia(m) == (fixed + pairs, pairs, size - len(pts))
        pts2, wts2 = [], []
        for k in range(fixed):
            pts2.append(G(k + 1, k + 1))
            wts2.append(F(k + 1))
        for k in range(pairs):
            a = G(k + 1, -(k + 2))
            pts2 += [a, swap(a)]
            wts2 += [F(k + 2), F(k + 2)]
        m2 = mirror_form(pts2, swap, wts2, size)
        ok_swap &= is_hermitian(m2) and inertia(m2) == (fixed + pairs, pairs, size - len(pts2))
    out['D4_mirror_form_on_the_half_line'] = ok_half
    # the two mirrors are one: z = (1 - iota)(s - 1/2) carries s -> 1 - s-dagger into J(z) = iota z-dagger
    chart = lambda s: G(1, -1)*(s - F(1, 2))
    out['D4_the_two_mirrors_are_one_in_the_chart_of_the_diagonal_vector'] = all(
        chart(half(G(F(p, 7), F(q, 5)))) == IOTA*chart(G(F(p, 7), F(q, 5))).dag()
        for p in range(-9, 10) for q in range(-9, 10)) and all(
        swap(chart(G(F(1, 2), F(q, 3)))) == chart(G(F(1, 2), F(q, 3))) for q in range(-6, 7))
    out['D4_mirror_form_on_the_45_degree_line'] = ok_swap
    # one pair by hand:  m (f(x) g + g f(x))  =  (m/2) (|f + g|^2 - |f - g|^2)
    ok = True
    for _ in range(50):
        fx, fy = gi(5), gi(5)
        lhs = fx*fy.dag() + fy*fx.dag()
        ok &= lhs.t == 0 and 2*lhs.r == (fx + fy).norm() - (fx - fy).norm()
    out['D4_a_pair_is_one_plus_and_one_minus'] = ok
    # moving a point off the line adds exactly one to the count
    on = [G(F(1, 2), 1), G(F(1, 2), 2), G(F(1, 2), 3)]
    off = [G(F(1, 2), 1), G(F(1, 2), 2), G(F(7, 10), 3), half(G(F(7, 10), 3))]
    out['D4_one_point_off_the_line_is_one_count'] = (
        inertia(mirror_form(on, half, [F(1)]*3, 6))[1] == 0 and inertia(mirror_form(off, half, [F(1)]*4, 6))[1] == 1)
    out['D4_unequal_weights_on_a_pair_are_refused'] = not is_hermitian(
        mirror_form([G(F(7, 10), 3), half(G(F(7, 10), 3))], half, [F(1), F(2)], 4))

    # D5 -------------------------------------------------------------------------------------------------
    deg = 16
    ea2 = {(2*n, 0): F(1, factorial(n)) for n in range(deg//2 + 1)}
    eb2 = {(0, 2*n): F(1, factorial(n)) for n in range(deg//2 + 1)}
    eab = {(n, n): F(1, factorial(n)) for n in range(deg//2 + 1)}
    gauss_w = bexp({(2, 0): F(-1), (1, 1): F(2), (0, 2): F(-1)}, deg)          # Exp(-(a - b)^2)
    out['D5_seam_measure_identity'] = bmul(eab, eab, deg) == bmul(bmul(ea2, eb2, deg), gauss_w, deg)
    sos = {}
    for i in range(deg//2 + 1):
        for jx in range(i + 1, deg + 1):
            if 2*(i + jx) > deg:
                break
            sq = bmul({(i, jx): F(1), (jx, i): F(-1)}, {(i, jx): F(1), (jx, i): F(-1)}, deg)
            for k, v in sq.items():
                sos[k] = sos.get(k, F(0)) + v/(factorial(i)*factorial(jx))
    gram = dict(bmul(ea2, eb2, deg))
    for k, v in bmul(eab, eab, deg).items():
        gram[k] = gram.get(k, F(0)) - v
    out['D5_gram_defect_is_a_sum_of_squares'] = (
        {k: v for k, v in gram.items() if v != 0} == {k: v for k, v in sos.items() if v != 0})
    out['D5_returned_share_termwise_bounds'] = all(
        F(4**n, factorial(2*n + 1)) <= F(1, factorial(n)**2) <= F(4**n, factorial(2*n)) for n in range(0, 200))
    rw1, pm1 = load('rw1', 'rw1_returned_winding_count'), load('pm1', 'pm1_prime_turn_series')
    ok, shares = True, {}
    for a, b in ((F(1), F(1)), (F(2), F(2)), (F(3, 2), F(1, 2)), (F(3), F(1)), (F(5, 2), F(2))):
        (zl, zh), _ = rw1.enclose(0, a*a*b*b)
        ea, eb = pm1.exp_neg(a*a), pm1.exp_neg(b*b)
        share = (zl*ea[0]*eb[0], zh*ea[1]*eb[1])
        upper = (pm1.exp_neg((a - b)**2)[0] + pm1.exp_neg((a + b)**2)[0])/2
        x = 2*a*b
        sinh_over_x = sum(x**(2*n)/factorial(2*n + 1) for n in range(60))     # a lower bound: positive terms
        lower = sinh_over_x*ea[0]*eb[0]
        ok &= lower <= share[0] and share[1] <= upper
        shares[f'a={a},b={b}'] = [round(float(lower), 6), round(float(share[0]), 6), round(float(upper), 6)]
    out['D5_returned_share_between_its_bounds'] = ok
    num['D5_lower_share_upper'] = shares
    out['D5_measure_is_one_on_the_diagonal_only'] = (
        pm1.exp_neg(F(0)) == (F(1), F(1)) and pm1.exp_neg(F(1, 100))[1] < 1)

    # D6 -------------------------------------------------------------------------------------------------
    order = 30
    ok = True
    for nu in range(0, 5):
        w = rw1.mean_series(nu, order + 1)
        v = rw1.theta(w)
        x = rw1.smul(w, w, order + 1)[1:]                              # w^2 / y
        lhs = rw1.theta(x)[:order]                                     # y dx/dy
        rhs = rw1.smul(w, rw1.sadd(rw1.sscale(v, 2), w, -1), order + 1)[1:][:order]   # w (2V - w)/y
        ok &= lhs == rhs
    out['D6_rate_of_approach'] = ok
    ok, xs = True, []
    for nu in (0, 1, 3):
        prev = F(0)
        for y in [F(1, 10), F(1), F(10), F(100), F(1000), F(10000)]:
            wl, wh = rw1.mean_interval(nu, y)
            lo, hi = wl*wl/y, wh*wh/y
            ok &= prev < lo and hi < 1
            prev = hi
            if nu == 0:
                xs.append(round(float(lo), 6))
    out['D6_always_closer_never_across'] = ok
    num['D6_x_for_nu0'] = xs

    # D7 -------------------------------------------------------------------------------------------------
    ok_blind = ok_four = True
    for cell in (1, 2, 3, 5, 6):
        marks = 4*cell
        for k in range(marks):
            zero = is_zero_sum_of_marks([k*a for a in range(cell)], marks)
            ok_blind &= zero == (k % 4 == 0 and k != 0)
        # sum_k |W(k)|^2 = sum over k, a, b of zeta^(k(a-b)): reduced exactly, it is the whole number 4 cell^2
        total = reduce_sum_of_marks([k*(a - b) for k in range(marks) for a in range(cell) for b in range(cell)], marks)
        ok_four &= total[0] == 4*cell*cell and all(c == 0 for c in total[1:])
    out['D7_cell_is_blind_to_unit_free_harmonics_except_the_mean'] = ok_blind
    out['D7_four_cells'] = ok_four
    out['D7_cyclotomic_reduction_is_exact'] = (
        cyclotomic(8) == [1, 0, 0, 0, 1] and cyclotomic(12) == [1, 0, -1, 0, 1]
        and is_zero_sum_of_marks([0, 4], 8) and not is_zero_sum_of_marks([0, 2], 8))

    # controls ----------------------------------------------------------------------------------------------
    out['control_inertia_on_known_forms'] = (
        inertia([[G(0), G(1)], [G(1), G(0)]]) == (1, 1, 0)
        and inertia([[G(0), IOTA], [-IOTA, G(0)]]) == (1, 1, 0)
        and inertia([[G(2), G(1, 1)], [G(1, -1), G(1)]]) == (1, 0, 1)
        and inertia([[G(0)]*3]*3) == (0, 0, 3))
    out['control_root_count_with_multiplicity'] = (
        roots_with_multiplicity([F(1, 4), F(-1), F(1)], F(0), F(1)) == 2          # (z - 1/2)^2
        and roots_with_multiplicity([F(-1), F(1)], F(0), F(1)) == 0               # root at the end point 1
        and roots_with_multiplicity([F(3, 8), F(-2), F(3), F(-1)], F(0), F(1)) == 2)   # roots 1/2, 0.349..., 2.151...
    # the count needs the seen form: with R = diag(4, 1, 1) and v = (3/2, 0, 0), lost/seen = 9/16 (no count),
    # while v* v alone is 9/4 (would count one)
    r3 = [[G(4 if i == j == 0 else (1 if i == j else 0)) for j in range(3)] for i in range(3)]
    v3 = [[G(F(3, 2))], [G(0)], [G(0)]]
    b_right = mmul(dagger(v3), mmul(inverse(r3), v3))
    b_wrong = mmul(dagger(v3), v3)
    out['control_count_needs_the_seen_form'] = (
        b_right == [[G(F(9, 16))]] and inertia(madd(r3, mmul(v3, dagger(v3)), -1))[1] == 0
        and inertia(madd(b_right, eye(1), -1))[0] == 0 and inertia(madd(b_wrong, eye(1), -1))[0] == 1)

    out['pass'] = all(v for v in out.values() if isinstance(v, bool))
    return out, num


if __name__ == '__main__':
    out, num = run()
    with open(os.path.join(HERE, 'DS1_RESULT.json'), 'w') as fh:
        json.dump({'checks': out, 'numbers': num}, fh, indent=1, sort_keys=True)
        fh.write('\n')
    for k, v in out.items():
        print(f'{k}: {v}')
    print(num)
