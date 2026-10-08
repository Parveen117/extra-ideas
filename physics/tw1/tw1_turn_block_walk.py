"""TW1 -- the turn block's walk: one diagonal step is a boost by the pair's flip.

Record (AG1-A6): a closed surface of the turn block with the heat weight, q = Exp(-pi u):

    K+ = sum m^2 q^(m^2)                 source          (m = 1, 2, 3, ...)
    K- = sum (-1)^(m-1) m^2 q^(m^2)      flip            rho = K-/K+
    a = sum q^(n^2) ,  b = sum (-1)^n q^(n^2)   (n whole) :  one turn of the pair ;  S = a^2, F = b^2 of AG1,
    r = b/a ;  lost part L, L^2 = S^2 - F^2 = a^4 - b^4 ;  after n diagonal steps L_n = L at q^(2^n).

W1  one diagonal step (q -> q^2) is exact:
        4 a(q^2) K+(q^2) = a K+ - b K-          4 b(q^2) K-(q^2) = a K- - b K+
    so  rho' r' = (rho - r)/(1 - r rho) :  rapidities subtract.
W2  the record is written by the lost parts along the walk:
        16 K- = b (L_0^2 + sum_(n>=1) 2^n L_n^2)        16 K+ = a (L_0^2 - sum_(n>=1) 2^n L_n^2)
    invariant of the walk:  (K-/b + K+/a)/L^2 = 1/8.
W3  positivity is carried by the walk from the strong side (no mirror description):  0 < r < rho < 1 at every
    coupling;  (3r + r^3)/(1 + 3r^2) < rho < (r + r')/(1 + r r').
W4  the weak end by the backward walk:  beta = rho/r obeys  beta + 1 = 2(beta' + 1)/(1 + delta),
    delta = 2 r^2 beta'/(1 + r^2) ;  the limit of (beta + 1)/2^n agrees with the mirror's value.
W5  rho'(1 + rho) > 2 sqrt(rho) at every coupling:
        strong side  q <= 7/10 :  k(q^2)^2 k(q^4)^2 > o(q)^2 o(q^2) k(q^8)   (series with one sign, a grid)
        weak side    q >= 7/10 :  beta >= 99/10                              (base boxes, backward walk)
W6  consequences and controls.

Exact rational arithmetic; directed fixed point only in the enclosures of W2-W4 and W6.  Stdlib only.
"""
from fractions import Fraction as F
import importlib.util
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def load(folder, name):
    path = os.path.join(HERE, '..', folder, name + '.py')
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ag1 = load('ag1', 'ag1_diagonal_walk')
pm1 = ag1.pm1

ORDER = 600
TERMS = 16                      # contents kept in every rational sum; the rest is bounded
Q_FIRST = F(3, 5)               # strong side: one-sign bound up to here
Q_MEET = F(7, 10)               # strong side up to here, weak side from here
GRID = 400                      # strong side grid step 1/GRID between Q_FIRST and Q_MEET
BOXES = 100                     # weak side: boxes of the base [Q_MEET^4, Q_MEET^2]
BETA_MIN = F(99, 10)


# ---------------------------------------------------------------- sparse integer series in q
def sparse(order, coef, parity=None):
    """sum over m >= 1 of coef(m) q^(m^2), as {exponent: integer}"""
    out, m = {}, 1
    while m*m <= order:
        if parity is None or m % 2 == parity:
            out[m*m] = coef(m)
        m += 1
    return out


def smul(a, b, order=ORDER):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= order:
                out[i + j] = out.get(i + j, 0) + x*y
    return {k: v for k, v in out.items() if v}


def sadd(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c*v
    return {k: v for k, v in out.items() if v}


def sscale(a, c):
    return {k: c*v for k, v in a.items()}


def dil(a, j, order=ORDER):
    return {j*k: v for k, v in a.items() if j*k <= order}


def series(order=ORDER):
    one = {0: 1}
    a = sadd(one, sparse(order, lambda m: 2))
    b = sadd(one, sparse(order, lambda m: 2*(-1)**m))
    kp = sparse(order, lambda m: m*m)
    km = sparse(order, lambda m: (-1)**(m - 1)*m*m)
    od = sparse(order, lambda m: m*m, parity=1)
    return a, b, kp, km, od


def lattice_even(order=ORDER):
    """sum over whole (m, n), m + n even, of m^2 q^(m^2+n^2), and the same with the mark (-1)^m"""
    plain, marked = {}, {}
    top = math.isqrt(order)
    for m in range(-top, top + 1):
        for n in range(-top, top + 1):
            e = m*m + n*n
            if e <= order and (m + n) % 2 == 0 and m:
                plain[e] = plain.get(e, 0) + m*m
                marked[e] = marked.get(e, 0) + (-1)**m*m*m
    return plain, marked


# ---------------------------------------------------------------- rational sums with a bounded rest
def rest(q, power):
    """upper bound of sum over m > TERMS of m^power q^(m^2): first term over (1 - ratio), ratio decreasing in m"""
    m = TERMS + 1
    ratio = F(m + 1, m)**power*q**(2*m + 1)
    assert ratio < 1
    return m**power*q**(m*m)/(1 - ratio)


def part(q, power, parity=None):
    return sum(m**power*q**(m*m) for m in range(1, TERMS + 1) if parity is None or m % 2 == parity)


def k_lo(x):
    """k(x) = K+(x)/x = 1 + 4x^3 + 9x^8 + ... from below"""
    return part(x, 2)/x


def k_hi(x):
    return (part(x, 2) + rest(x, 2))/x


def o_hi(x):
    """o(x) = (odd contents of K+)/x = 1 + 9x^8 + 25x^24 + ... from above"""
    return (part(x, 2, 1) + rest(x, 2))/x


def left_lo(q):
    return k_lo(q*q)**2*k_lo(q**4)**2


def right_hi(q):
    return o_hi(q)**2*o_hi(q*q)*k_hi(q**8)


def box(q1, q2):
    """enclosures of r = b/a and beta = rho/r for every q in [q1, q2]; each partial sum rises with q"""
    od0 = (part(q1, 0, 1), part(q2, 0, 1) + rest(q2, 0))
    ev0 = (part(q1, 0, 0), part(q2, 0, 0) + rest(q2, 0))
    od2 = (part(q1, 2, 1), part(q2, 2, 1) + rest(q2, 2))
    ev2 = (part(q1, 2, 0), part(q2, 2, 0) + rest(q2, 2))
    a = (1 + 2*(ev0[0] + od0[0]), 1 + 2*(ev0[1] + od0[1]))
    b = (1 + 2*(ev0[0] - od0[1]), 1 + 2*(ev0[1] - od0[0]))
    kp = (od2[0] + ev2[0], od2[1] + ev2[1])
    km = (od2[0] - ev2[1], od2[1] - ev2[0])
    assert b[0] > 0 and km[0] > 0
    r = (b[0]/a[1], b[1]/a[0])
    beta = (km[0]*a[0]/(kp[1]*b[1]), km[1]*a[1]/(kp[0]*b[0]))
    return r, beta


def loss(r_hi, beta_hi, steps=6):
    """upper bound of the product of (1 + delta_j) over every backward step j >= 1 below a box.

    r_j <= R_j with R_j = R_(j-1)^2 (1 + R_(j-1)^4)/2 ;  beta_j <= B_j = 2 B_(j-1) + 1 ;
    delta_j <= D_j = 2 R_j^2 B_(j-1) ;  D_(j+1)/D_j <= 3 R_j^2.
    """
    assert 0 < r_hi < 1 and beta_hi >= 1
    big_r, big_b, prod, d = pm1.up(r_hi), pm1.up(beta_hi), F(1), None       # every line rounded upward
    for _ in range(steps):
        big_r = pm1.up(big_r**2*(1 + big_r**4)/2)
        d = pm1.up(2*big_r**2*big_b)
        prod = pm1.up(prod*(1 + d))
        big_b = 2*big_b + 1
    ratio = pm1.up(3*big_r**2)
    assert ratio <= F(1, 2)
    tail = pm1.up(d*ratio/(1 - ratio))
    assert tail < 1
    return pm1.up(prod/(1 - tail))


# ---------------------------------------------------------------- enclosures (W4)
dn, up = pm1.down, pm1.up


def half_step(r, beta):
    """(r, beta) at q, as enclosures  ->  (r, beta) at sqrt(q)"""
    r4 = (dn(r[0]**4), up(r[1]**4))
    root = (ag1.isqrt_iv((1 - r4[1], 1 - r4[1]))[0], ag1.isqrt_iv((1 - r4[0], 1 - r4[0]))[1])
    x = (dn(r[0]**2/(1 + root[1])), up(r[1]**2/(1 + root[0])))            # x = r^2/(1 + sqrt(1 - r^4)), rising in r

    def new(b, y):
        return 2*(b + 1)/(1 + 2*y*y*b/(1 + y*y)) - 1                      # rising in b, falling in y
    return x, (dn(new(beta[0], x[1])), up(new(beta[1], x[0])))


def at_strong(u, pi_iv):
    """enclosures of r and beta at heat time u from the description by contents"""
    s, f, _ = ag1.sfl(F(u), pi_iv)
    r = ag1.isqrt_iv(ag1.idiv(f, s))
    kp, km = ag1.k_contents(F(u), pi_iv)
    rho = ag1.idiv(km, kp)
    return r, ag1.idiv(rho, r), rho


def run():
    out, num = {}, {}
    rnd = random.Random(20261009)
    a, b, kp, km, od = series()
    a2, b2, kp2, km2, od2 = (dil(x, 2) for x in (a, b, kp, km, od))

    # ---- W1: the step
    out['W1 4 a(q^2) K+(q^2) = a K+ - b K-  (series to order 600)'] = \
        sscale(smul(a2, kp2), 4) == sadd(smul(a, kp), smul(b, km), -1)
    out['W1 4 b(q^2) K-(q^2) = a K- - b K+  (series to order 600)'] = \
        sscale(smul(b2, km2), 4) == sadd(smul(a, km), smul(b, kp), -1)
    plain, marked = lattice_even()
    out['W1 lattice: a K+ - b K- is the sum of m^2 over the diagonal sublattice'] = \
        sadd(smul(a, kp), smul(b, km), -1) == plain
    out['W1 lattice: a K- - b K+ is minus that sum with the mark (-1)^m'] = \
        sadd(smul(a, km), smul(b, kp), -1) == sscale(marked, -1)
    out['W1 the pair itself: 2 a(q^2)^2 = a^2 + b^2 and b(q^2)^2 = a b (AG1-A4)'] = \
        sscale(smul(a2, a2), 2) == sadd(smul(a, a), smul(b, b)) and smul(b2, b2) == smul(a, b)

    def moebius(r, rho):                     # the step in ratios, from the two lines of W1
        return (rho - r)/(1 - r*rho)
    ok = True
    for _ in range(200):
        x, y, s, t = (F(rnd.randint(1, 99), rnd.randint(100, 200)) for _ in range(4))   # a, b, K+, K- as numbers
        if x*s == y*t or x*t == y*s:
            continue
        lhs = ((x*t - y*s)/4)/((x*s - y*t)/4)                 # (b' K-')/(a' K+') = rho' r'
        ok &= lhs == moebius(y/x, t/s)
    out['W1 in ratios: rho\' r\' = (rho - r)/(1 - r rho)'] = ok

    # ---- W2: the record written by the lost parts along the walk
    def fourth(x):
        return smul(smul(x, x), smul(x, x))
    lost = []                                                   # L_n^2 = a^4 - b^4 at q^(2^n), while it reaches the order
    n = 0
    while 2**n <= ORDER:
        lost.append(dil(sadd(fourth(a), fourth(b), -1), 2**n))
        n += 1
    later = {}
    for n in range(1, len(lost)):
        later = sadd(later, lost[n], 2**n)
    ab = smul(a, b)
    out['W2 invariant: 8 (a K- + b K+) = a b L_0^2  (series to order 600)'] = \
        sscale(sadd(smul(a, km), smul(b, kp)), 8) == smul(ab, lost[0])
    out['W2 8 (a K- - b K+) = a b (2 L_1^2 + 4 L_2^2 + 8 L_3^2 + ...)  (series to order 600)'] = \
        sscale(sadd(smul(a, km), smul(b, kp), -1), 8) == smul(ab, later)
    out['W2 16 K- = b (L_0^2 + later) and 16 K+ = a (L_0^2 - later)'] = \
        sscale(km, 16) == smul(b, sadd(lost[0], later)) and sscale(kp, 16) == smul(a, sadd(lost[0], later, -1))
    out['W2 every lost part is a series of one sign, beginning 16 q^(2^n)'] = \
        all(min(x.values()) > 0 and x[2**i] == 16 and min(x) == 2**i for i, x in enumerate(lost))
    out['W2 one step of the lost part: 4 L_1^2 = (S - F)^2'] = \
        sscale(lost[1], 4) == smul(sadd(smul(a, a), smul(b, b), -1), sadd(smul(a, a), smul(b, b), -1))
    ok = True
    for _ in range(200):
        r, beta = F(rnd.randint(1, 99), 100), F(rnd.randint(101, 900), 100)
        if r*r*beta >= 1:
            continue
        nxt = (1 + r*r)*(beta - 1)/(2*(1 - r*r*beta))                       # beta' from W1 and r'^2 = 2r/(1+r^2)
        lam, lam_next = (beta - 1)/(beta + 1), (nxt - 1)/(nxt + 1)
        ok &= lam == (1 - r*r)/(2*(1 + r*r))*(1 + lam_next)                 # 2 L_1^2/L_0^2 = (S - F)/(2 (S + F))
    out['W2 the step in Lambda = (rho - r)/(rho + r) is affine: Lambda = (S - F)/(2 (S + F)) (1 + Lambda at q^2)'] = ok

    # ---- W3: positivity from the strong side, and the two bounds
    half = F(1, 2)
    out['W3 base: for q <= 1/2 the terms of b and of K- fall (2q <= 1, 4q^3 <= 1/2), so b > 1 - 2q >= 0, K- > q - 4q^4 > 0'] = \
        2*half <= 1 and 4*half**3 <= half and half - 4*half**4 > 0 and \
        all(F(m + 1, m)**2*half**(2*m + 1) <= 4*half**3 for m in range(1, 40))
    ok = True
    for _ in range(200):
        r, beta = F(rnd.randint(1, 99), 100), F(rnd.randint(101, 900), 100)
        if r*r*beta >= 1:
            continue
        nxt = (1 + r*r)*(beta - 1)/(2*(1 - r*r*beta))                       # beta' from W1 and r'^2 = 2r/(1+r^2)
        ok &= (nxt > 1) == (beta*(1 + 3*r*r) > 3 + r*r)                     # rho' > r'  <=>  rho > tanh 3 gamma
        ok &= (3*r + r**3)*(1 + r*r) - 2*r*(1 + 3*r*r) == r*(1 - r*r)**2    # tanh 3 gamma >= pair flip at q^2
        ok &= beta + 1 == 2*(nxt + 1)/(1 + 2*r*r*nxt/(1 + r*r))             # W4's backward step
    out['W3 rho\' > r\' is rho > (3r + r^3)/(1 + 3r^2), which is above the pair\'s flip at q^2'] = ok

    pi_iv = pm1.pi_interval()
    ok, rows = True, []
    for u in (F(1, 4), F(1, 2), F(1), F(2)):
        r, beta, rho = at_strong(u, pi_iv)
        r2, _, rho2 = at_strong(2*u, pi_iv)
        low = ((3*r[1] + r[1]**3)/(1 + 3*r[0]**2))
        high = (r[0] + r2[0])/(1 + r[1]*r2[1])
        step = ag1.imul(rho2, r2)
        law = ag1.idiv(ag1.isub(rho, r), ag1.isub(ag1.iv(1), ag1.imul(r, rho)))
        ok &= 0 < r[0] and r[1] < rho[0] and rho[1] < 1 and low < rho[0] and rho[1] < high and ag1.overlap(step, law)
        rows.append([str(u), float(r[0]), float(rho[0]), float(low), float(high)])
    out['W3 at u = 1/4, 1/2, 1, 2: r < (3r+r^3)/(1+3r^2) < rho < (r+r\')/(1+r r\') < 1 and the step law (enclosures)'] = ok
    num['u, r, rho, lower bound, upper bound'] = rows

    # W2 at couplings: seen R = L_0^2, memory D = sum 2^n L_n^2 ; lost/seen = Lambda < 1 from finitely many terms
    def lam_of(u, steps):
        terms = []
        for n in range(steps + 1):
            s_, f_, _ = ag1.sfl(u*2**n, pi_iv)
            terms.append(ag1.iscale(ag1.isub(ag1.imul(s_, s_), ag1.imul(f_, f_)), 2**n))
        total = ag1.iv(0)
        for t in terms[1:]:
            total = ag1.iadd(total, t)
        total = (total[0], total[1] + terms[-1][1])                          # the rest is below the last term (ratio <= 1/2)
        return ag1.idiv(total, terms[0]), [float(t[0]/total[0]) for t in terms[1:5]]
    u = F(1, 4)
    r, beta, _ = at_strong(u, pi_iv)
    lam, shares = lam_of(u, 12)
    from_lost = ag1.idiv(ag1.iadd(ag1.iv(1), lam), ag1.isub(ag1.iv(1), lam))
    out['W2 at u = 1/4: (1 + Lambda)/(1 - Lambda) from twelve lost parts is rho/r (60 digits)'] = \
        ag1.overlap(from_lost, beta) and ag1.width(from_lost) < F(1, 10**60) and 0 < lam[0] and lam[1] < 1
    num['u = 1/4: Lambda, and the shares of its first four terms'] = [float(lam[0])] + shares
    ok, rows = True, []
    for u in (F(1), F(1, 4), F(1, 16), F(1, 64)):
        lam, _ = lam_of(u, 16)
        ok &= 0 < lam[0] and lam[1] < 1 and ag1.width(lam) < F(1, 10**40)
        rows.append([str(u), float(lam[0]), float(1 - lam[1]), float(4*u/pi_iv[0])])
    out['W2 lost/seen = Lambda is below 1 at u = 1, 1/4, 1/16, 1/64 (sixteen terms and a bounded rest)'] = ok
    num['u, Lambda, 1 - Lambda, 4u/pi'] = rows

    # ---- W4: the weak end by the backward walk, against the mirror description
    r, beta, _ = at_strong(1, pi_iv)
    track = []
    for n in range(1, 9):
        r, beta = half_step(r, beta)
        track.append((r, beta))
    r4, beta4 = track[3]                                                     # u = 1/16
    walked = ag1.imul(r4, beta4)
    kpw, kmw = ag1.k_windings(F(1, 16), F(1, 64), pi_iv)
    mirror = ag1.idiv(kmw, kpw)
    out['W4 the flip at u = 1/16 walked back from u = 1 agrees with the mirror description (80 digits)'] = \
        ag1.overlap(walked, mirror) and ag1.width(walked) < F(1, 10**80) and ag1.width(mirror) < F(1, 10**80)
    num['flip at u = 1/16'] = float(walked[0])
    r8, beta8 = track[7]                                                     # u = 1/256
    top = (beta8[1] + 1)/256
    bottom = (beta8[0] + 1)/256/loss(r8[1], beta8[1])
    half_pi = (pi_iv[0]/2, pi_iv[1]/2)
    out['W4 the limit of (beta + 1)/2^n from u = 1 is pi/2 to 80 digits (the mirror\'s value)'] = \
        bottom <= half_pi[1] and half_pi[0] <= top and top - bottom < F(1, 10**80)
    num['limit of u (beta + 1)'] = float(top)
    ok = all(track[i][0][1] <= track[i - 1][0][0]**2 for i in range(1, 8))
    out['W4 the pair\'s ratio at least squares at every backward step: r(u/2) <= r(u)^2'] = ok

    # ---- W5: the inequality as one line of series
    big = sadd(smul(smul(km2, km2), smul(sadd(kp, km), sadd(kp, km))),
               sscale(smul(smul(kp2, kp2), smul(kp, km)), 4), -1)
    kp4, kp8 = dil(kp, 4), dil(kp, 8)
    prod_form = sadd(smul(smul(kp2, kp2), smul(kp4, kp4)), smul(smul(od, od), smul(od2, kp8)), -1)
    out['W5 K-(q^2)^2 (K+ + K-)^2 - 4 K+(q^2)^2 K+ K- = 64 (K+(q^2)^2 K+(q^4)^2 - O^2 O(q^2) K+(q^8))'] = \
        big == sscale(prod_form, 64)
    first = min(big)
    negative = sum(1 for v in big.values() if v < 0)
    out['W5 the difference begins 512 q^18 - 1152 q^20 : not of one sign term by term'] = \
        first == 18 and big[18] == 512 and big[20] == -1152 and negative > 0
    num['negative coefficients of the difference up to order 600'] = negative
    pair = sadd(smul(smul(smul(b2, b2), smul(b2, b2)), smul(sadd(smul(a, a), smul(b, b)), sadd(smul(a, a), smul(b, b)))),
                sscale(smul(smul(smul(a2, a2), smul(a2, a2)), smul(smul(a, a), smul(b, b))), 4), -1)
    out['W5 for the pair the same difference is zero (the mean step is exact)'] = pair == {}

    # ---- W5 strong side
    q = Q_FIRST
    out['W5 strong, first piece: right - 1 < 8 q^6 at q = 3/5 (so for every q <= 3/5: right - 1 < 8 q^6 <= left - 1)'] = \
        right_hi(q) - 1 < 8*q**6 and k_lo(q*q) >= 1 + 4*q**6
    steps = int((Q_MEET - Q_FIRST)*GRID)
    grid = [Q_FIRST + F(i, GRID) for i in range(steps + 1)]
    margins = [left_lo(grid[i])/right_hi(grid[i + 1]) for i in range(steps)]
    out['W5 strong, grid 3/5 ... 7/10 in steps of 1/400: left at the lower end > right at the upper end (40 pieces)'] = \
        grid[-1] == Q_MEET and len(margins) == 40 and all(m > 1 for m in margins)
    num['smallest margin on the grid'] = float(min(margins))

    # ---- W5 weak side
    lo, hi = Q_MEET**4, Q_MEET**2
    cuts = [lo + (hi - lo)*F(i, BOXES) for i in range(BOXES + 1)]
    ok, worst, r_top = True, None, F(0)
    for i in range(BOXES):
        r, beta = box(cuts[i], cuts[i + 1])
        bound = 4*(beta[0] + 1)/loss(r[1], beta[1]) - 1                      # beta at every level >= 2 below the box
        ok &= bound >= BETA_MIN
        worst = bound if worst is None else min(worst, bound)
        r_top = max(r_top, r[1])
    out['W5 weak: beta >= 99/10 for every q >= 7/10 (100 boxes of the base, backward walk with bounded loss)'] = ok
    num['smallest certified beta for q >= 7/10'] = float(worst)
    out['W5 weak: beta >= 99/10 gives the inequality: (beta - 1)^2 > 8 beta, and it stays so above'] = \
        (BETA_MIN - 1)**2 > 8*BETA_MIN and 2*(BETA_MIN - 1) > 8
    ok = True
    for _ in range(300):
        r, beta = F(rnd.randint(1, 60), 1000), F(rnd.randint(99, 4000), 10)
        if r*r*beta >= 1:
            continue
        rho = r*beta
        rho_next_sq = r*(1 + r*r)*(beta - 1)**2/(2*(1 - r*r*beta)**2)        # rho'^2 from W1
        ok &= (rho_next_sq*(1 + rho)**2 > 4*rho) == ((1 + r*r)*(beta - 1)**2*(1 + r*beta)**2 > 8*beta*(1 - r*r*beta)**2)
        ok &= rho_next_sq*(1 + rho)**2 > 4*rho
    out['W5 weak: the inequality in (r, beta) form, and true at 300 random points with beta >= 99/10, r^2 beta < 1'] = ok

    # ---- W6 consequences, numbers, controls
    ok = True
    for _ in range(200):
        x = F(rnd.randint(1, 999), 1000)                                     # x = sqrt(rho)
        rho = x*x
        mean_step = 2*x/(1 + rho)
        eps, eps_next = (1 - rho)/(1 + rho), (1 - mean_step)/(1 + mean_step)
        ok &= mean_step > rho and eps_next == eps**2*((1 + rho)/(1 + x)**2)**2 and eps_next < eps**2
    out['W6 algebra: the mean step is above rho, and takes odd/even below its square'] = ok
    rows, ok = [], True
    for u in (F(1, 8), F(1, 4), F(1, 2), F(1)):
        rho, rho2 = at_strong(u, pi_iv)[2], at_strong(2*u, pi_iv)[2]
        eps = ((1 - rho[1])/(1 + rho[1]), (1 - rho[0])/(1 + rho[0]))
        eps2 = ((1 - rho2[1])/(1 + rho2[1]), (1 - rho2[0])/(1 + rho2[0]))
        ok &= eps2[1] < eps[0]**2 and rho[1] < rho2[0]
        rows.append([str(u), float(rho[0]), float(eps[0]), float(eps2[0])])
    out['W6 at u = 1/8, 1/4, 1/2, 1: rho(q^2) > rho(q) and odd/even at q^2 below its square at q (enclosures)'] = ok
    num['u, rho, odd/even, odd/even at 2u'] = rows

    def torus(q):
        """flip of the torus record (weights m^0): (1 - b)/(a - 1), enclosed"""
        od0 = (part(q, 0, 1), part(q, 0, 1) + rest(q, 0))
        ev0 = (part(q, 0, 0), part(q, 0, 0) + rest(q, 0))
        top_ = (2*(od0[0] - ev0[1]), 2*(od0[1] - ev0[0]))
        bot = (2*(od0[0] + ev0[0]), 2*(od0[1] + ev0[1]))
        return (top_[0]/bot[1], top_[1]/bot[0])
    ok = True
    for q in (F(1, 5), F(1, 3), F(1, 2)):
        t, t2 = torus(q), torus(q*q)
        ok &= t2[1]**2*(1 + t[1])**2 < 4*t[0]
    out['W6 control: for the torus record (weights 1 in place of m^2) the inequality is reversed at q = 1/5, 1/3, 1/2'] = ok
    out['W6 control: a claim 5 per cent stronger is not certified by the grid'] = \
        not all(m > F(105, 100) for m in margins)
    out['W6 control: the wrong step (plus in place of minus) fails as a series'] = \
        sscale(smul(a2, kp2), 4) != sadd(smul(a, kp), smul(b, km))
    out['W6 control: a bound of 12.87 is not certified for q >= 7/10 (the value at q = 7/10 is 12.84)'] = \
        worst < BETA_MIN*F(13, 10)

    # floating numbers for the table of the note (not certificates)
    def fl(u):
        qf = math.exp(-math.pi*u)
        if u >= 0.3:
            kpf = sum(m*m*qf**(m*m) for m in range(1, 60))
            kmf = sum((-1)**(m - 1)*m*m*qf**(m*m) for m in range(1, 60))
            return kmf/kpf
        v = 1/u
        t3 = 1 + 2*sum(math.exp(-math.pi*v*n*n) for n in range(1, 30) if math.pi*v*n*n < 700)
        kpv = sum(n*n*math.exp(-math.pi*v*n*n) for n in range(1, 30) if math.pi*v*n*n < 700)
        t2 = 2*sum(math.exp(-math.pi*v*(j + .5)**2) for j in range(30) if math.pi*v*(j + .5)**2 < 700)
        mm = sum((j + .5)**2*math.exp(-math.pi*v*(j + .5)**2) for j in range(30) if math.pi*v*(j + .5)**2 < 700)
        return (mm/u - t2/(4*math.pi))/(t3/(4*math.pi) - kpv/u)
    num['floating: u, rho\'(1+rho)/(2 sqrt rho)'] = \
        [[u, round(fl(2*u)*(1 + fl(u))/(2*math.sqrt(fl(u))), 6)] for u in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)]

    def pair_flip(t):                                                        # (b/a)^2 of the pair at Exp(-pi t)
        if t >= 0.5:
            qf = math.exp(-math.pi*t)
            a_ = 1 + 2*sum(qf**(n*n) for n in range(1, 40))
            b_ = 1 + 2*sum((-1)**n*qf**(n*n) for n in range(1, 40))
            return (b_/a_)**2
        v = 1/t
        t3 = 1 + 2*sum(math.exp(-math.pi*v*n*n) for n in range(1, 30) if math.pi*v*n*n < 700)
        t2 = 2*sum(math.exp(-math.pi*v*(j + .5)**2) for j in range(30) if math.pi*v*(j + .5)**2 < 700)
        return (t2/t3)**2

    def same_flip(u):                                                        # U with pair_flip(U) = rho(u)
        target, lo_, hi_ = fl(u), u, 4*u
        for _ in range(100):
            mid = (lo_ + hi_)/2
            lo_, hi_ = (mid, hi_) if pair_flip(mid) < target else (lo_, mid)
        return lo_
    num['floating: u, U/u (the pair with the same flip sits at Exp(-pi U))'] = \
        [[u, round(same_flip(u)/u, 4)] for u in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5)]
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out.items():
        print('PASS' if v else 'FAIL', k)
    for k, v in num.items():
        print(k, v)
    with open(os.path.join(HERE, 'TW1_RESULT.json'), 'w') as fh:
        json.dump({'checks': {k: bool(v) for k, v in out.items()}, 'numbers': num}, fh, indent=1)
        fh.write('\n')
