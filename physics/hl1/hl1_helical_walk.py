"""HL1 -- the walk is a helical orbit: UGD digits of the pair, the record's own clock, and the three-sector
ledger of the turn block.

Sources matched here (each re-derived on this side, none imported):
    Publications emk-ugd-algebra  EMK-1 (K^2 = 1, det = D_par + D_perp), EMK-2 T2 (Cayley law y1 (+) y2,
                                  seam channel multiplicative), UGD-1 (digit = phase exponent, scale index, seam charge)
    Publications emk-recognition-geometry  EMK-G3 (helical lift: n circuits displace by (n alpha, n beta, n q);
                                  visible return with lifted non-return), UGD-G1 (three closure rules; lossy projection)
    RKF F00-G  Theorem 5.1 (Exp(A(y)) = (1+y)/(1-y)), Theorem 7.1 (Log(xy) = Log x + Log y)
    RKF theorum/43  bilateral strands, jet parity

Record (AG1, TW1), q = Exp(-pi u):
    a = sum q^(n^2), b = sum (-1)^n q^(n^2) ; S = a^2, F = b^2 ; L the lost part, L^2 = S^2 - F^2 ; L_n = L at q^(2^n)
    K+ = sum m^2 q^(m^2), K- = sum (-1)^(m-1) m^2 q^(m^2) ; rho = K-/K+ ; r = b/a

H1  TW1's step is the UGD multiplicative law on the block K+ + K- K (K the seam reflection):
        a'K+' + b'K-' K = (1/4)(a - bK)(K+ + K-K) ;  seam determinant channel multiplicative ;
        in the native chart x = (1+y)/(1-y):  x(rho' r') = x(rho)/x(r).
H2  the diagonal step is a deck transformation: every whole cut-complex number is uniquely
        iota^phi (1+iota)^n beta,  beta primary of odd norm :  a digit (phase mod 4, scale index, odd part);
    sheets of the pair:  S = 1 + L_1 + L_2 + ... ,  F = 1 - L_1 + L_2 + ... ;  the lost part is the sheet before 0.
H3  the record's own clock:  M(S', L') = M(S, L)/2  next to  M(S', F') = M(S, F) = 1 ;  u M(S, L) = 1 (agreement).
H4  three sectors in every coefficient:  K-/b = sum 2^n(N) sigma(N_odd) q^N ,
        K+/a = sum (-1)^(N-1) 2^n(N) sigma(N_odd) q^N ;  scale and count free, the mark tied to sheet 0.
H5  what a power of u cannot see: (1 - Lambda) pi/(4u) - 1 is not zero, is below high powers of u,
        and is the first sheet term (8 - 4 pi/u) Exp(-pi/u).
H6  the two strands (forward and backward walk) on the seam u = 1:  8 pi K+ = a ,  1 - Lambda = 4/(pi S^2).
H7  the toron core: the commutator of two unit blocks sees only the three cut components,
        c = 1 - 2 |u x v|^2 ; it is blind to both centre marks ; two turns give TW1's torus record
        (<chi_m> = 1/m) ; three turns: <c12 c23 c31> = 5/72.

Exact integer series (order 400), exact rationals, directed enclosures (tools/pm1, AG1). Stdlib only.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, isqrt
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


tw1 = load('tw1', 'tw1_turn_block_walk')
ag1 = tw1.ag1
pm1 = tw1.pm1
smul, sadd, sscale, dil = tw1.smul, tw1.sadd, tw1.sscale, tw1.dil
iv, iadd, isub, imul, iscale, idiv = ag1.iv, ag1.iadd, ag1.isub, ag1.imul, ag1.iscale, ag1.idiv

ORDER = 400


def mul(a, b):
    return smul(a, b, ORDER)


def dl(a, j):
    return dil(a, j, ORDER)


# ---------------------------------------------------------------- H2: digits of a whole cut-complex number
def times(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])


def is_primary(x):
    """odd norm, x = 1 mod (1 + iota)^3 : real part odd, iota part even, their sum 1 mod 4"""
    return x[0] % 2 == 1 and x[1] % 2 == 0 and (x[0] + x[1]) % 4 == 1


def digit(x):
    """x = iota^phi (1 + iota)^n beta with beta primary; returns (phi, n, beta)"""
    assert x != (0, 0)
    n = 0
    while (x[0] + x[1]) % 2 == 0:
        x = ((x[0] + x[1])//2, (x[1] - x[0])//2)                 # divide by 1 + iota
        n += 1
    found = []
    for phi in range(4):
        unit = [(1, 0), (0, 1), (-1, 0), (0, -1)][phi]
        inv = [(1, 0), (0, -1), (-1, 0), (0, 1)][phi]
        beta = times(x, inv)
        if is_primary(beta):
            found.append((phi, n, beta))
            assert times(unit, beta) == x
    assert len(found) == 1
    return found[0]


def undigit(phi, n, beta):
    x = beta
    for _ in range(n):
        x = times(x, (1, 1))
    return times(x, [(1, 0), (0, 1), (-1, 0), (0, -1)][phi])


# ---------------------------------------------------------------- H3: the mean of two enclosures
def mean(x, y, steps=60):
    """enclosure of the arithmetic-geometric mean of two positive enclosures (the mean rises with both)"""
    def run_one(p, q, low):
        for _ in range(steps):
            m = (p + q)/2
            g = ag1.isqrt_iv((p*q, p*q))[0 if low else 1]
            p, q = (pm1.down(m), g) if low else (pm1.up(m), g)
        return min(p, q) if low else max(p, q)
    return run_one(x[0], y[0], True), run_one(x[1], y[1], False)


# ---------------------------------------------------------------- H4: sectors of a whole number
def split2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k, n


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


# ---------------------------------------------------------------- H7: unit blocks (rational quaternions)
def qmul(p, q):
    a0, a1, a2, a3 = p
    b0, b1, b2, b3 = q
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)


def qconj(p):
    return (p[0], -p[1], -p[2], -p[3])


def unit_block(w):
    """a rational point of the unit sphere from a rational vector w: (1 - |w|^2, 2w)/(1 + |w|^2)"""
    n = sum(x*x for x in w)
    return ((1 - n)/(1 + n),) + tuple(2*x/(1 + n) for x in w)


def cross_sq(u, v):
    c = (u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])
    return sum(x*x for x in c)


def odd_fact(n):
    """(2n - 1)!!"""
    out = 1
    for k in range(1, 2*n, 2):
        out *= k
    return out


def sphere_moment(exps):
    """Haar mean of x0^e0 x1^e1 x2^e2 x3^e3 on the unit sphere of four dimensions (exact)"""
    if any(e % 2 for e in exps):
        return F(0)
    half = [e//2 for e in exps]
    num = 1
    for h in half:
        num *= odd_fact(h)
    den = 1
    for j in range(sum(half)):
        den *= 4 + 2*j
    return F(num, den)


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


def x_poly(i, j, turns):
    """|u_i x u_j|^2 = |u_i|^2 |u_j|^2 - (u_i . u_j)^2 as a polynomial in the vector parts of `turns` blocks"""
    def var(block, comp):
        e = [0]*(3*turns)
        e[3*block + comp] = 1
        return {tuple(e): 1}
    ni, nj, dot = {}, {}, {}
    for c in range(3):
        ni = padd(ni, pmul(var(i, c), var(i, c)))
        nj = padd(nj, pmul(var(j, c), var(j, c)))
        dot = padd(dot, pmul(var(i, c), var(j, c)))
    return padd(pmul(ni, nj), pmul(dot, dot), -1)


def haar_mean(poly, turns):
    total = F(0)
    for e, c in poly.items():
        term = F(c)
        for t in range(turns):
            term *= sphere_moment((0,) + e[3*t:3*t + 3])
        total += term
    return total


def x_moment(k):
    """mean of |u x v|^(2k) over two unit blocks: <s^(2k)>^2 <sin^(2k) of the angle between the axes>"""
    radial = F(2*comb(2*k + 2, k + 1), 4**(k + 1))
    even = 1
    for j in range(2, 2*k + 1, 2):
        even *= j
    return radial*radial*F(even, odd_fact(k + 1))


def cheb_u(m):
    """coefficients of U_m(c), the content with m + 1 marks as a polynomial in c = scalar part"""
    p0, p1 = [F(1)], [F(0), F(2)]
    if m == 0:
        return p0
    for _ in range(m - 1):
        nxt = [F(0)] + [2*x for x in p1]
        for i, x in enumerate(p0):
            nxt[i] -= x
        p0, p1 = p1, nxt
    return p1


def run():
    out, num = {}, {}
    rnd = random.Random(20261009)
    a, b, kp, km, od = tw1.series(ORDER)
    a2, b2, kp2, km2 = dl(a, 2), dl(b, 2), dl(kp, 2), dl(km, 2)
    s, f = mul(a, a), mul(b, b)

    # ---- H1: the step as a product of blocks x + yK, K^2 = 1
    def block(p, q):                                      # (x1 + y1 K)(x2 + y2 K)
        return sadd(mul(p[0], q[0]), mul(p[1], q[1])), sadd(mul(p[0], q[1]), mul(p[1], q[0]))
    nxt = (sscale(mul(a2, kp2), 4), sscale(mul(b2, km2), 4))
    out['H1 4(a\'K+\' + b\'K-\' K) = (a - bK)(K+ + K-K) as blocks with K^2 = 1 (series)'] = \
        block((a, sscale(b, -1)), (kp, km)) == nxt
    det = lambda p: sadd(mul(p[0], p[0]), mul(p[1], p[1]), -1)
    out['H1 the seam determinant channel is multiplicative: 16(S\'K+\'^2 - F\'K-\'^2) = (S - F)(K+^2 - K-^2)'] = \
        det(nxt) == mul(det((a, b)), det((kp, km))) and det((a, b)) == sadd(s, f, -1)
    out['H1 each seam channel steps by itself: 4(a\'K+\' +- b\'K-\') = (a -+ b)(K+ +- K-)'] = \
        sadd(nxt[0], nxt[1]) == mul(sadd(a, b, -1), sadd(kp, km)) and \
        sadd(nxt[0], nxt[1], -1) == mul(sadd(a, b), sadd(kp, km, -1))
    ok = True
    for _ in range(300):
        y1, y2 = (F(rnd.randint(-99, 99), 100) for _ in range(2))
        comp = (y1 + y2)/(1 + y1*y2)                                           # Cayley coordinate of a product block
        chart = lambda y: (1 + y)/(1 - y)
        ok &= chart(comp) == chart(y1)*chart(y2) and abs(comp) < 1
        blk = (1 + y1*y2, y1 + y2)                                             # (1 + y1 K)(1 + y2 K)
        ok &= blk[1]/blk[0] == comp
    out['H1 the Cayley law (y1 + y2)/(1 + y1 y2) is multiplication in the chart (1+y)/(1-y) (EMK-2 T2, F00-G 7.1)'] = ok
    pi_iv = pm1.pi_interval()
    ok, rows = True, []
    for u in (F(1, 4), F(1, 2), F(1)):
        r, _, rho = tw1.at_strong(u, pi_iv)
        r2, _, rho2 = tw1.at_strong(2*u, pi_iv)
        chart = lambda y: idiv(iadd(iv(1), y), isub(iv(1), y))
        ok &= ag1.overlap(imul(chart(imul(rho2, r2)), chart(r)), chart(rho))
        rows.append([str(u), float(chart(r)[0]), float(chart(rho)[0]), float(chart(imul(rho2, r2))[0])])
    out['H1 at u = 1/4, 1/2, 1: x(rho\' r\') x(r) = x(rho) (enclosures)'] = ok
    num['u, x(r), x(rho), x(rho\' r\')'] = rows

    # ---- H2: the deck transformation and the digits
    box = [(x, y) for x in range(-20, 21) for y in range(-20, 21) if (x, y) != (0, 0)]
    digits = {x: digit(x) for x in box}
    out['H2 every whole cut-complex number of the box is uniquely iota^phi (1+iota)^n beta, beta primary (1680 numbers)'] = \
        len(digits) == 1680 and all(undigit(*d) == x for x, d in digits.items()) and \
        all(2**d[1]*(d[2][0]**2 + d[2][1]**2) == x[0]**2 + x[1]**2 for x, d in digits.items())
    ok = True
    for _ in range(400):
        x, y = rnd.choice(box), rnd.choice(box)
        dx, dy, dz = digit(x), digit(y), digit(times(x, y))
        ok &= dz == ((dx[0] + dy[0]) % 4, dx[1] + dy[1], times(dx[2], dy[2]))
    out['H2 under multiplication phases add mod 4, scale indices add, odd parts multiply (400 pairs)'] = ok
    out['H2 one step is times 1 + iota: scale index + 1, norm doubled; two steps are a doubling and a quarter turn'] = \
        all(digit(times(x, (1, 1))) == (d[0], d[1] + 1, d[2]) for x, d in digits.items()) and \
        times((1, 1), (1, 1)) == (0, 2) and digit((2, 0)) == (3, 2, (1, 0))
    p8 = (1, 0)
    powers = []
    for _ in range(8):
        p8 = times(p8, (1, 1))
        powers.append(p8)
    out['H2 eight steps: the phase is back, the scale is not ((1+iota)^8 = 16) ; no earlier power is a whole real'] = \
        p8 == (16, 0) and digit(p8) == (0, 8, (1, 0)) and all(p[1] != 0 or p[0] < 0 for p in powers[:7])
    sheet = {}
    for x in range(-ORDER, ORDER + 1):
        for y in range(-isqrt(max(ORDER - x*x, 0)), isqrt(max(ORDER - x*x, 0)) + 1):
            n = x*x + y*y
            if 0 < n <= ORDER:
                k = split2(n)[0]
                sheet.setdefault(k, {})
                sheet[k][n] = sheet[k].get(n, 0) + 1
    lost1 = sscale(sadd(s, f, -1), F(1, 2))                                  # L_1 = (S - F)/2 (AG1-A4)
    lost1 = {k: int(v) for k, v in lost1.items()}
    lost = {n: dl(lost1, 2**(n - 1)) for n in range(1, 10)}                  # L_n as a series in q, n >= 1
    out['H2 the sheet of scale index n is L_(n+1): contents of norm 2^n x odd (lattice count, order 400)'] = \
        all(sheet.get(n, {}) == lost[n + 1] for n in range(0, 8))
    later = {}
    for n in range(2, 10):
        later = sadd(later, lost[n])
    out['H2 S = 1 + L_1 + L_2 + ... and F = 1 - L_1 + L_2 + ... : the cut reverses sheet 0 only'] = \
        s == sadd(sadd({0: 1}, lost[1]), later) and f == sadd(sadd({0: 1}, lost[1], -1), later)
    half = {}
    for x in range(-30, 30):
        for y in range(-30, 30):
            n2 = (2*x + 1)**2 + (2*y + 1)**2                                 # 4 x norm of (x + 1/2) + (y + 1/2) iota
            if n2//2 <= ORDER:
                half[n2//2] = half.get(n2//2, 0) + 1                         # exponent of p = q^(1/2)
    out['H2 the lost part is the sheet before 0: the half lattice is (1+iota)^(-1) x odd norms, L_0(p) = L_1 at p'] = \
        half == lost1 and all(((x + y) % 2 == 1) for x, y in box if (x*x + y*y) % 2 == 1)
    out['H2 L_0^2 = 4 L_1 S\' (in p = q^(1/2))'] = \
        mul(lost1, lost1) == sscale(mul(dl(lost1, 2), dl(s, 4)), 4)

    # ---- H3: the record's own clock
    ok = True
    for x, y in ((5, 3), (13, 5), (17, 8), (25, 7)):
        root = isqrt(x*x - y*y)
        ok &= root*root == x*x - y*y and F((x + y) + (x - y), 2) == x and (x + y)*(x - y) == root*root
    out['H3 one step of the mean takes (x + y, x - y) to (x, sqrt(x^2 - y^2)) (exact triples)'] = ok
    out['H3 S\' = (S + F)/2 and L\' = (S - F)/2 (series): the next source and lost part are half sum and half difference'] = \
        sscale(dl(s, 2), 2) == sadd(s, f) and sscale(lost[2], 2) == sadd(dl(s, 2), dl(f, 2), -1)
    ok, rows = True, []
    for u in (F(1, 4), F(1), F(3)):
        s1, f1, l1 = ag1.sfl(u, pi_iv)
        s2, f2, l2 = ag1.sfl(2*u, pi_iv)
        m_sf, m_sl, m_sl2 = mean(s1, f1), mean(s1, l1), mean(s2, l2)
        ok &= ag1.overlap(m_sf, iv(1)) and ag1.width(m_sf) < F(1, 10**100)
        ok &= ag1.overlap(iscale(m_sl2, 2), m_sl) and ag1.width(m_sl) < F(1, 10**100)
        ok &= ag1.overlap(iscale(m_sl, u), iv(1))
        rows.append([str(u), float(m_sf[0]), float(m_sl[0]), float(m_sl2[0])])
    out['H3 at u = 1/4, 1, 3: M(S, F) = 1, M(S\', L\') = M(S, L)/2 and u M(S, L) = 1 (enclosures, 100 digits)'] = ok
    num['u, M(S,F), M(S,L), M(S\',L\')'] = rows

    # ---- H4: three sectors in every coefficient
    c_minus = {n: 2**split2(n)[0]*sigma(split2(n)[1]) for n in range(1, ORDER + 1)}
    c_plus = {n: (1 if n % 2 else -1)*c for n, c in c_minus.items()}
    out['H4 K- = b sum 2^n(N) sigma(N_odd) q^N and K+ = a sum (-1)^(N-1) 2^n(N) sigma(N_odd) q^N (order 400)'] = \
        km == mul(b, c_minus) and kp == mul(a, c_plus)
    fourth = lambda x: mul(mul(x, x), mul(x, x))
    lost0_sq = sadd(fourth(a), fourth(b), -1)
    out['H4 the count sector: L_0^2 = 16 sum over odd m of sigma(m) q^m (order 400)'] = \
        lost0_sq == {m: 16*sigma(m) for m in range(1, ORDER + 1, 2)}
    count4 = {}
    for x in product(range(-6, 7), repeat=4):
        n = sum(t*t for t in x)
        if n <= 36:
            count4[n] = count4.get(n, 0) + 1
    out['H4 it is the block of four: 2 x (whole quaternions of odd norm m) = 16 sigma(m), m <= 35 (lattice count)'] = \
        all(2*count4[m] == 16*sigma(m) == lost0_sq[m] for m in range(1, 36, 2))
    ok = all(c_minus[m*n] == c_minus[m]*c_minus[n]
             for m in range(1, 20) for n in range(1, 20) if m*n <= ORDER and tw1.math.gcd(m, n) == 1)
    out['H4 the coefficient is a product over the sectors, and multiplicative in N'] = ok
    sect = {n: ((n - 1) % 2, split2(n)[0], sigma(split2(n)[1])) for n in range(1, 61)}
    free = True
    for i, j in ((1, 2), (2, 1), (0, 2), (2, 0), (0, 1)):   # agree on sector i, differ on sector j
        free &= any(sect[m][i] == sect[n][i] and sect[m][j] != sect[n][j] for m in sect for n in sect)
    tied = all((v[0] == 1) == (v[1] > 0) for v in sect.values())
    out['H4 scale and count are free of each other, and of the mark; but the mark is tied: it is - exactly off sheet 0'] = \
        free and tied and sect[3] == (0, 0, 4) and sect[6] == (1, 1, 4) and sect[12] == (1, 2, 4) and sect[7] == (0, 0, 8)
    drop_scale = {}
    for n in range(0, 9):
        drop_scale = sadd(drop_scale, dl({m: sigma(m) for m in range(1, ORDER + 1, 2)}, 2**n))
    out['H4 a reading that drops the scale sector (2^n -> 1) is another series: it differs from K-/b at q^2 (1 against 2)'] = \
        drop_scale != c_minus and drop_scale[2] == 1 and c_minus[2] == 2 and drop_scale[1] == c_minus[1]

    # ---- H5: what a power of u cannot see
    def lam_of(u, steps=20):
        terms = []
        for n in range(steps + 1):
            s_, f_, _ = ag1.sfl(u*2**n, pi_iv)
            terms.append(iscale(isub(imul(s_, s_), imul(f_, f_)), 2**n))
        total = iv(0)
        for t in terms[1:]:
            total = iadd(total, t)
        total = (total[0], total[1] + terms[-1][1])                          # rest below the last term (TW1-W2)
        return idiv(total, terms[0])
    ok, rows = True, []
    for u, power in ((F(1, 4), 5), (F(1, 8), 9), (F(1, 16), 15), (F(1, 64), 40)):
        lam = lam_of(u)
        delta = isub(imul(isub(iv(1), lam), iscale(pi_iv, 1/(4*u))), iv(1))  # (1 - Lambda) pi/(4u) - 1
        e = pm1.exp_neg_pi(1/u, pi_iv)
        first = imul(isub(iv(8), iscale(pi_iv, 4/u)), e)                     # (8 - 4 pi/u) Exp(-pi/u)
        ratio = idiv(iscale(delta, -1), iscale(first, -1))
        ok &= delta[1] < 0 and -delta[0] < u**power and F(999, 1000) < ratio[0] and ratio[1] < F(1001, 1000)
        rows.append([str(u), float(delta[0]), float(first[0]), power])
    out['H5 at u = 1/4, 1/8, 1/16, 1/64: (1 - Lambda) pi/(4u) - 1 is negative, below u^5, u^9, u^15, u^40, '
        'and within a thousandth of (8 - 4 pi/u) Exp(-pi/u)'] = ok
    num['u, residue, first sheet term, power of u it is below'] = rows

    # ---- H6: the two strands on the seam
    s1, f1, l1 = ag1.sfl(F(1), pi_iv)
    a1 = ag1.isqrt_iv(s1)
    kp1, km1 = ag1.k_contents(F(1), pi_iv)
    lam1 = lam_of(F(1))
    seam = idiv(iv(4), imul(pi_iv, imul(s1, s1)))
    out['H6 on the seam u = 1: F = L, 8 pi K+ = a and 1 - Lambda = 4/(pi S^2) (enclosures, 100 digits)'] = \
        ag1.overlap(f1, l1) and ag1.overlap(imul(iscale(pi_iv, 8), kp1), a1) and \
        ag1.width(imul(iscale(pi_iv, 8), kp1)) < F(1, 10**100) and ag1.overlap(isub(iv(1), lam1), seam)
    num['seam: a, 8 pi K+, 1 - Lambda'] = [float(a1[0]), float(imul(iscale(pi_iv, 8), kp1)[0]), float(1 - lam1[1])]
    off = F(11, 10)
    so, _, _ = ag1.sfl(off, pi_iv)
    kpo, _ = ag1.k_contents(off, pi_iv)
    out['H6 control: off the seam (u = 11/10) 8 pi K+ is not a'] = \
        not ag1.overlap(imul(iscale(pi_iv, 8), kpo), ag1.isqrt_iv(so))

    # ---- H7: the toron core
    ok, parallel = True, True
    for _ in range(200):
        u = unit_block([F(rnd.randint(-9, 9), rnd.randint(1, 9)) for _ in range(3)])
        v = unit_block([F(rnd.randint(-9, 9), rnd.randint(1, 9)) for _ in range(3)])
        comm = qmul(qmul(u, v), qmul(qconj(u), qconj(v)))
        x = cross_sq(u, v)
        ok &= sum(t*t for t in u) == 1 and comm[0] == 1 - 2*x and 0 <= x <= 1
        for su, sv in product((1, -1), repeat=2):
            uu, vv = tuple(su*t for t in u), tuple(sv*t for t in v)
            ok &= qmul(qmul(uu, vv), qmul(qconj(uu), qconj(vv)))[0] == comm[0]
        k = F(rnd.randint(-5, 5), 3)
        w = unit_block([F(rnd.randint(1, 9), 7) for _ in range(3)])
        t = unit_block([k*c for c in (w[1], w[2], w[3])])                    # an axis parallel to w's
        parallel &= qmul(qmul(w, t), qmul(qconj(w), qconj(t)))[0] == 1
    out['H7 the commutator of two unit blocks has scalar part 1 - 2|u x v|^2: only the three cut components enter; '
        'blind to both centre marks; 1 on a common axis (200 exact pairs)'] = ok and parallel
    x12 = x_poly(0, 1, 2)
    direct = [haar_mean(x12, 2), haar_mean(pmul(x12, x12), 2), haar_mean(pmul(pmul(x12, x12), x12), 2)]
    out['H7 moments of |u x v|^2 over two unit blocks: 3/8, 5/24, 35/256 by two routes'] = \
        direct == [x_moment(1), x_moment(2), x_moment(3)] == [F(3, 8), F(5, 24), F(35, 256)]
    ok = True
    for m in range(1, 13):
        coef = cheb_u(m - 1)                                                 # chi_m(c), c = 1 - 2X
        total = F(0)
        for k, ck in enumerate(coef):
            total += ck*sum(comb(k, j)*(-2)**j*x_moment(j) for j in range(k + 1))
        ok &= total == F(1, m)
    out['H7 two turns: the mean of the content with m marks on the commutator is 1/m, m <= 12: TW1\'s torus record'] = ok
    out['H7 one turn at its centre points: the content with m marks reads m and (-1)^(m-1) m (the next jet; TW1\'s K+, K-)'] = \
        all(sum(cheb_u(m - 1)) == m and sum(c*(-1)**i for i, c in enumerate(cheb_u(m - 1))) == (-1)**(m - 1)*m
            for m in range(1, 13))
    x3 = [x_poly(0, 1, 3), x_poly(1, 2, 3), x_poly(2, 0, 3)]
    one = {tuple([0]*9): 1}
    c3 = [padd(one, p, -2) for p in x3]
    triple = haar_mean(pmul(pmul(c3[0], c3[1]), c3[2]), 3)
    pair = haar_mean(pmul(c3[0], c3[1]), 3)
    def sphere2(e):                                           # mean of x^e0 y^e1 z^e2 on the sphere of directions
        if any(t % 2 for t in e):
            return F(0)
        h = [t//2 for t in e]
        return F(odd_fact(h[0])*odd_fact(h[1])*odd_fact(h[2]), odd_fact(sum(h) + 1))
    # second route: radial parts times the directions; the middle axis fixed along z, n1 and n3 free
    ang = F(0)
    for e1 in product(range(0, 5), repeat=3):
        for e3 in product(range(0, 5), repeat=3):
            if sum(e1) > 4 or sum(e3) > 4:
                continue
            # coefficient of n1^e1 n3^e3 in (1 - z1^2)(1 - z3^2)(1 - (n1.n3)^2)
            coef = 0
            for t1, t3, td in product((0, 1), repeat=3):       # choose 1 or the subtracted square in each factor
                sign = (-1)**(t1 + t3 + td)
                base1, base3 = [0, 0, 2*t1], [0, 0, 2*t3]
                if td == 0:
                    if tuple(base1) == e1 and tuple(base3) == e3:
                        coef += sign
                else:
                    for i in range(3):
                        for j in range(3):
                            m1, m3 = base1[:], base3[:]
                            m1[i] += 1; m1[j] += 1; m3[i] += 1; m3[j] += 1
                            if tuple(m1) == e1 and tuple(m3) == e3:
                                coef += sign
            if coef:
                ang += coef*sphere2(e1)*sphere2(e3)
    radial4 = F(2*comb(6, 3), 4**3)                           # mean of s^4 = 5/8
    out['H7 three turns: <c12> = 1/4, <c12 c23> = 1/8 (= mean of c^4, CT1), <c12 c23 c31> = 5/72 ; '
        'the triple by two routes (directions: 64/225)'] = \
        haar_mean(c3[0], 3) == F(1, 4) and pair == F(1, 8) == F(comb(4, 2), 3*16) and triple == F(5, 72) and \
        ang == F(64, 225) and haar_mean(pmul(pmul(x3[0], x3[1]), x3[2]), 3) == radial4**3*ang
    num['three turns: <c12 c23 c31>'] = str(triple)
    out['H7 the closed triple is above the open pair times a face (1/32) and above three free faces (1/64)'] = \
        triple > pair*F(1, 4) > F(1, 64)
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out.items():
        print('PASS' if v else 'FAIL', k)
    for k, v in num.items():
        print(k, v)
    with open(os.path.join(HERE, 'HL1_RESULT.json'), 'w') as fh:
        json.dump({'checks': {k: bool(v) for k, v in out.items()}, 'numbers': num}, fh, indent=1)
        fh.write('\n')
