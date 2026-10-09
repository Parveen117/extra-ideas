"""SL1 -- the seam law from the walk: the mirror t -> 1/t of the pair record is proved from the diagonal step,
the two means, and one count at the weak end.

Record (AG1):  q = Exp(-pi u),
    a = sum q^(n^2),  b = sum (-1)^n q^(n^2),  c = sum q^((n+1/2)^2)   (n whole) ;   S = a^2, F = b^2, L = c^2.

M1  three lattice identities of one diagonal step (q -> q^2), each by the substitution (x + y, x - y):
        a^2 = a'^2 + c'^2 ,   b^2 = a'^2 - c'^2 ,   c^2 = 2 a' c' .
    So S = S' + L', F = S' - L', L^2 = 4 S' L', and S^2 = F^2 + L^2 follows.
M2  two means.  One step of the mean on (S, F) is the pair at 2u:  M(S, F) = lim S(2^n u) = 1.
    One step of the mean on (S, L) is half the pair at u/2:  2^-n L(u/2^n) <= M(S, L)(u) <= 2^-n S(u/2^n).
M3  the weak-end count (the one input: the area under Exp(-pi x^2) is 1):
        1 - sqrt t <= sqrt t a(t) <= 1 + sqrt t ,   1 - sqrt t <= sqrt t c(t) <= 1 + 2 sqrt t ,
    so with t = u/2^n:  (1 - sqrt t)^2 <= u M(S, L)(u) <= (1 + sqrt t)^2 for every n:   u M(S, L)(u) = 1.
M4  a right triangle is fixed by its two means: x^2 = y^2 + z^2, M(x, y) = 1, M(x, z) = w has one solution.
    (S, F, L) at 1/u and (uS, uL, uF) at u are both solutions for w = u.  Hence
        S(1/u) = u S(u) ,   F(1/u) = u L(u) ,   L(1/u) = u F(u) .
M5  controls and what follows.

Exact integer series in r = q^(1/4) (order 1600), lattice enumerations, directed enclosures (tools/pm1, AG1).
Stdlib only.
"""
from fractions import Fraction as F
from math import isqrt
import importlib.util
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(folder, name):
    path = os.path.join(HERE, '..', folder, name + '.py')
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


hl1 = load('hl1', 'hl1_helical_walk')
ag1, pm1 = hl1.ag1, hl1.pm1
iv, iadd, isub, imul, iscale, idiv = ag1.iv, ag1.iadd, ag1.isub, ag1.imul, ag1.iscale, ag1.idiv

ORDER = 1600                      # in r = q^(1/4)


def series(expo, sign=lambda n: 1):
    """sum over whole n of sign(n) r^expo(n), as {exponent: integer}"""
    out, n = {}, 0
    while True:
        done = True
        for m in ((n, -n) if n else (0,)):
            e = expo(m)
            if e <= ORDER:
                out[e] = out.get(e, 0) + sign(m)
                done = False
        if done and n > 2:
            return {k: v for k, v in out.items() if v}
        n += 1


def smul(x, y):
    out = {}
    for i, p in x.items():
        for j, q in y.items():
            if i + j <= ORDER:
                out[i + j] = out.get(i + j, 0) + p*q
    return {k: v for k, v in out.items() if v}


def sadd(x, y, c=1):
    out = dict(x)
    for k, v in y.items():
        out[k] = out.get(k, 0) + c*v
    return {k: v for k, v in out.items() if v}


def sscale(x, c):
    return {k: c*v for k, v in x.items()}


def pair(e, e4):
    """enclosures of S, F, L from enclosures e of q and e4 of q^(1/4)"""
    t3, t4, t2 = ag1.theta3(e), ag1.theta4(e), ag1.theta2(e4)
    return imul(t3, t3), imul(t4, t4), imul(t2, t2)


def run():
    out, num = {}, {}
    a = series(lambda n: 4*n*n)
    b = series(lambda n: 4*n*n, lambda n: -1 if n % 2 else 1)
    c = {(2*n + 1)**2: 2 for n in range(0, 40) if (2*n + 1)**2 <= ORDER}       # n and -n-1 give the same power
    a1 = series(lambda n: 8*n*n)                                 # a at q^2
    c1 = {2*(2*n + 1)**2: 2 for n in range(0, 40) if 2*(2*n + 1)**2 <= ORDER}

    # ---- M1: the three identities of one step
    out['M1 a^2 = a\'^2 + c\'^2 , b^2 = a\'^2 - c\'^2 , c^2 = 2 a\' c\'  (series, order 400 in q)'] = \
        smul(a, a) == sadd(smul(a1, a1), smul(c1, c1)) and smul(b, b) == sadd(smul(a1, a1), smul(c1, c1), -1) and \
        smul(c, c) == sscale(smul(a1, c1), 2)
    box = range(-12, 13)
    even = sorted((m, n) for m in box for n in box if (m + n) % 2 == 0)
    odd = sorted((m, n) for m in box for n in box if (m + n) % 2 == 1)
    ok = True
    for m, n in even:                                            # (m, n) = (x + y, x - y): m^2 + n^2 = 2 (x^2 + y^2)
        x, y = (m + n)//2, (m - n)//2
        ok &= (x + y, x - y) == (m, n) and m*m + n*n == 2*(x*x + y*y)
    for m, n in odd:                                             # (m, n) = (j + k + 1, j - k): 2 (j+1/2)^2 + 2 (k+1/2)^2
        j, k = (m + n - 1)//2, (m - n - 1)//2
        ok &= (j + k + 1, j - k) == (m, n) and 4*(m*m + n*n) == 2*(2*j + 1)**2 + 2*(2*k + 1)**2
    out['M1 lattice: pairs with m + n even are (x + y, x - y), with m + n odd are (j + k + 1, j - k) (625 pairs)'] = \
        ok and len(even) + len(odd) == 625
    ok = True
    for m in range(-10, 10):
        for n in range(-10, 10):                                 # c^2: (m + 1/2)^2 + (n + 1/2)^2 = (x^2 + y^2)/2, x + y odd
            x, y = m + n + 1, m - n
            ok &= (x + y) % 2 == 1 and 2*((2*m + 1)**2 + (2*n + 1)**2) == 4*(x*x + y*y)
    out['M1 lattice: on the half lattice (m+1/2)^2 + (n+1/2)^2 = (x^2 + y^2)/2 with x + y odd, which gives c^2 = 2 a\' c\' (400 points)'] = ok
    s, f, l = smul(a, a), smul(b, b), smul(c, c)
    s1, l1 = smul(a1, a1), smul(c1, c1)
    out['M1 so S = S\' + L\', F = S\' - L\', L^2 = 4 S\' L\', and the triangle S^2 = F^2 + L^2 follows (series)'] = \
        s == sadd(s1, l1) and f == sadd(s1, l1, -1) and smul(l, l) == sscale(smul(s1, l1), 4) and \
        smul(s, s) == sadd(smul(f, f), smul(l, l))

    # ---- M2: the two means
    pi_iv = pm1.pi_interval()
    ok = True
    for x, y in ((F(5), F(3)), (F(7, 2), F(1, 3))):              # one step of the mean, exact: squares only
        am, gm_sq = (x + y)/2, x*y
        ok &= am*am - gm_sq == ((x - y)/2)**2 and am > 0
    out['M2 one step of the mean keeps the order (mean^2 - product = (half difference)^2 >= 0)'] = ok
    u = F(1)
    s0, f0, l0 = ag1.sfl(u, pi_iv)
    m_sf, m_sl = hl1.mean(s0, f0), hl1.mean(s0, l0)
    ok, rows = True, []
    prev = None
    for n in range(0, 7):
        t = u/2**n
        sn, _, ln = ag1.sfl(t, pi_iv)
        low, high = iscale(ln, F(1, 2**n)), iscale(sn, F(1, 2**n))
        ok &= low[0] <= m_sl[1] and m_sl[0] <= high[1]                         # the mean lies between them
        if prev is not None:
            ok &= prev[0][0] <= low[1] and high[0] <= prev[1][1]               # and they close in
        prev = (low, high)
        rows.append([n, float(low[0]), float(high[1])])
    out['M2 at u = 1: 2^-n L(u/2^n) <= M(S, L) <= 2^-n S(u/2^n) for n = 0 .. 6, closing in (enclosures)'] = \
        ok and ag1.width(m_sl) < F(1, 10**100)
    num['n, 2^-n L(1/2^n), 2^-n S(1/2^n)'] = rows
    s2, f2, _ = ag1.sfl(2*u, pi_iv)
    out['M2 one step of the mean on (S, F) is the pair at 2u, and M(S, F) = 1 (enclosures)'] = \
        ag1.overlap(iscale(iadd(s0, f0), F(1, 2)), s2) and ag1.overlap(imul(f2, f2), imul(s0, f0)) and \
        ag1.overlap(m_sf, iv(1)) and ag1.width(m_sf) < F(1, 10**100)

    # ---- M3: the weak-end count
    ok, rows = True, []
    for root in (F(1, 2), F(1, 4), F(1, 8), F(1, 16)):           # sqrt t, t = 1/4 .. 1/256
        t = root*root
        st, _, lt = ag1.sfl(t, pi_iv)
        ra, rc = iscale(ag1.isqrt_iv(st), root), iscale(ag1.isqrt_iv(lt), root)
        ok &= 1 - root <= ra[0] and ra[1] <= 1 + root and 1 - root <= rc[0] and rc[1] <= 1 + 2*root
        ok &= (1 - root)**2 <= iscale(lt, t)[0] and iscale(st, t)[1] <= (1 + root)**2
        rows.append([str(t), float(ra[0] - 1), float(rc[0] - 1)])
    out['M3 at t = 1/4, 1/16, 1/64, 1/256: 1 - sqrt t <= sqrt t a <= 1 + sqrt t and 1 - sqrt t <= sqrt t c <= 1 + 2 sqrt t'] = ok
    num['t, sqrt t a - 1, sqrt t c - 1'] = rows
    ok = True
    for n in range(0, 7):
        root_sq = u/2**n                                         # (1 - sqrt t)^2 <= u M <= (1 + sqrt t)^2 ; use t >= sqrt t^2
        lo_root = F(isqrt(int(root_sq*10**20)), 10**10)
        hi_root = lo_root + F(1, 10**10)
        ok &= (1 - hi_root)**2 <= iscale(m_sl, u)[0] if hi_root < 1 else True
        ok &= iscale(m_sl, u)[1] <= (1 + hi_root)**2
    out['M3 so u M(S, L)(u) lies in [(1 - sqrt t)^2, (1 + sqrt t)^2] for every t = u/2^n: it is 1 (u = 1, n = 0 .. 6)'] = \
        ok and ag1.overlap(iscale(m_sl, u), iv(1))
    ok = True
    for u in (F(1, 4), F(3)):
        su, _, lu = ag1.sfl(u, pi_iv)
        ok &= ag1.overlap(iscale(hl1.mean(su, lu), u), iv(1))
    out['M3 u M(S, L) = 1 at u = 1/4 and 3 as well (enclosures, 100 digits)'] = ok

    # ---- M4: the triangle is fixed by its two means; the seam law
    ok = True
    ks = [F(1, 5), F(2, 5), F(3, 5), F(4, 5)]
    means = [hl1.mean(iv(1), iv(k)) for k in ks]
    ok &= all(means[i][1] < means[i + 1][0] for i in range(3))   # M(1, k) rises with k
    # Phi(k) = M(1, k')/M(1, k): k = 3/5 has k' = 4/5 and k = 4/5 has k' = 3/5
    phi_35, phi_45 = idiv(means[3], means[2]), idiv(means[2], means[3])
    ok &= phi_45[1] < 1 < phi_35[0]
    out['M4 M(1, k) rises with k, so M(1, k\')/M(1, k) falls: a right triangle is fixed by its two means'] = ok
    ok, rows = True, []
    for u in (F(2), F(3, 2), F(5), F(1, 3)):
        su, fu, lu = ag1.sfl(u, pi_iv)
        sv, fv, lv = ag1.sfl(1/u, pi_iv)
        ok &= ag1.overlap(sv, iscale(su, u)) and ag1.overlap(fv, iscale(lu, u)) and ag1.overlap(lv, iscale(fu, u))
        ok &= max(ag1.width(sv), ag1.width(iscale(su, u))) < F(1, 10**100)
        # both triples solve the same problem: M(x, y) = 1, M(x, z) = u
        ok &= ag1.overlap(hl1.mean(sv, fv), iv(1)) and ag1.overlap(hl1.mean(sv, lv), iv(u))
        ok &= ag1.overlap(hl1.mean(iscale(su, u), iscale(lu, u)), iv(1)) and \
            ag1.overlap(hl1.mean(iscale(su, u), iscale(fu, u)), iv(u))
        rows.append([str(u), float(sv[0]), float(fv[0]), float(lv[0])])
    out['M4 the seam law S(1/u) = u S, F(1/u) = u L, L(1/u) = u F at u = 2, 3/2, 5, 1/3; both triples have means 1 and u'] = ok
    num['u, S(1/u), F(1/u), L(1/u)'] = rows
    ok = True
    for root in (F(2), F(3, 2), F(1, 3)):                        # sqrt u rational: one turn
        u = root*root
        su, fu, lu = ag1.sfl(u, pi_iv)
        sv, fv, lv = ag1.sfl(1/u, pi_iv)
        ok &= ag1.overlap(ag1.isqrt_iv(sv), iscale(ag1.isqrt_iv(su), root))
        ok &= ag1.overlap(ag1.isqrt_iv(fv), iscale(ag1.isqrt_iv(lu), root))
    out['M4 for one turn: a(1/u) = sqrt u a(u), b(1/u) = sqrt u c(u) (u = 4, 9/4, 1/9)'] = ok

    # ---- M5: controls
    u = F(1)
    e, e4 = pm1.exp_neg(F(3)*u), pm1.exp_neg(F(3, 4)*u)          # the weight Exp(-3u) in place of Exp(-pi u)
    s3, f3, l3 = pair(e, e4)
    m3 = hl1.mean(s3, l3)
    third = iscale(pi_iv, F(1, 3))
    out['M5 control: with Exp(-3u) in place of Exp(-pi u) the same walk gives u M(S, L) = pi/3, not 1'] = \
        ag1.overlap(m3, third) and m3[0] > 1 and ag1.width(m3) < F(1, 10**100) and ag1.overlap(hl1.mean(s3, f3), iv(1))
    num['with Exp(-3u): M(S, L) at u = 1'] = float(m3[0])
    wrong = hl1.mean(iscale(s0, 2), iscale(f0, 2))
    out['M5 control: a triple with the wrong hypotenuse (2S, 2F) does not have mean 1'] = not ag1.overlap(wrong, iv(1))
    k1p, _ = ag1.k_contents(F(1), pi_iv)
    out['M5 what follows on the seam: 8 pi K+ = a at u = 1 (the first change of sqrt u a(u) = a(1/u))'] = \
        ag1.overlap(imul(iscale(pi_iv, 8), k1p), ag1.isqrt_iv(s0))
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out.items():
        print('PASS' if v else 'FAIL', k)
    for k, v in num.items():
        print(k, v)
    with open(os.path.join(HERE, 'SL1_RESULT.json'), 'w') as fh:
        json.dump({'checks': {k: bool(v) for k, v in out.items()}, 'numbers': num}, fh, indent=1)
        fh.write('\n')
