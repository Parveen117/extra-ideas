"""GC1 -- the gap of the core as a count: at most one reading of the core lies under a certified line.

The core (TC1, CR1): C a free 3 x d matrix with columns c_1 .. c_d (the turns), X the sum over pairs of
|c_i cross c_j|^2, h = -(1/2) Laplacian + X.  Readings unchanged by turning all the cuts together (the row
condition) are the physical ones.  CR1-R4 gave the comparison
        h  >=  k_1 + .. + k_d ,      k = -(1/4) Laplacian + beta |c| ,   beta^2 = (d - 1)/2 ,
one copy of k for each turn.  This stage counts with it.

G1  one turn.  With r = l s, l^3 = 1/(4 beta), k is  unit x ( -u'' + [n(n+1)/s^2] u + s u )  on u = s f,
    unit = ((d - 1)/8)^(1/3).  The solution v of v'' = (s - lam) v, v(0) = 0, v'(0) = 1 is an exact power series.
        lam = 2.3381 :  v > 0 on the whole half line            =>  first level  >= 2.3381
        lam = 4.0879 :  v > 0, one node, v < 0 after it         =>  second level >= 4.0879
    (a positive v with -v'' + (s - lam) v = 0 gives  form >= lam  on every interval where it keeps one sign).
    Readings of one turn with zero average over every sphere: local rate of s^2 exp(-k s^(3/2)), k^2 = 2/9,
        level >= (9/4) 3^(1/3) = 3.2451 .
G2  d turns under the row condition.  Split every turn into its sphere average and the rest.  A row-unchanged
    reading has no part with exactly one turn outside its sphere average.  Hence at most ONE reading of the
    comparison lies under
        z = unit x min( (d - 1) e0 + e1 , (d - 2) e0 + 2 p0 ) :      5.5210 (three turns) ,  8.0060 (four) .
G3  the gap.  The lowest rate is under CR1's ladder: 5.1868, 8.0030.  So
        second rate - lowest rate  >=  0.334 (three turns) ,  0.003 (four turns) .
G4  the lowest rate from both sides.  With no second reading under z, for any reading psi
        lowest rate  >=  eta - (mean of h^2 - eta^2)/(z - eta) ,      eta = mean of h in psi :
        three turns  [5.1865, 5.1868] ,      four turns  [7.994, 8.0030] .
G5  with CM1's turning readings (7.5787, 10.9738) above:  gap in [0.334, 2.393] and [0.003, 2.98].
G6  the walk in the number of turns.  The core of n turns is the sum, over its n sub-cores of n - 1 turns, of
    a T + b X with a = 1/(n - 1), b = 1/(n - 2); a dilation makes each (a^2 b)^(1/3) x the core of n - 1 turns.
    Hence, for the lowest rate E(d),
        rho_d = E(d) / ( d (d - 1)^(1/3) )   never falls as d grows,   and   rho_d^3 <= 729/256  (Gaussian reading).
        rho_3 in [1.3721, 1.3723] ,  rho_4 in [1.3857, 1.3873] ,  every d >= 4:  1.3857 <= rho_d <= 1.4175 .

Exact integers and rationals; the trial vector is found in floats and then used as an exact rational vector.
Stdlib only.
"""
from fractions import Fraction as F
from math import factorial, sqrt
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


cm1 = load('cm1', 'cm1_core_sector_memory')
cr1 = cm1.cr

E0 = F(23381, 10000)                      # certified floor of the first level of -u'' + s u on the half line
E1 = F(40879, 10000)                      # of the second
P0 = F(3245, 1000)                        # of readings with zero sphere average: (9/4) 3^(1/3) = 3.24507..
UNIT = {3: F(62996, 100000), 4: F(72112, 100000)}     # floors of ((d - 1)/8)^(1/3)
STEP = 256                                # grid of the half line: 1/256
REACH = {1: 8, 2: 9}                      # where each walk on the half line ends
TOP = 10                                  # degree of the trial reading (CR1's ladder)
SHIFT = {3: 5.18, 4: 7.99}                # float shifts for the inverse iteration
DIGITS = 12


# ---------------------------------------------------------------- one turn: the solution on the half line
def series(lam, reach):
    """Integers A_n with c_n = A_n/(n! q^n), lam = p/q, for v'' = (s - lam) v, v(0) = 0, v'(0) = 1; and a
    bound for what the series beyond the last kept term can add to v and to v' on [0, reach]."""
    p, q = lam.numerator, lam.denominator
    a = [0, q, 0]
    reach = F(reach)

    def term(n):
        return abs(F(a[n], factorial(n)*q**n))*reach**n

    n = 1
    while True:
        a.append(n*q**3*a[n - 1] - p*q*a[n])          # A_(n+2)
        theta = (reach**3 + lam*reach**2)/((n + 3)*(n + 2))
        if theta <= F(1, 8) and len(a) > n + 2 and n >= 3:
            m0 = max(term(n), term(n + 1), term(n + 2))
            tail = 3*m0*theta/(1 - theta)
            if tail < F(1, 10**(DIGITS + 8)):
                slope = 3*m0*((n + 2)*theta/(1 - theta) + 3*theta/(1 - theta)**2)/reach
                return a[:n + 3], tail, slope
        n += 1


def residual_degree(a, lam):
    """Lowest degree left in p'' - (s - lam) p for the kept polynomial p: it must be the last two only."""
    q = lam.denominator
    c = [F(x, factorial(n)*q**n) for n, x in enumerate(a)]
    top = len(c) - 1
    out = {}
    for n in range(2, top + 1):
        out[n - 2] = out.get(n - 2, 0) + n*(n - 1)*c[n]
    for n in range(top + 1):
        out[n + 1] = out.get(n + 1, 0) - c[n]
        out[n] = out.get(n, 0) + lam*c[n]
    left = [n for n, x in out.items() if x]
    return min(left), top


def grid_values(lam, reach, step=STEP):
    """Enclosures (lo, hi) of v and of v' at s = k/step, k = 0 .. reach*step."""
    a, tail, slope = series(lam, reach)
    q = lam.denominator
    top = len(a) - 1
    big = q*step
    b = [a[n]*(factorial(top)//factorial(n))*big**(top - n) for n in range(top + 1)]
    den = factorial(top)*big**top
    scale = 10**DIGITS
    pad_v = F(int(tail*scale) + 1, scale)
    pad_d = F(int(slope*scale) + 1, scale)
    values = []
    for k in range(reach*step + 1):
        s = 0
        for n in range(top, -1, -1):
            s = s*k + b[n]
        t = 0
        for n in range(top, 0, -1):
            t = t*k + n*b[n]
        lo = (s*scale)//den
        dlo = (t*step*scale)//den
        values.append(((F(lo, scale) - pad_v, F(lo + 1, scale) + pad_v),
                       (F(dlo, scale) - pad_d, F(dlo + 1, scale) + pad_d)))
    return values


def steps(lam, reach, step=STEP):
    """For each grid step: enclosure of v and of v' on the whole step, from v'' = (s - lam) v."""
    values = grid_values(lam, reach, step)
    h = F(1, step)
    out = []
    for k in range(len(values) - 1):
        (vl, vu), (dl, du) = values[k]
        bend = max(abs(F(k, step) - lam), abs(F(k + 1, step) - lam))
        if not bend*h*h < 1:
            raise ValueError('grid too coarse for the step bound')
        size = (max(abs(vl), abs(vu)) + max(abs(dl), abs(du))*h)/(1 - bend*h*h/2)
        err = bend*size*h*h/2
        out.append(dict(v=(min(vl, vl + dl*h) - err, max(vu, vu + du*h) + err),
                        d=(dl - bend*size*h, du + bend*size*h), left=(vl, vu)))
    return out, values[-1]


def level_one(lam, reach=REACH[1], step=STEP):
    """True when v > 0 on the whole half line is certified: then the first level is at least lam."""
    st, ((vl, vu), (dl, du)) = steps(lam, reach, step)
    if not st[0]['d'][0] > 0:
        return False
    if any(not x['v'][0] > 0 for x in st[1:]):
        return False
    return reach >= lam and vl > 0 and dl > 0


def level_two(lam, reach=REACH[2], step=STEP):
    """(True, node) when v > 0, then one node, then v < 0 for ever is certified: second level >= lam."""
    st, ((vl, vu), (dl, du)) = steps(lam, reach, step)
    if not st[0]['d'][0] > 0:
        return False, None
    state, first, last = 'pos', None, None
    for k in range(1, len(st)):
        x = st[k]
        if state == 'pos':
            if x['v'][0] > 0:
                continue
            if x['d'][1] < 0 and x['left'][0] > 0:
                state, first = 'fall', k
            else:
                return False, None
        if state == 'fall':
            if x['v'][1] < 0:
                state, last = 'neg', k
            elif not x['d'][1] < 0:
                return False, None
            continue
        if state == 'neg' and not x['v'][1] < 0:
            return False, None
    ok = state == 'neg' and reach >= lam and vu < 0 and du < 0
    return ok, (F(first, step), F(last, step)) if ok else None


# ---------------------------------------------------------------- d turns: the line under which one reading lies
def line(d, rows=True):
    """z of G2.  With rows=False the row condition is dropped (a single turn may leave its sphere average)."""
    cases = [(d - 1)*E0 + E1, (d - 2)*E0 + 2*P0]
    if not rows:
        cases.append((d - 1)*E0 + P0)
    return UNIT[d]*min(cases)


# ---------------------------------------------------------------- the trial reading and its two means
def solve_float(a, b):
    n = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        piv = max(range(c, n), key=lambda i: abs(m[i][c]))
        m[c], m[piv] = m[piv], m[c]
        for i in range(c + 1, n):
            f = m[i][c]/m[c][c]
            if f:
                for j in range(c, n + 1):
                    m[i][j] -= f*m[c][j]
    x = [0.0]*n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][j]*x[j] for j in range(i + 1, n)))/m[i][i]
    return x


def near_lowest(s, h, n, shift, rounds=6):
    """A rational vector near the lowest direction of the ladder (floats; only its exact means are used)."""
    sc = [1/sqrt(float(s[i][i])) for i in range(n)]
    sf = [[float(s[i][j])*sc[i]*sc[j] for j in range(n)] for i in range(n)]
    af = [[float(h[i][j])*sc[i]*sc[j] - shift*sf[i][j] for j in range(n)] for i in range(n)]
    x = [1.0] + [0.0]*(n - 1)
    for _ in range(rounds):
        y = solve_float(af, [sum(sf[i][j]*x[j] for j in range(n)) for i in range(n)])
        top = max(abs(t) for t in y)
        x = [t/top for t in y]
    return [F(x[i]*sc[i]) for i in range(n)]


def quad(m, x, y=None):
    y = x if y is None else y
    return sum(x[i]*sum(m[i][j]*y[j] for j in range(len(y)) if y[j]) for i in range(len(x)) if x[i])


def trial(d, top=TOP):
    """eta = mean of h, and the mean of h^2, in one exact rational reading of CR1's ladder of degree top."""
    w = cr1.RATE[d]
    basis, s, h = cr1.ladder(d, w, top + 2)
    n = len(cm1.inv_basis(top))
    assert basis[:n] == cm1.inv_basis(top)
    index = {e: i for i, e in enumerate(basis)}
    x = near_lowest(s, h, n, SHIFT[d])
    y = [F(0)]*len(basis)                               # h applied to the reading, in the ladder two higher
    for j in range(n):
        for e, v in cm1.scalar_h({basis[j]: F(1)}, d, w).items():
            y[index[e]] += v*x[j]
    sxx = quad(s, x)
    hxx = quad(h, x)
    txx = quad(s, y)
    return dict(eta=hxx/sxx, square=txx/sxx, same=quad(s, x, y) == hxx, norm=sxx > 0, size=n)


def floor_of_lowest(eta, square, rho):
    """Lowest rate >= eta - (square - eta^2)/(rho - eta), when no second reading lies under rho > eta."""
    assert rho > eta
    return eta - (square - eta*eta)/(rho - eta)


def turning_upper(d):
    with open(os.path.join(HERE, '..', 'cm1', 'CM1_RESULT.json')) as f:
        rows = json.load(f)['numbers']['turning_sector_d%d' % d]
    return F(rows[-1]['ritz_upper'])


def cut(x, places=4, up=False):
    t = x*10**places
    n = t.numerator//t.denominator
    if up and n != t:
        n += 1
    return F(n, 10**places)


def root3(x, places=4, up=False):
    """Cube root of a positive rational, cut to the given places from below (or from above)."""
    t = 10**places
    n = int(round(float(x)**(1/3)*t))
    while F(n, t)**3 > x:
        n -= 1
    while F(n + 1, t)**3 <= x:
        n += 1
    return F(n + 1, t) if up and F(n, t)**3 != x else F(n, t)


# ---------------------------------------------------------------- the certificate
def run():
    checks, num = {}, {}

    # G1: one turn
    for name, lam, reach in (('first', E0, REACH[1]), ('second', E1, REACH[2])):
        a, tail, slope = series(lam, reach)
        low, top = residual_degree(a, lam)
        checks['G1 %s: the kept polynomial solves v\'\' = (s - lam) v up to its last two degrees' % name] = \
            low >= top - 1
        checks['G1 %s: what the series can still add is under 10^-20' % name] = \
            tail < F(1, 10**20) and slope < F(1, 10**17)
        num['terms kept, %s walk' % name] = top
    checks['G1 v > 0 on the whole half line at lam = 2.3381'] = level_one(E0)
    checks['G1 control: at lam = 2.3382 the same certificate is refused'] = not level_one(E0 + F(1, 10000))
    ok, node = level_two(E1)
    checks['G1 v > 0, one node, v < 0 for ever at lam = 4.0879'] = ok
    checks['G1 control: at lam = 4.0880 the same certificate is refused'] = not level_two(E1 + F(1, 10000))[0]
    num['node of the second walk'] = [float(node[0]), float(node[1])] if ok else None
    # zero sphere average: A/t + B t^2 >= 3 (A^2 B/4)^(1/3), A^2 = 729 k^2/16, B = 1 - 9k^2/4, k^2 = 2/9
    k2 = F(2, 9)
    a2, bb = F(729, 16)*k2, 1 - F(9, 4)*k2
    checks['G1 zero sphere average: p0^3 <= 27 A^2 B/4 = 2187/64'] = \
        27*a2*bb/4 == F(2187, 64) and P0**3 <= F(2187, 64)
    checks['G1 order of the floors: e0 <= e1, e0 <= p0'] = E0 <= E1 and E0 <= P0
    num['floors of one turn (unit 1)'] = dict(first=float(E0), second=float(E1), zero_average=float(P0))

    # G2 .. G5
    window = {}
    for d in (3, 4):
        checks['G2 unit^3 <= (d - 1)/8, d = %d' % d] = UNIT[d]**3 <= F(d - 1, 8) < (UNIT[d] + F(1, 10**5))**3
        z = line(d)
        t = trial(d)
        eta, square = t['eta'], t['square']
        checks['G4 the reading is exact: <psi, h psi> two ways; norm positive, d = %d' % d] = \
            t['same'] and t['norm']
        ladder_top = {3: F(51868, 10000), 4: F(80030, 10000)}[d]
        checks['G3 mean of h in the reading under CR1\'s bound, d = %d' % d] = eta < ladder_top
        checks['G3 the line is above the reading: z > eta, d = %d' % d] = z > eta
        gap_low = z - eta
        checks['G3 gap floor, d = %d' % d] = gap_low > {3: F(334, 1000), 4: F(3, 1000)}[d]
        low = floor_of_lowest(eta, square, z)
        checks['G4 spread of the reading is positive and small, d = %d' % d] = \
            0 < square - eta*eta < F(1, 10**4)
        checks['G4 floor of the lowest rate, d = %d' % d] = low > {3: F(51865, 10000), 4: F(7994, 1000)}[d]
        checks['G4 above the earlier floors (CR1-R4; d copies of one turn), d = %d' % d] = \
            low > d*E0*UNIT[d] > {3: F(414, 100), 4: F(632, 100)}[d]
        free = line(d, rows=False)
        checks['G2 control: without the row condition the count stops under the lowest rate, d = %d' % d] = \
            free < low
        up = turning_upper(d)
        checks['G5 the turning reading of CM1 is above the line, d = %d' % d] = up > z
        gap_high = up - low
        checks['G5 gap window, d = %d' % d] = gap_low < gap_high < {3: F(2393, 1000), 4: F(298, 100)}[d]
        window[d] = (low, eta)
        num['d = %d' % d] = dict(
            line=float(cut(z)), line_without_rows=float(cut(free)),
            lowest_rate=[float(cut(low)), float(cut(eta, up=True))],
            second_rate=[float(cut(z)), float(cut(up, up=True))],
            gap=[float(cut(gap_low)), float(cut(gap_high, up=True))],
            spread=float(square - eta*eta), readings=t['size'])

    # G6: the walk in the number of turns
    for d in (3, 4, 5):
        mean = cr1.free_means(d, F(1))
        checks['G6 Gaussian mean of X is 3d(d - 1)/4 at w = 1, d = %d' % d] = \
            mean({(0, 1, 0): F(1)}) == F(3*d*(d - 1), 4)
        # value(w) = 3dw/4 + 3d(d-1)/(4w^2) = u + u + t >= 3 (u u t)^(1/3), u = 3dw/8 ; 27 u u t = 729 d^3 (d-1)/256
        w = root3(F(2*(d - 1)), 6)
        u, t = 3*d*w/8, F(3*d*(d - 1), 4)/(w*w)
        best = F(729*d**3*(d - 1), 256)
        checks['G6 Gaussian reading: value^3 >= 729 d^3 (d - 1)/256, met within 10^-9 at w^3 = 2(d - 1), d = %d' % d] = \
            27*u*u*t == best and 0 <= (2*u + t)**3 - best < F(1, 10**9)
    checks['G6 one step: n^3 a^2 b = [n (n-1)^(1/3)]^3 / [(n-1)(n-2)^(1/3)]^3 , n = 3 .. 12'] = all(
        n**3*F(1, (n - 1)**2*(n - 2)) == F(n**3*(n - 1), (n - 1)**3*(n - 2)) for n in range(3, 13))
    (low3, up3), (low4, up4) = window[3], window[4]
    checks['G6 the certified windows stand in the order of the walk: rho_3 < rho_4'] = up3**3/54 < low4**3/192
    checks['G6 the step from three turns gives a floor for four under its ceiling, and under the floor of G4'] = \
        F(192, 54)*low3**3 < low4**3 < up4**3
    checks['G6 every rho_d, d >= 4, between rho_4 and the Gaussian value'] = low4**3/192 < F(729, 256)
    num['rho'] = {'d = 3': [float(root3(low3**3/54)), float(root3(up3**3/54, up=True))],
                  'd = 4': [float(root3(low4**3/192)), float(root3(up4**3/192, up=True))],
                  'd >= 4': [float(root3(low4**3/192)), float(root3(F(729, 256), up=True))]}
    num['five turns, lowest rate'] = [float(root3(500*low4**3/192)), float(root3(F(729*500, 256), up=True))]
    num['four turns from three by one step'] = float(root3(F(192, 54)*low3**3))
    num['two turns, lowest rate at most'] = float(root3(8*up3**3/54, up=True))

    out = dict(stage='GC1', checks=len(checks), passed=sum(1 for v in checks.values() if v),
               all_pass=all(checks.values()), detail=checks, numbers=num)
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out['detail'].items():
        print('PASS' if v else 'FAIL', k)
    print(json.dumps(num, indent=1))
    print(out['passed'], '/', out['checks'])
    with open(os.path.join(HERE, 'GC1_RESULT.json'), 'w') as f:
        json.dump(out, f, indent=1)
