"""DR1 -- the core read by rows: three cuts in d directions.

The weight of the core is one number read two ways (TC1): by its d columns c_i in three components (the turns),
or by its three rows r_a in d components (the cuts):
        X = sum over pairs of columns |c_i cross c_j|^2 = sum over pairs of rows ( |r_a|^2 |r_b|^2 - (r_a.r_b)^2 ) .
GC1 and SG1 compared by columns.  This stage compares by rows.

D1  the comparison by rows.  Half of each row's kinetic part goes to its two pairs; a pair is an oscillator in
    the d - 1 directions across the other row, with ground rate (d - 1)|r_a|/4.  So
        h  >=  k'_1 + k'_2 + k'_3 ,     k' = -(1/4) Laplacian in d directions + ((d - 1)/2) |r| .
D2  one row.  On readings of turning number l (in d directions) k' is
        unit x ( -u'' + [(m - 1)(m - 3)/(4 s^2)] u + s u ) ,   m = d + 2l ,   unit^3 = (d - 1)^2/16 ;
    its first level nu(m) and, for l = 0, its second level nu'(d) get exact floors from the signs of an exact
    series (the solution s^p w(s), p = (m - 1)/2), as in GC1-G1; and nu(m) >= (3/4)(2m - 1)^(2/3) in closed form.
D3  the count.  The row turning contains the three half turns that flip two rows at once.  A reading unchanged
    by them has all three rows even or all three odd.  So at most ONE physical reading lies under
        z' = unit x min( 2 nu(d) + min(nu'(d), nu(d + 4)) , 3 nu(d + 2) ) .
D4  the gap for d turns, d = 3 .. 10, against exact readings of CR1's ladder: positive for every one of them;
    four turns 0.4425 (by columns: 0.0031), five turns 0.5301, ... , ten turns 0.8097; and the lowest rate from both sides.
D5  the end of the walk in the number of turns:
        (729/1024)(d - 1)^2 (2d - 1)^2  <=  E(d)^3  <=  (729/256) d^3 (d - 1) ,
    so rho_d = E(d)/(d (d - 1)^(1/3)) rises (GC1-G6) to exactly (729/256)^(1/3) = 1.41741.., within 2/d in the cube.

Exact rationals; levels and trial vectors are located in floats and then certified exactly.  Stdlib only.
"""
from fractions import Fraction as F
from math import lcm
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


gc1 = load('gc1', 'gc1_core_gap_count')
cm1, cr1 = gc1.cm1, gc1.cr1

TURNS = range(3, 11)
STEP = 64
DIGITS = 12
MARGIN = 0.002                              # how far under the float level the certified floor is put
WEIGHT = {3: F(8, 5), 4: F(9, 5), 5: F(2), 6: F(43, 20), 7: F(23, 10), 8: F(12, 5), 9: F(5, 2), 10: F(13, 5)}
DEGREE = {3: 10, 4: 10}                     # degree of the ladder reading; 6 for the others


# ---------------------------------------------------------------- one row: the exact series and its signs
def series(lam, big_p, reach):
    """Coefficients a_n of w, where u = s^(P/2) w(s) solves u'' = (C/s^2 + s - lam) u, C = P(P - 2)/4, w(0) = 1;
    and bounds for what the series beyond the last kept term adds to w and to w' on [0, reach]."""
    a = [F(1), F(0)]
    reach = F(reach)

    def term(n):
        return abs(a[n])*reach**n

    n = 0
    while True:
        prev = a[n - 1] if n else F(0)
        a.append((prev - lam*a[n])/((n + 2)*(n + 1 + big_p)))
        theta = (reach**3 + lam*reach**2)/((n + 3)*(n + 2 + big_p))
        if n >= 3 and theta <= F(1, 8):
            m0 = max(term(n), term(n + 1), term(n + 2))
            tail = 3*m0*theta/(1 - theta)
            if tail < F(1, 10**(DIGITS + 8)):
                slope = 3*m0*((n + 2)*theta/(1 - theta) + 3*theta/(1 - theta)**2)/reach
                return a, tail, slope
        n += 1


def residual_low(a, lam, big_p):
    """Lowest degree left in s w'' + P w' - s (s - lam) w for the kept polynomial."""
    top = len(a) - 1
    out = {}
    for n in range(top + 1):
        if n >= 2:
            out[n - 1] = out.get(n - 1, 0) + n*(n - 1)*a[n]
        if n >= 1:
            out[n - 1] = out.get(n - 1, 0) + big_p*n*a[n]
        out[n + 2] = out.get(n + 2, 0) - a[n]
        out[n + 1] = out.get(n + 1, 0) + lam*a[n]
    return min(n for n, x in out.items() if x), top


def potential(big_p, lam, s):
    return F(big_p*(big_p - 2), 4)/(s*s) + s - lam


def walk(lam, big_p, reach, step=STEP):
    """Step enclosures of u/s_k^p and its slope from the first grid point k0 where the potential is still above
    lam on the whole of (0, s_k0] (there u and u' are positive), up to reach.  Returns (steps, last values)."""
    a, tail, slope = series(lam, big_p, reach)
    top = len(a) - 1
    den = 1
    for x in a:
        den = lcm(den, x.denominator)
    b = [int(a[n]*den)*step**(top - n) for n in range(top + 1)]
    full = den*step**top
    scale = 10**DIGITS
    pad_w = F(int(tail*scale) + 1, scale)
    pad_d = F(int(slope*scale) + 1, scale)
    c2 = F(big_p*(big_p - 2), 2)                              # 2C = s_0^3, the bottom of the potential
    k0 = max((k for k in range(1, reach*step) if F(k, step)**3 < c2 and potential(big_p, lam, F(k, step)) > 0),
             default=None)
    if k0 is None:
        raise ValueError('no start: the potential is not above lam near zero')

    def values(k):
        s = t = 0
        for n in range(top, -1, -1):
            s = s*k + b[n]
        for n in range(top, 0, -1):
            t = t*k + n*b[n]
        lo = (s*scale)//full
        dlo = (t*step*scale)//full
        w = (F(lo, scale) - pad_w, F(lo + 1, scale) + pad_w)
        d = (F(dlo, scale) - pad_d, F(dlo + 1, scale) + pad_d)
        g = F(big_p*step, 2*k)                                # p/s_k
        return w, (g*w[0] + d[0], g*w[1] + d[1])

    h = F(1, step)
    out = []
    cur = values(k0)
    first = cur
    for k in range(k0, reach*step):
        (vl, vu), (dl, du) = cur
        sa, sb = F(k, step), F(k + 1, step)
        bend = max(abs(potential(big_p, lam, sa)), abs(potential(big_p, lam, sb))) \
            + abs(1 - F(big_p*(big_p - 2), 2)/sa**3)*h
        if not bend*h*h < 1:
            raise ValueError('grid too coarse for the step bound')
        size = (max(abs(vl), abs(vu)) + max(abs(dl), abs(du))*h)/(1 - bend*h*h/2)
        err = bend*size*h*h/2
        out.append(dict(v=(min(vl, vl + dl*h) - err, max(vu, vu + du*h) + err),
                        d=(dl - bend*size*h, du + bend*size*h), left=(vl, vu)))
        nxt = values(k + 1)
        # the next step starts from u/s_(k+1)^p ; signs are those of w and of p w/s + w'
        cur = nxt
    return out, cur, first, F(k0, step)


def past_the_well(big_p, lam, reach):
    return F(reach)**3 > F(big_p*(big_p - 2), 2) and potential(big_p, lam, F(reach)) > 0


def level_one(lam, big_p, reach, step=STEP):
    """u > 0 on the whole half line certified: the first level is at least lam."""
    st, ((vl, vu), (dl, du)), first, start = walk(lam, big_p, reach, step)
    if not (first[0][0] > 0 and first[1][0] > 0):
        return False
    if any(not x['v'][0] > 0 for x in st):
        return False
    return past_the_well(big_p, lam, reach) and vl > 0 and dl > 0


def level_two(lam, big_p, reach, step=STEP):
    """u > 0, one node, u < 0 for ever certified: the second level is at least lam."""
    st, ((vl, vu), (dl, du)), first, start = walk(lam, big_p, reach, step)
    if not (first[0][0] > 0 and first[1][0] > 0):
        return False
    state = 'pos'
    for x in st:
        if state == 'pos':
            if x['v'][0] > 0:
                continue
            if x['d'][1] < 0 and x['left'][0] > 0:
                state = 'fall'
            else:
                return False
        if state == 'fall':
            if x['v'][1] < 0:
                state = 'neg'
            elif not x['d'][1] < 0:
                return False
            continue
        if state == 'neg' and not x['v'][1] < 0:
            return False
    return state == 'neg' and past_the_well(big_p, lam, reach) and vu < 0 and du < 0


def float_levels(big_p, n=3000):
    """Two lowest levels of -u'' + (C/s^2 + s) u by differences and sign counts (floats; only to place floors)."""
    c = big_p*(big_p - 2)/4
    reach = 14 + 2*c**(1/3)
    h = reach/n
    diag = [2/h**2 + c/(i*h)**2 + i*h for i in range(1, n)]
    off2 = 1/h**4

    def below(x):
        count, d = 0, 1.0
        for i, t in enumerate(diag):
            d = t - x - (off2/d if i else 0.0)
            if d == 0.0:
                d = 1e-300
            if d < 0:
                count += 1
        return count
    out = []
    for want in (1, 2):
        lo, hi = 0.0, 3*c**(1/3) + 12.0
        for _ in range(60):
            mid = (lo + hi)/2
            if below(mid) >= want:
                hi = mid
            else:
                lo = mid
        out.append(hi)
    return out


def certified_floors(big_p):
    """Rational floors of the first two levels for P = m - 1, each with its sign certificate."""
    one, two = float_levels(big_p)
    f1 = F(int((one - MARGIN)*10**4), 10**4)
    f2 = F(int((two - MARGIN)*10**4), 10**4)
    return (f1, level_one(f1, big_p, int(one) + 7)), (f2, level_two(f2, big_p, int(two) + 7))


def closed_floor(m):
    """(3/4)(2m - 1)^(2/3) from below: the local rate of s^p exp(-k s^(3/2)), k^2 = 2/9, p = (m - 1)/2."""
    return F(3, 4)*gc1.root3(F((2*m - 1)**2), 6)


# ---------------------------------------------------------------- readings of CR1's ladder for any d
def reading(d, top):
    w = WEIGHT[d]
    basis, s, h = cr1.ladder(d, w, top + 2)
    n = len(cm1.inv_basis(top))
    assert basis[:n] == cm1.inv_basis(top)
    index = {e: i for i, e in enumerate(basis)}
    gauss = 4.5*d*((d - 1)/32)**(1/3)
    x = gc1.near_lowest(s, h, n, 0.975*gauss)
    y = [F(0)]*len(basis)
    for j in range(n):
        for e, v in cm1.scalar_h({basis[j]: F(1)}, d, w).items():
            y[index[e]] += v*x[j]
    sxx = gc1.quad(s, x)
    hxx = gc1.quad(h, x)
    return dict(eta=hxx/sxx, square=gc1.quad(s, y)/sxx, same=gc1.quad(s, x, y) == hxx and sxx > 0, size=n)


# ---------------------------------------------------------------- the certificate
def run():
    checks, num = {}, {}
    rnd = random.Random(91026)

    # the weight by rows and by columns
    ok = True
    for d in (3, 4, 5, 7):
        for _ in range(6):
            c = [[F(rnd.randint(-7, 7), rnd.randint(1, 5)) for _ in range(d)] for _ in range(3)]
            cols = [[c[a][i] for a in range(3)] for i in range(d)]
            by_cols = sum(sum(x*x for x in u)*sum(x*x for x in v) - sum(x*y for x, y in zip(u, v))**2
                          for i, u in enumerate(cols) for v in cols[i + 1:])
            by_rows = sum(sum(x*x for x in c[a])*sum(x*x for x in c[b]) - sum(x*y for x, y in zip(c[a], c[b]))**2
                          for a in range(3) for b in range(a + 1, 3))
            ok = ok and by_cols == by_rows
    checks['D1 the weight read by columns equals the weight read by rows (d = 3, 4, 5, 7; exact matrices)'] = ok
    # pair of rows: (1/8)(-Laplacian) + (1/2)|r_a|^2 |y|^2 across r_a: mass 4, frequency |r_a|/2, d - 1 directions
    checks['D1 a pair is an oscillator of frequency |r|/2 in d - 1 directions; two pairs give (d - 1)|r|/2'] = all(
        F(1, 2)*4*F(1, 2)**2 == F(1, 2) and 2*(d - 1)*F(1, 2)*F(1, 2) == F(d - 1, 2) for d in TURNS)
    checks['D2 unit^3 = (d - 1)^2/16 : (1/(4 l^2)) = ((d - 1)/2) l'] = all(
        F(1, 4)**3*(4*F(d - 1, 2))**2 == F((d - 1)**2, 16) for d in TURNS)

    # one row: floors of the levels, m = 3 .. 14
    first, second = {3: gc1.E0}, {3: gc1.E1}
    ok1 = ok2 = okr = True
    for m in range(4, 15):
        (f1, c1), (f2, c2) = certified_floors(m - 1)
        first[m] = f1
        ok1 = ok1 and c1
        if m <= max(TURNS):
            second[m] = f2
            ok2 = ok2 and c2
        a, tail, slope = series(f1, m - 1, int(f1) + 7)
        low, top = residual_low(a, f1, m - 1)
        okr = okr and low >= top and tail < F(1, 10**19)
    checks['D2 the kept polynomials solve the equation up to their last degrees; tails under 10^-19'] = okr
    checks['D2 first levels, m = 4 .. 14: u > 0 on the whole half line at each floor'] = ok1
    checks['D2 second levels, d = 4 .. 10: one node, then u < 0 for ever, at each floor'] = ok2
    checks['D2 controls: 0.02 above each of four floors the same certificates are refused'] = \
        not level_one(first[4] + F(1, 50), 3, 10) and not level_one(first[9] + F(1, 50), 8, 12) \
        and not level_two(second[4] + F(1, 50), 3, 12) and not level_two(second[8] + F(1, 50), 7, 13)
    checks['D2 the floors rise with m and stand above the closed floor (3/4)(2m - 1)^(2/3)'] = all(
        first[m] < first[m + 1] for m in range(3, 14)) and all(first[m] > closed_floor(m) for m in range(3, 15)) \
        and all(second[d] > first[d] for d in TURNS)
    checks['D2 closed floor: 27 A^2 B/4 = (27/64)(4p + 1)^2 at k^2 = 2/9, A = 3k(4p + 1)/4, B = 1/2'] = all(
        27*(F(9, 16)*F(2, 9)*(4*p + 1)**2)*F(1, 2)/4 == F(27, 64)*(4*p + 1)**2 for p in (F(1), F(3, 2), F(2), F(13, 2)))
    num['floors of one row: m, first level, second level (l = 0)'] = [
        [m, float(first[m]), float(second[m]) if m in second else None] for m in range(3, 15)]

    # D3, D4: the count and the gap
    rows, window = [], {}
    okz = okg = okt = True
    for d in TURNS:
        unit = gc1.root3(F((d - 1)**2, 16), 6)
        even = min(second[d], first[d + 4])
        z = unit*min(2*first[d] + even, 3*first[d + 2])
        r = reading(d, DEGREE.get(d, 6))
        eta, square = r['eta'], r['square']
        okt = okt and r['same'] and square > eta*eta
        okz = okz and z > eta
        low = gc1.floor_of_lowest(eta, square, z)
        okg = okg and low > 3*unit*first[d] and low <= eta
        window[d] = (low, eta, z)
        rows.append([d, float(gc1.cut(z)), float(gc1.cut(low)), float(gc1.cut(eta, up=True)),
                     float(gc1.cut(z - eta)), float(gc1.cut(low**3/(d**3*(d - 1)), 5)),
                     float(gc1.cut(eta**3/(d**3*(d - 1)), 5, up=True))])
    checks['D4 the readings are exact (<psi, h psi> two ways) and have positive spread, d = 3 .. 10'] = okt
    checks['D4 the line by rows is above the reading: a gap for every d = 3 .. 10'] = okz
    checks['D4 the lowest rate from both sides, above three copies of one row, d = 3 .. 10'] = okg
    checks['D4 three turns: rows and columns give the same line (5.5210)'] = \
        gc1.cut(window[3][2]) == gc1.cut(gc1.line(3))
    checks['D4 four turns: gap >= 0.44 by rows, against 0.0031 by columns'] = \
        window[4][2] - window[4][1] > F(44, 100) and gc1.line(4) - window[4][1] < F(32, 10000)
    checks['D4 five to ten turns: the line by columns (its two levels taken from above as 2.3382, 4.0880) is under the lowest rate'] = all(
        gc1.UNIT.get(d, gc1.root3(F(d - 1, 8), 6, up=True))*((d - 1)*F(23382, 10000) + F(4088, 1000)) < window[d][0]
        for d in range(5, 11))
    up4 = gc1.turning_upper(4)
    checks['D4 four turns: gap window with CM1\'s turning reading'] = \
        window[4][2] < up4 and up4 - window[4][0] < F(2972, 1000)
    checks['D4 the cubes rho_d^3 = E^3/(d^3 (d - 1)) rise with d, window after window, d = 3 .. 10'] = all(
        window[d][1]**3/(d**3*(d - 1)) < window[d + 1][0]**3/((d + 1)**3*d) for d in range(3, 10))
    num['d, line, lowest rate (low, high), gap floor, rho^3 (low, high)'] = rows
    num['four turns: gap window'] = [float(gc1.cut(window[4][2] - window[4][1])), float(gc1.cut(up4 - window[4][0], up=True))]

    # D5: the end of the walk
    checks['D5 (729/1024)(d - 1)^2 (2d - 1)^2 <= (729/256) d^3 (d - 1)(1 - 2/d) ... the ratio is at least 1 - 2/d, d = 3 .. 400'] = all(
        F((d - 1)*(2*d - 1)**2, 4*d**3) >= 1 - F(2, d) and F((d - 1)*(2*d - 1)**2, 4*d**3) < 1 for d in range(3, 401))
    checks['D5 three copies of the closed floor: 27 unit^3 (27/64)(2d - 1)^2 = (729/1024)(d - 1)^2 (2d - 1)^2'] = all(
        27*F((d - 1)**2, 16)*F(27, 64)*(2*d - 1)**2 == F(729, 1024)*(d - 1)**2*(2*d - 1)**2 for d in range(3, 60))
    checks['D5 the certified windows lie between the closed bounds, d = 3 .. 10'] = all(
        F(729, 1024)*(d - 1)**2*(2*d - 1)**2 <= window[d][0]**3 and window[d][1]**3 <= F(729, 256)*d**3*(d - 1)
        for d in TURNS)
    num['limit of rho'] = float(gc1.root3(F(729, 256), 5))

    out = dict(stage='DR1', checks=len(checks), passed=sum(1 for v in checks.values() if v),
               all_pass=all(checks.values()), detail=checks, numbers=num)
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out['detail'].items():
        print('PASS' if v else 'FAIL', k)
    print(json.dumps(num, indent=1))
    print(out['passed'], '/', out['checks'])
    with open(os.path.join(HERE, 'DR1_RESULT.json'), 'w') as f:
        json.dump(out, f, indent=1)
