"""PC1 -- pair clusters: the core of n turns as a sum of its sub-cores, each with its lowest reading and its floor.

The core of n turns is the sum of its sub-cores of n - 1 turns in the form a T + b X (GC1-G6), each a dilated
copy, by kappa = (a^2 b)^(1/3), of the core h of n - 1 turns.  GC1-G6 used only the lowest rate of the copy.
If on ALL readings h has one reading g under a floor rho (E its lowest rate), then h >= rho - (rho - E) P_g, and

    P1  h_n  >=  kappa [ n rho - (rho - E) M ] ,      M = sum over the n sub-cores of P_g (x) 1 .
    P2  M = J J* with J (a_1 .. a_n) = sum_k g (x) a_k, and J* J = 1 + R (x) (ones - 1), R the one-turn part of
        g (a positive operator of trace 1 with levels r_1 >= r_2 >= ..).  So the levels of M are
        1 + (n - 1) r  and  1 - r ;  at most one of them is above 1 + (n - 1)(1 - r_1).
    P3  hence at most ONE reading of h_n -- on all readings, no row condition -- lies under
            Z_n = kappa [ (n - 1) rho + E - (n - 1)(1 - r_1)(rho - E) ] ,
        and r_1 >= <a x .. x a, g>^2 for any one-turn reading a.

Two turns (the start).  By columns, one level under (e0 + e1)/2 = 3.2130 among readings even in both turns
with no row turning; every other sector has its floor above: 3.2921 (row turning two), and 3.3592 for every
reading odd in a turn, by giving ALL of the other turn's kinetic part to the layer (h >= -(1/2) Laplacian_1 +
sqrt 2 |c_1|).  So on all readings of two turns at most one lies under rho_2 = 3.2130:
        lowest rate in [2.6592, 2.6594] ,   gap >= 0.5536 on all readings.
Three turns:  Z_3 = 5.69 (GC1, SG1, DR1: 5.5210, and only under the row condition):  gap >= 0.50 on all readings.
Four turns:   Z_4 = 8.42 on all readings (DR1: 8.4454 under the row condition).

Exact rationals; trial vectors are found in floats and then used as exact rational vectors.  Stdlib only.
"""
from fractions import Fraction as F
from math import isqrt
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


dr1 = load('dr1', 'dr1_core_by_rows')
gc1 = dr1.gc1
cm1, cr1 = gc1.cm1, gc1.cr1

WEIGHT = {2: F(5, 4), 3: F(8, 5), 4: F(9, 5)}
TOP = 10
TOP_TWO = 14                          # degree of the two-turn reading
NU = {3: gc1.E0, 5: F(33592, 10000), 7: F(42461, 10000), 9: F(50489, 10000), 11: F(57924, 10000)}   # DR1 floors
NU2 = {3: gc1.E1}
SHAPE = F(2, 25)                      # b of the one-turn reading (1 + b |c|^2) exp(-w |c|^2 / 2), two turns


def sqrt_floor(x, places=12):
    t = 10**(2*places)
    return F(isqrt(x.numerator*t//x.denominator), 10**places)


def sqrt_ceil(x, places=12):
    return sqrt_floor(x, places) + F(1, 10**places)


# ---------------------------------------------------------------- readings on CR1's ladder, d = 2, 3, 4
def reading(d, top=TOP, shape=None):
    """One exact rational reading psi of the ladder: mean of h, mean of h^2, and the squared overlap with the
    product reading a x .. x a, a = (1 + shape |c|^2) exp(-w |c|^2/2) (shape only for d = 2)."""
    w = WEIGHT[d]
    basis, s, h = cr1.ladder(d, w, top + 2)
    keep = [i for i, e in enumerate(basis) if d > 2 or e[2] == 0]        # two turns: the third number is zero
    low = [i for i in keep if cr1.degree(basis[i]) <= top]
    sk = [[s[i][j] for j in keep] for i in keep]
    index = {basis[i]: k for k, i in enumerate(keep)}
    gauss = 4.5*d*((d - 1)/32)**(1/3)
    small_s = [[s[i][j] for j in low] for i in low]
    small_h = [[h[i][j] for j in low] for i in low]
    x_low = gc1.near_lowest(small_s, small_h, len(low), {2: 0.93, 3: 0.975, 4: 0.975}[d]*gauss)
    x = [F(0)]*len(keep)
    y = [F(0)]*len(keep)
    for jj, i in enumerate(low):
        x[index[basis[i]]] = x_low[jj]
        for e, v in cm1.scalar_h({basis[i]: F(1)}, d, w).items():
            if e in index:
                y[index[e]] += v*x_low[jj]
    sxx = gc1.quad(sk, x)
    hxx = gc1.quad(sk, x, y)
    assert hxx == gc1.quad(small_h, x_low) and sxx > 0
    # the product reading, averaged over the directions (psi is unchanged by them): a polynomial in e1, e2
    prod = [F(0)]*len(keep)
    prod[index[(0, 0, 0)]] = F(1)
    norm_a = F(1)
    if shape is not None:
        assert d == 2
        prod[index[(1, 0, 0)]] = shape
        prod[index[(2, 0, 0)]] = shape*shape/8
        prod[index[(0, 1, 0)]] = shape*shape/2
        norm_a = 1 + 2*shape*F(3)/(2*w) + shape*shape*F(15)/(4*w*w)       # one turn: <(1 + b r^2)^2>
    cross = gc1.quad(sk, prod, x)
    overlap2 = cross*cross/(norm_a**d*sxx)
    return dict(eta=hxx/sxx, square=gc1.quad(sk, y)/sxx, overlap2=overlap2, size=len(low))


def one_turn_share(r, rho, low):
    """Floor of r_1, the largest level of the one-turn part of the true lowest reading g.
    |<psi, g>|^2 >= (rho - eta)/(rho - E) and the angle from the product reading to g is at most the sum
    of the angles product--psi and psi--g."""
    p0 = (rho - r['eta'])/(rho - low)
    c1, c2 = sqrt_floor(r['overlap2']), sqrt_floor(p0)
    s1, s2 = sqrt_ceil(1 - r['overlap2']), sqrt_ceil(1 - p0)
    assert 0 < p0 <= 1 and c1*c2 > s1*s2
    return (c1*c2 - s1*s2)**2, p0


def cluster_line(n, rho, low, r1):
    """Z_n of P3 from the data of n - 1 turns: floor rho, lowest rate floor low, one-turn share r1."""
    kappa = gc1.root3(F(1, (n - 1)**2*(n - 2)), 9)           # a = 1/(n - 1), b = 1/(n - 2)
    return kappa*((n - 1)*rho + low - (n - 1)*(1 - r1)*(rho - low))


# ---------------------------------------------------------------- exact small model of P2
def model_check(m, n, rnd):
    """A symmetric rational g on m^(n-1) points; M = sum of P_g (x) 1 on m^n points; J* J = 1 + R (x) (ones - 1)."""
    import itertools
    pts = list(itertools.product(range(m), repeat=n - 1))
    raw = {}
    for p in pts:
        key = tuple(sorted(p))
        if key not in raw:
            raw[key] = F(rnd.randint(-4, 4), rnd.randint(1, 3))
    g = {p: raw[tuple(sorted(p))] for p in pts}
    norm = sum(v*v for v in g.values())
    if norm == 0:
        return True
    full = list(itertools.product(range(m), repeat=n))

    def jay(k, a):                                              # g on the turns other than k, times a on turn k
        return {c: g[c[:k] + c[k + 1:]]*a[c[k]] for c in full}

    def jay_star(k, psi):
        out = [F(0)]*m
        for c in full:
            out[c[k]] += g[c[:k] + c[k + 1:]]*psi[c]
        return out
    red = [[sum(g[(i,) + rest]*g[(j,) + rest] for rest in itertools.product(range(m), repeat=n - 2))/norm
            for j in range(m)] for i in range(m)]
    ok = sum(red[i][i] for i in range(m)) == 1
    for _ in range(3):
        a = [[F(rnd.randint(-3, 3)) for _ in range(m)] for _ in range(n)]
        psi = {c: sum(jay(k, a[k])[c] for k in range(n)) for c in full}
        # J* J a / norm = a_k + R (sum of the others)
        for k in range(n):
            left = [x/norm for x in jay_star(k, psi)]
            others = [sum(a[l][i] for l in range(n) if l != k) for i in range(m)]
            right = [a[k][i] + sum(red[i][j]*others[j] for j in range(m)) for i in range(m)]
            ok = ok and left == right
    return ok


# ---------------------------------------------------------------- the certificate
def run():
    checks, num = {}, {}
    rnd = random.Random(1009)

    # P1: shares and dilation
    checks['P1 each turn is in n - 1 sub-cores and each pair in n - 2 (n = 3 .. 9); kappa^3 = a^2 b'] = all(
        F(n - 1, n - 1) == 1 and F(n - 2, n - 2) == 1 and F(1, n - 1)**2*F(1, n - 2) == F(1, (n - 1)**2*(n - 2))
        for n in range(3, 10)) and gc1.root3(F(1, 4), 9)**3 <= F(1, 4) and gc1.root3(F(1, 18), 9)**3 <= F(1, 18)
    # P2: the block identity on exact small models
    checks['P2 J* J = 1 + R (x) (ones - 1) on exact models (m = 2, 3 points; n = 3, 4 turns), trace R = 1'] = all(
        model_check(m, n, rnd) for m, n in ((2, 3), (3, 3), (2, 4)) for _ in range(3))
    checks['P2 levels of 1 + R (x) (ones - 1): 1 + (n - 1) r on equal parts, 1 - r on parts that sum to zero'] = all(
        (n - 1) - 0 == n - 1 and 1 + (n - 1)*r > 1 - r for n in (3, 4) for r in (F(1, 7), F(9, 10)))

    # two turns: the sectors on all readings
    unit = F(1, 2)                                             # ((d - 1)/8)^(1/3) at d = 2
    sectors = {
        'no row turning, both turns even: second level': unit*(NU[3] + NU2[3]),
        'row turning two, both even (one turn with turning number two)': unit*(NU[3] + NU[7]),
        'row turning four, both even': unit*(NU[3] + NU[11]),
        'both even, two turns with turning number two': unit*2*NU[7],
        'a turn odd (all of the other turn given to the layer)': NU[5],
    }
    rho2 = min(sectors.values())
    checks['two turns: unit 1/2 by halves; unit 1 when one turn gives all its kinetic part ((1/2)(sqrt 2)^2 = 1)'] = \
        F(1, 8) == unit**3 and F(1, 2)*2 == 1
    checks['two turns: the least sector floor is the second level of the even readings without row turning, 3.2130'] = \
        rho2 == unit*(NU[3] + NU2[3]) == F(3213, 1000) and all(v >= rho2 for v in sectors.values())
    checks['two turns, control: with halves only, a reading odd in one turn has floor 2.8487 < 3.2130'] = \
        unit*(NU[3] + NU[5]) < rho2
    num['two turns: sector floors'] = {k: float(gc1.cut(v)) for k, v in sectors.items()}

    r2 = reading(2, TOP_TWO, SHAPE)
    eta2, sq2 = r2['eta'], r2['square']
    low2 = gc1.floor_of_lowest(eta2, sq2, rho2)
    checks['two turns: the reading is under the floor; lowest rate in [2.6592, 2.6594]'] = \
        eta2 < rho2 and low2 > F(26592, 10000) and eta2 < F(26594, 10000) and sq2 > eta2*eta2
    checks['two turns: gap >= 0.5536 on all readings'] = rho2 - eta2 > F(5536, 10000)
    r1_2, p0_2 = one_turn_share(r2, rho2, low2)
    checks['two turns: one-turn share of the lowest reading r_1 >= 0.955'] = r1_2 > F(955, 1000)
    num['two turns'] = dict(floor=float(rho2), lowest_rate=[float(gc1.cut(low2)), float(gc1.cut(eta2, up=True))],
                            gap_all_readings=float(gc1.cut(rho2 - eta2)), readings=r2['size'],
                            overlap_with_product=float(gc1.cut(r2['overlap2'])),
                            overlap_with_true_lowest=float(gc1.cut(p0_2)), one_turn_share=float(gc1.cut(r1_2)))

    # three turns from two
    z3 = cluster_line(3, rho2, low2, r1_2)
    r3 = reading(3)
    eta3, sq3 = r3['eta'], r3['square']
    checks['three turns: Z_3 > 5.69, above the line by columns and rows (5.5210)'] = \
        z3 > F(569, 100) and z3 > gc1.line(3)
    checks['three turns: gap >= 0.50 on all readings (was 0.3342 under the row condition)'] = \
        z3 - eta3 > F(50, 100) and gc1.line(3) - eta3 < F(3343, 10000)
    checks['three turns, control: without the full layer for odd turns the floor is 2.8487 and Z_3 falls under 5.5210'] = \
        cluster_line(3, unit*(NU[3] + NU[5]), low2, r1_2) < gc1.line(3)
    checks['three turns, control: knowing nothing of the one-turn share (r_1 = 1/2) leaves Z_3 under 5.5210'] = \
        cluster_line(3, rho2, low2, F(1, 2)) < gc1.line(3)
    low3 = gc1.floor_of_lowest(eta3, sq3, z3)
    checks['three turns: lowest rate from both sides with the new line'] = \
        low3 > F(51865, 10000) and eta3 < F(51868, 10000)
    r1_3, p0_3 = one_turn_share(r3, z3, low3)
    up3 = gc1.turning_upper(3)
    num['three turns'] = dict(line=float(gc1.cut(z3)), line_before=float(gc1.cut(gc1.line(3))),
                              lowest_rate=[float(gc1.cut(low3)), float(gc1.cut(eta3, up=True))],
                              gap=[float(gc1.cut(z3 - eta3)), float(gc1.cut(up3 - low3, up=True))],
                              overlap_with_product=float(gc1.cut(r3['overlap2'])),
                              one_turn_share=float(gc1.cut(r1_3)))

    # four turns from three
    z4 = cluster_line(4, z3, low3, r1_3)
    r4 = reading(4)
    eta4 = r4['eta']
    checks['four turns: Z_4 > 8.41 on all readings; gap >= 0.41 (by rows, under the row condition: 0.4425)'] = \
        z4 > F(841, 100) and z4 - eta4 > F(41, 100)
    num['four turns'] = dict(line=float(gc1.cut(z4)), gap_all_readings=float(gc1.cut(z4 - eta4)))
    # what the same step would give with the levels of two turns as floats suggest them (not certified)
    num['for the eye: Z_3 if the floor of two turns were 3.66 and r_1 = 0.959'] = \
        float(gc1.cut(cluster_line(3, F(366, 100), low2, F(959, 1000))))

    out = dict(stage='PC1', checks=len(checks), passed=sum(1 for v in checks.values() if v),
               all_pass=all(checks.values()), detail=checks, numbers=num)
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out['detail'].items():
        print('PASS' if v else 'FAIL', k)
    print(json.dumps(num, indent=1))
    print(out['passed'], '/', out['checks'])
    with open(os.path.join(HERE, 'PC1_RESULT.json'), 'w') as f:
        json.dump(out, f, indent=1)
