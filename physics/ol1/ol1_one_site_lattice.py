"""OL1 -- the one-site lattice: the count of GC1 on compact links, at strong and at arbitrarily weak coupling.

Three links U_1, U_2, U_3 in SU(2) on the periodic lattice of one site (three plaquettes).  A link is a unit
block (cos psi, sin psi m), m a unit direction, 0 <= psi <= pi; its three cut components are u = sin psi m.
    O1  1 - W_ij = 2 |u_i cross u_j|^2 ,  W_ij = (1/2) trace of U_i U_j U_i^-1 U_j^-1          (HL1-H7)
    so  H = - sum of the sphere Laplacians + 2 theta sum over pairs |u_i cross u_j|^2   on three unit spheres.
        The Wilson operator of the YM line is A = H/4 - 3 theta_YM with theta = 4 theta_YM.
    O2  the layer on the sphere:  - Laplacian + mu (x_1^2 + x_2^2)  >=  2 (sqrt(4 + mu) - 2)
        (local rate of exp(-kappa t), t = x_1^2 + x_2^2, with mu = 8 kappa + 4 kappa^2).
    O3  the comparison:  H >= k_1 + k_2 + k_3 ,   k = -(1/2) Laplacian + sqrt(4 + 4 theta sin^2 psi) - 2 .
    O4  the count: among readings unchanged by the row turning (one conjugation of all links) and by the
        three centre flips U_i -> -U_i, at most one lies under  z = min(2 e0 + e1, e0 + 2 p0),
        e0, e1 the first two levels of k on class readings even about psi = pi/2, p0 its floor on readings
        with zero average over conjugation.
    O5  strong side, 0 <= theta <= 15/4 :  e0 from the local rate of exp(-kappa sin^2 psi), e1 >= 4,
        p0 from the local rate of sin psi exp(-kappa sin^2 psi), the lowest rate from the reading
        1 + b cos 2 psi on each link:  gap > 0 on the whole window (steps of 1/8, monotone in theta).
    O6  weak side, every theta >= 10^7 :  sin psi >= sigma min(psi, a), sigma = sin a / a, turns each link into
        the half line of GC1-G1 cut at s_a >= 9; the same two sign certificates and the same zero-average
        reading (continued by a hyperbolic cosine) give
            z >= (2 e0 + e1)(2 theta sigma^2)^(1/3) - 15/2 ,   lowest rate <= (27/2)(theta/2)^(1/3) + 3 ,
        so   gap >= 0.1312 theta^(1/3) - 21/2 > 0 ,  and  gap/theta^(1/3) is at least 0.327 in the limit.

Exact rationals; parameters of the readings are found in floats and then used as exact rationals.  Stdlib only.
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


gc1 = load('gc1', 'gc1_core_gap_count')
hl1 = load('hl1', 'hl1_helical_walk')

E0, E1, P0 = gc1.E0, gc1.E1, gc1.P0          # floors of one turn on the half line (GC1-G1)
STRONG_END = F(15, 4)
STRONG_STEP = F(1, 8)
CELLS = 1000
WEAK_START = 10**7
CUT = F(2, 5)                                 # the angle a of the weak side
REACH = 9                                     # s_a must pass the reach of both sign certificates of GC1
JOIN, RATE = 5, F(6, 5)                       # the zero-average reading is continued at s = 5 by cosh(RATE (b - s))


# ---------------------------------------------------------------- small exact tools
def sqrt_floor(x, places=12):
    """A rational not above the square root of the rational x >= 0."""
    t = 10**(2*places)
    return F(isqrt(x.numerator*t//x.denominator), 10**places)


def dfact(n):
    """n!! with (-1)!! = 0!! = 1"""
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


def link_mean(a, b):
    """(1/pi) x integral over [0, pi] of cos^a psi sin^b psi, a and b even."""
    return F(dfact(a - 1)*dfact(b - 1), dfact(a + b))


# ---------------------------------------------------------------- O1, O2: the weight and the layer on the sphere
def plaquette(u, v):
    """Scalar part of U V U^-1 V^-1 for unit blocks."""
    return hl1.qmul(hl1.qmul(u, v), hl1.qmul(hl1.qconj(u), hl1.qconj(v)))[0]


def sphere_laplacian_of_power(m):
    """Laplacian on the unit 3-sphere of t^m, t = x_1^2 + x_2^2, as a polynomial in t (harmonic route)."""
    out = {m: F(-4*m*(m + 1))}
    if m:
        out[m - 1] = F(4*m*m)
    return {e: c for e, c in out.items() if c}


def layer_formula_of_power(m):
    """4 [ t (1 - t) f'' + (1 - 2t) f' ] on f = t^m."""
    out = {}
    for e, c in ((m - 1, 4*m*(m - 1) + 4*m), (m, -4*m*(m - 1) - 8*m)):
        if c and e >= 0:
            out[e] = out.get(e, 0) + F(c)
    return out


def layer_floor(mu):
    """2 (sqrt(4 + mu) - 2), from below."""
    return 2*(sqrt_floor(4 + mu) - 2)


# ---------------------------------------------------------------- O5: local rates and readings on one link
def local_floor(theta, kappa, lead, slope, cells=CELLS):
    """Least value over y = sin^2 psi in [0, 1], from below, of
           lead - slope y - 2 kappa^2 y (1 - y) + 2 sqrt(1 + theta y) - 2 ."""
    low = None
    for i in range(cells):
        ya, yb = F(i, cells), F(i + 1, cells)
        if yb <= F(1, 2):
            bump = yb*(1 - yb)
        elif ya >= F(1, 2):
            bump = ya*(1 - ya)
        else:
            bump = F(1, 4)
        value = lead - slope*yb - 2*kappa*kappa*bump + 2*sqrt_floor(1 + theta*ya) - 2
        low = value if low is None else min(low, value)
    return low


def class_floor(theta, kappa):
    """Floor of the first level of k on class readings: local rate of exp(-kappa sin^2 psi)."""
    return local_floor(theta, kappa, 3*kappa, 4*kappa)


def rest_floor(theta, kappa):
    """Floor of k on readings with zero conjugation average: local rate of sin psi exp(-kappa sin^2 psi)."""
    return local_floor(theta, kappa, F(3, 2) + 5*kappa, 6*kappa)


def best_kappa(theta, lead, slope, top=3.0, n=300, grid=200):
    """A float search for the kappa of a local rate; only the exact floor at the result is used."""
    th = float(theta)
    best, arg = None, 0.0
    for i in range(n + 1):
        k = top*i/n
        low = min(lead(k) - slope(k)*y - 2*k*k*y*(1 - y) + 2*(1 + th*y)**0.5 - 2
                  for y in (j/grid for j in range(grid + 1)))
        if best is None or low > best:
            best, arg = low, k
    return F(round(arg*1000), 1000)


def two_term(theta, b):
    """Mean of H in the reading 1 + b cos 2 psi on each of the three links (exact)."""
    norm = 2 - 2*b + b*b
    kinetic = 4*b*b/norm
    spread = (12 - 16*b + 7*b*b)/(8*norm)
    return 3*kinetic + 4*theta*spread*spread


def two_term_by_integrals(b):
    """Kinetic mean and mean of sin^2 psi for 1 + b cos 2 psi = (1 - b) + 2b cos^2 psi, from link_mean."""
    c0, c1 = 1 - b, 2*b
    norm = c0*c0*link_mean(0, 2) + 2*c0*c1*link_mean(2, 2) + c1*c1*link_mean(4, 2)
    spread = c0*c0*link_mean(0, 4) + 2*c0*c1*link_mean(2, 4) + c1*c1*link_mean(4, 4)
    kinetic = 16*b*b*link_mean(2, 4)                      # derivative -4b cos psi sin psi, weight sin^2 psi
    return kinetic/norm, spread/norm


def best_b(theta, n=400):
    best, arg = None, F(0)
    for i in range(n):
        b = F(i, n)
        e = two_term(theta, b)
        if best is None or e < best:
            best, arg = e, b
    return arg


def power_reading(theta, n):
    """Mean of H in the reading cos^(2n) psi on each link: 36 n^2/(4n - 1) + (9 theta/4)/(n + 1)^2."""
    return F(36*n*n, 4*n - 1) + F(9, 4)*theta/(n + 1)**2


def power_reading_by_integrals(n):
    kinetic = 4*n*n*link_mean(4*n - 2, 4)/link_mean(4*n, 2)
    spread = link_mean(4*n, 4)/link_mean(4*n, 2)
    return kinetic, spread


# ---------------------------------------------------------------- the certificate
def run():
    checks, num = {}, {}
    rnd = random.Random(20261009)

    # O1
    ok = True
    for _ in range(40):
        u = hl1.unit_block([F(rnd.randint(-9, 9), rnd.randint(1, 9)) for _ in range(3)])
        v = hl1.unit_block([F(rnd.randint(-9, 9), rnd.randint(1, 9)) for _ in range(3)])
        ok = ok and plaquette(u, v) == 1 - 2*hl1.cross_sq(u, v) \
            and sum(x*x for x in u) == 1 and sum(x*x for x in u[1:]) == 1 - u[0]*u[0]
    checks['O1 1 - W = 2 |u cross v|^2 on 40 rational pairs of links; |u|^2 = sin^2 psi'] = ok

    # O2
    checks['O2 Laplacian of the 3-sphere on functions of t is 4[t(1-t) f\'\' + (1-2t) f\'] (t^m, m = 0 .. 8)'] = all(
        sphere_laplacian_of_power(m) == layer_formula_of_power(m) for m in range(9))
    ok = True
    for kappa in (F(1, 7), F(1, 2), F(3, 2), F(5), F(40)):
        mu = 8*kappa + 4*kappa*kappa
        # local rate of exp(-kappa t): 4 kappa + (mu - 8 kappa - 4 kappa^2) t + 4 kappa^2 t^2
        for t in (F(0), F(1, 3), F(1, 2), F(9, 10), F(1)):
            rate = -(4*kappa*kappa*t*(1 - t) - 4*kappa*(1 - 2*t)) + mu*t
            ok = ok and rate == 4*kappa + 4*kappa*kappa*t*t and rate >= 4*kappa
        ok = ok and (2*kappa + 2)**2 == 4 + mu and 0 <= 4*kappa - layer_floor(mu) < F(1, 10**9)
    checks['O2 local rate of exp(-kappa t) is 4 kappa + 4 kappa^2 t^2 at mu = 8 kappa + 4 kappa^2; 4 kappa = 2(sqrt(4 + mu) - 2)'] = ok

    # O3: three links, each pair gets 1/4 of a Laplacian: (d - 1) a W(theta s^2/a) with a = 1/4 is sqrt(4 + 4 theta s^2) - 2
    ok = True
    for kappa in (F(1, 5), F(2), F(11)):
        s2 = F(3, 7)
        theta = (8*kappa + 4*kappa*kappa)*F(1, 4)/s2          # so that theta s^2/a = 8 kappa + 4 kappa^2
        ok = ok and 2*F(1, 4)*4*kappa == 2*kappa and (2*kappa + 2)**2 == 4 + 4*theta*s2
    checks['O3 the comparison potential of one link is sqrt(4 + 4 theta sin^2 psi) - 2'] = ok

    # readings: exact integrals
    checks['O5 the reading 1 + b cos 2 psi: kinetic 4b^2/(2 - 2b + b^2), spread (12 - 16b + 7b^2)/(8(2 - 2b + b^2))'] = all(
        two_term_by_integrals(b) == (4*b*b/(2 - 2*b + b*b), (12 - 16*b + 7*b*b)/(8*(2 - 2*b + b*b)))
        for b in (F(0), F(1, 9), F(1, 3), F(1, 2), F(7, 8), F(2)))
    checks['O5 at b = 2 the reading is the class function of weight two: kinetic 8'] = \
        two_term_by_integrals(F(2))[0] == 8
    checks['O6 the reading cos^(2n) psi: kinetic 12 n^2/(4n - 1), spread 3/(4n + 4), n = 1 .. 40'] = all(
        power_reading_by_integrals(n) == (F(12*n*n, 4*n - 1), F(3, 4*n + 4)) for n in range(1, 41))

    # O5: the strong window
    steps = int(STRONG_END/STRONG_STEP)
    margins, table = [], []
    for k in range(steps):
        theta, nxt = k*STRONG_STEP, (k + 1)*STRONG_STEP
        e0 = max(F(0), class_floor(theta, best_kappa(theta, lambda x: 3*x, lambda x: 4*x)))
        p0 = max(F(3, 2), rest_floor(theta, best_kappa(theta, lambda x: 1.5 + 5*x, lambda x: 6*x)))
        z = min(2*e0 + 4, e0 + 2*p0)
        up = two_term(nxt, best_b(nxt))
        margins.append(z - up)
        if k % 6 == 0 or k == steps - 1:
            table.append([float(theta), float(gc1.cut(e0)), float(gc1.cut(p0)), float(gc1.cut(z)),
                          float(gc1.cut(up, up=True)), float(gc1.cut(z - up))])
    checks['O5 strong side: the line at theta_k is above the reading at theta_k + 1/8, every step up to 15/4'] = \
        all(m > 0 for m in margins)
    checks['O5 at theta = 0 the count gives 3 (the line) against 0 (the lowest rate)'] = \
        min(2*class_floor(F(0), F(0)) + 4, class_floor(F(0), F(0)) + 2*rest_floor(F(0), F(0))) == 3
    num['strong side'] = dict(window=[0.0, float(STRONG_END)], least_margin=float(gc1.cut(min(margins))),
                              rows_theta_e0_p0_line_reading_margin=table)

    # O6: the weak side
    a = CUT
    sin_low = a - a**3/6 + a**5/120 - a**7/5040
    sigma = sin_low/a
    checks['O6 sin a >= a - a^3/6 + a^5/120 - a^7/5040 (alternating, falling terms) ; sigma = sin a/a > 0.9735'] = \
        a**9/362880 < a**7/5040 < a**5/120 < a**3/6 < a and F(9735, 10000) < sigma < 1
    reach3 = 4*sqrt_floor(F(WEAK_START))*a*a*sin_low
    checks['O6 the cut of the half line is past the reach of both certificates: s_a^3 = 4 sqrt(theta) a^2 sin a >= 9^3'] = \
        reach3 >= REACH**3 and REACH >= max(gc1.REACH.values()) and a < F(3, 2)
    checks['O6 the two sign certificates of GC1-G1 (v > 0 for ever; one node, v < 0 for ever)'] = \
        gc1.level_one(E0) and gc1.level_two(E1)[0]
    # zero average: s^2 exp(-k s^(3/2)) up to s = 5, then cosh(RATE (b - s)), b >= s_a >= 9
    exp_low = sum(F(48, 5)**j/gc1.factorial(j) for j in range(40))        # exp(9.6) from below
    tanh_low = 1 - F(2)/exp_low                                            # tanh(4.8) from below
    slope_need = F(31624, 20000) - F(2, 5)                                 # (3k/2) sqrt 5 - 2/5 = sqrt 10/2 - 2/5, from above
    checks['O6 zero average: past s = 5 the rate 5 - (6/5)^2 is above p0; the join bends downward'] = \
        JOIN - RATE*RATE >= P0 and F(31624, 10000)**2 >= 10 and RATE*tanh_low >= slope_need \
        and RATE*(REACH - JOIN) >= F(24, 5) and P0**3 <= F(2187, 64) and REACH >= JOIN
    two_low = gc1.root3(F(2), 6)
    half_up = gc1.root3(F(1, 2), 6, up=True)
    sig_low = gc1.root3(sigma*sigma, 6)
    count = min(2*E0 + E1, E0 + 2*P0)
    c_line = count*two_low*sig_low
    c_read = F(27, 2)*half_up
    checks['O6 line: z >= c theta^(1/3) - 15/2 with c > 10.846 ; reading: lowest rate <= c\' theta^(1/3) + 3 with c\' < 10.7151'] = \
        c_line > F(10846, 1000) and c_read < F(107151, 10000) and count == 2*E0 + E1
    checks['O6 36 n^2 = (9n + 9/4)(4n - 1) + 9/4 and (9/4)/(4n - 1) <= 3/4 , n = 1 .. 60'] = all(
        36*n*n == (9*n + F(9, 4))*(4*n - 1) + F(9, 4) and F(9, 4)/(4*n - 1) <= F(3, 4) for n in range(1, 61))
    slope = c_line - c_read
    third_low = gc1.root3(F(WEAK_START), 4)
    checks['O6 gap >= 0.1312 theta^(1/3) - 21/2 , positive from theta = 10^7 on, and >= 0.08 theta^(1/3) there'] = \
        slope > F(1312, 10000) and slope*third_low - F(21, 2) > 0 and (slope - F(8, 100))*third_low >= F(21, 2)
    spots = []
    ok = True
    for theta in (10**7, 10**9, 10**12, 10**18):
        x_low = gc1.root3(F(theta, 2), 6)
        n = int(x_low)
        reading = power_reading(F(theta), n)
        t_low, t_up = gc1.root3(F(theta), 6), gc1.root3(F(theta), 6, up=True)
        line = c_line*t_low - F(15, 2)
        ok = ok and n >= 1 and reading <= c_read*t_up + 3 and line > reading
        spots.append([theta, n, float(gc1.cut(line)), float(gc1.cut(reading, up=True)),
                      float(gc1.cut((line - reading)/t_up))])
    checks['O6 spot values: the line above the exact power reading at theta = 10^7, 10^9, 10^12, 10^18'] = ok
    limit = count*two_low - c_read
    checks['O6 in the limit gap/theta^(1/3) >= (2 e0 + e1) 2^(1/3) - (27/2) 2^(-1/3) > 0.327'] = limit > F(327, 1000)
    num['weak side'] = dict(start=WEAK_START, a=float(a), sigma=float(gc1.cut(sigma)),
                            line_coefficient=float(gc1.cut(c_line)), reading_coefficient=float(gc1.cut(c_read, up=True)),
                            gap_coefficient=float(gc1.cut(slope)), constant=10.5,
                            limit_coefficient=float(gc1.cut(limit)),
                            spots_theta_n_line_reading_gap_over_cube_root=spots)
    # the dictionary to the Wilson operator of the YM line: A = H/4 - 3 theta_YM, theta = 4 theta_YM
    num['in the units of the YM line'] = dict(
        strong_window_theta=[0.0, float(STRONG_END/4)], weak_start_theta=float(F(WEAK_START, 4)),
        weak_gap_coefficient=float(gc1.cut(slope*gc1.root3(F(4), 6)/4)))
    checks['dictionary: gap(A) = gap(H)/4 at theta = 4 theta_YM ; the weak coefficient is 0.1312 x 4^(1/3)/4 > 0.052'] = \
        slope*gc1.root3(F(4), 6)/4 > F(52, 1000)

    out = dict(stage='OL1', checks=len(checks), passed=sum(1 for v in checks.values() if v),
               all_pass=all(checks.values()), detail=checks, numbers=num)
    return out, num


if __name__ == '__main__':
    out, num = run()
    for k, v in out['detail'].items():
        print('PASS' if v else 'FAIL', k)
    print(json.dumps(num, indent=1))
    print(out['passed'], '/', out['checks'])
    with open(os.path.join(HERE, 'OL1_RESULT.json'), 'w') as f:
        json.dump(out, f, indent=1)
