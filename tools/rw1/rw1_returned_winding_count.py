"""RW1: the returned winding count.

Source object: the winding equation  dw/dh = K - w^2/h  of the "uncut horizons" draft (its "emotional
potential" (dw/dh)/sqrt(h) and "regret operator" are readings of this equation).
Sources read before building: RKF F00-E (native factorial exponential, tail lemma 2.1), RKF theorum/75-T1
(ladder law f_(c+1/2)/f_c <= kappa/(2(2c+2)), termwise on the positive series), physics CT1 (tanh k_c = I2/I1;
it quotes the face expansion f_j = 2 I_(2j+1)(kappa)/kappa from YM-3/4), physics SN1/SN2 (spread/mean).

Carrier: two factorial series  E(u) E(v) = sum u^a v^b/(a! b!)  (forward count a, returned count b).
Sector nu = a - b >= 0.  Weights on the sector:  t_n = y^n / (n! (n+nu)!) ,  y = u v ,  n = b.
    Z_nu(y) = sum t_n ,   w = <n> ,   V = <n^2> - <n>^2 .

R1  <n (n + nu)> = y  exactly (termwise).  With  y dw/dy = V :   y dw/dy = y - nu w - w^2 .
    nu = 0, y = K h is the draft's equation.  A power series solution has w(0) (w(0) + nu) = 0; the one with
    w(0) = 0 is unique and is w = <n>; for nu >= 1 the other start w(0) = -nu admits no power series.
R2  The sector weights depend on u, v only through the product y = u v  (<a b> = y is R1 read on both counts).
R3  Ladder:  Z_nu' = Z_(nu+1) ,  Z_nu - y Z_(nu+2) = (nu+1) Z_(nu+1) ,  hence  y = w_nu (nu + 1 + w_(nu+1)) .
R4  Spread law, two exact series with coefficients >= 0  (S_N = C(2N+2nu, N+nu)/(N! (N+2nu)!)):
        (V - w/2) Z^2 = sum S_N  N (2nu+1) / (4 (2N+2nu-1))  y^N
        (w - V)   Z^2 = sum S_N  N (N-1)   / (2 (2N+2nu-1))  y^N
    hence  w/2 < V < w  for every y > 0 and every nu >= 0.  The constants 1/2 and 1 cannot be improved:
    V/w is below 0.51 at y = 2500 and above 0.99 at y = 1/100 (certified points).
    The first coefficient is below S_N (2nu+1)/8 for nu >= 1 (below S_N/4 for nu = 0), so the remainder is bounded:
        0 < V - w/2 < (2nu+1)/8   (nu >= 1) ,      0 < V - w/2 < 1/4   (nu = 0) .
R5  Root-free enclosure of the ladder ratio:   w^2 + (nu + 1/2) w  <  y  <  w^2 + (nu + 1) w .
    With x = 2 sqrt(y), r = I_(nu+1)(x)/I_nu(x) = 2w/x :   2nu+1 < x (1 - r^2)/r < 2nu+2 .
R6  Face ladder (nu = 2c+1, x = kappa):  4c+3 < kappa (1 - r^2)/r < 4c+4  at every coupling;
    centre record (c = 0, r = tanh k_c):   sinh 2k_c = kappa/(1 + V/w) ,   kappa/2 < sinh 2k_c < 2 kappa/3 ,
    and for kappa > 2:   2 kappa/3 - kappa/(3 (kappa - 2)) < sinh 2k_c .
R7  K < 0: the weights alternate, there is no count; w leaves at the first zero of Z_0(-y), y0 in (36/25, 29/20)
    (Z_0(-y) >= 1 - y + y^2/4 - y^3/36 + y^4/576 - y^5/14400 > 0 on (0, 36/25], then a change of sign).
    On the regular branch the sign of dw/dh is the sign of K.

Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
from math import comb, factorial
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- formal series in y (lists of Fractions)
def coeff(nu, n):
    return F(1, factorial(n)*factorial(n + nu))


def zser(nu, order):
    return [coeff(nu, n) for n in range(order + 1)]


def smul(a, b, order):
    out = [F(0)]*(order + 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j > order:
                break
            out[i + j] += x*y
    return out


def sadd(a, b, s=1):
    n = max(len(a), len(b))
    a = a + [F(0)]*(n - len(a))
    b = b + [F(0)]*(n - len(b))
    return [x + s*y for x, y in zip(a, b)]


def sscale(a, c):
    return [c*x for x in a]


def sdiv(a, b, order):
    """a / b as a series, b[0] != 0"""
    out = []
    for n in range(order + 1):
        acc = a[n] if n < len(a) else F(0)
        for k in range(1, n + 1):
            if k < len(b):
                acc -= b[k]*out[n - k]
        out.append(acc/b[0])
    return out


def theta(a):
    """y d/dy"""
    return [n*x for n, x in enumerate(a)]


def shift(a, order):
    """multiplication by y"""
    return ([F(0)] + a)[:order + 1]


def mean_series(nu, order):
    z = zser(nu, order)
    return sdiv(theta(z), z, order)


def regular_branch(nu, order):
    """the power-series solution of  y w' = y - nu w - w^2  with w(0) = 0:
    (m + nu) a_m = [m = 1] - sum_(0<i<m) a_i a_(m-i)  fixes every coefficient"""
    a = [F(0)]*(order + 1)
    for m in range(1, order + 1):
        conv = sum(a[i]*a[m - i] for i in range(1, m))
        a[m] = ((1 if m == 1 else 0) - conv)/F(m + nu)
    return a


def other_start_obstruction(nu):
    """the start w(0) = -nu (nu >= 1):  (m - nu) a_m = [m = 1] - sum_(0<i<m) a_i a_(m-i).
    At m = nu the left side is zero; returns the right side there (non-zero = no power series)."""
    a = {0: F(-nu)}
    for m in range(1, nu + 1):
        rhs = (1 if m == 1 else 0) - sum(a[i]*a[m - i] for i in range(1, m))
        if m == nu:
            return rhs
        a[m] = rhs/F(m - nu)


def s_total(nu, n):
    return F(comb(2*n + 2*nu, n + nu), factorial(n)*factorial(n + 2*nu))


def spread_series(nu, order, lower=F(1, 2), upper=F(1)):
    """series of (V - lower*w) Z^2 and (upper*w - V) Z^2, using  w Z = theta Z  and  V = y - nu w - w^2"""
    z = zser(nu, order)
    tz = theta(z)
    yz2 = shift(smul(z, z, order), order)
    tz2 = smul(tz, tz, order)
    tzz = smul(tz, z, order)
    low = sadd(sadd(yz2, tz2, -1), sscale(tzz, nu + lower), -1)
    up = sadd(sadd(tz2, sscale(tzz, nu + upper)), yz2, -1)
    return low, up


# ---------------------------------------------------------------- certified values at a rational point
def enclose(nu, y, digits=40):
    """two-sided rational enclosures of Z_nu(y), theta Z_nu(y) for y > 0 (positive terms, geometric tail)"""
    eps = F(1, 10**digits)
    z = tz = F(0)
    n = 0
    while True:
        t = y**n*coeff(nu, n)
        z += t
        tz += n*t
        # from here on both term ratios are <= 1/2:  y/((n+1)(n+1+nu)) <= 1/2  and  y/(n(n+1+nu)) <= 1/2
        if n >= 1 and 2*y <= n*(n + 1 + nu) and (n + 1)*t <= eps:
            nxt = y**(n + 1)*coeff(nu, n + 1)
            return (z, z + 2*nxt), (tz, tz + 2*(n + 1)*nxt)
        n += 1


def mean_interval(nu, y, digits=40):
    (zl, zh), (tl, th) = enclose(nu, y, digits)
    return tl/zh, th/zl


def alternating_sign(y, terms=60):
    """sign of Z_0(-y) from two consecutive partial sums once the terms decrease (y < (n+1)^2)"""
    s, parts = F(0), []
    for n in range(terms):
        s += (-y)**n*coeff(0, n)
        parts.append(s)
    lo, hi = min(parts[-2:]), max(parts[-2:])
    return 1 if lo > 0 else (-1 if hi < 0 else 0)


def run():
    out, num = {}, {}
    order = 40
    nus = range(0, 7)

    # R1 -------------------------------------------------------------------------------------------------
    out['R1_second_moment_is_the_scale'] = all(
        n*(n + nu)*coeff(nu, n) == coeff(nu, n - 1) for nu in nus for n in range(1, 80))
    ok = True
    for nu in nus:
        w = mean_series(nu, order)
        lhs = theta(w)
        rhs = sadd(sadd(shift([F(1)], order), sscale(w, nu), -1), smul(w, w, order), -1)
        ok &= lhs[:order + 1] == rhs[:order + 1]
    out['R1_winding_equation_holds_in_every_sector'] = ok
    out['R1_series_with_w0_zero_is_the_mean'] = all(
        regular_branch(nu, order) == mean_series(nu, order) for nu in nus)
    obstruction = {nu: other_start_obstruction(nu) for nu in range(1, 7)}
    out['R1_the_other_start_has_no_power_series'] = all(v != 0 for v in obstruction.values())
    num['R1_obstruction_of_the_other_start'] = {str(k): str(v) for k, v in obstruction.items()}
    w0 = mean_series(0, 6)
    out['R1_first_terms_nu0'] = w0[:5] == [0, 1, F(-1, 2), F(1, 3), F(-11, 48)]
    # variance = theta(mean):  V Z^2 = (theta^2 Z) Z - (theta Z)^2
    ok = True
    for nu in nus:
        z = zser(nu, order)
        w = mean_series(nu, order)
        v_direct = sadd(sdiv(theta(theta(z)), z, order), smul(w, w, order), -1)
        ok &= v_direct == theta(w)
    out['R1_spread_is_the_scale_derivative_of_the_mean'] = ok

    # R2 -------------------------------------------------------------------------------------------------
    def sector(mu1, mu2, nu, nmax=25):
        ws = [mu1**(n + nu)*mu2**n/(factorial(n + nu)*factorial(n)) for n in range(nmax)]
        tot = sum(ws)
        return [x/tot for x in ws]
    out['R2_only_the_product_of_the_two_rates_enters'] = all(
        sector(F(3), F(1, 3), nu) == sector(F(1), F(1), nu) == sector(F(1, 7), F(7), nu) for nu in (0, 1, 4))

    # R3 -------------------------------------------------------------------------------------------------
    ok = True
    for nu in nus:
        z0, z1, z2 = zser(nu, order + 1), zser(nu + 1, order), zser(nu + 2, order)
        deriv = [(n + 1)*z0[n + 1] for n in range(order + 1)]
        ok &= deriv == z1
        ok &= sadd(z0[:order + 1], shift(z2, order), -1) == sscale(z1, nu + 1)
        wa, wb = mean_series(nu, order), mean_series(nu + 1, order)
        prod = smul(wa, sadd([F(nu + 1)], wb), order)
        ok &= prod == [F(0), F(1)] + [F(0)]*(order - 1)
    out['R3_ladder_identity'] = ok

    # R4 -------------------------------------------------------------------------------------------------
    big = 60
    ok_low = ok_up = ok_pos = True
    for nu in nus:
        low, up = spread_series(nu, big)
        for n in range(big + 1):
            s = s_total(nu, n)
            cl = s*F(n*(2*nu + 1), 4*(2*n + 2*nu - 1))
            cu = s*F(n*(n - 1), 2*(2*n + 2*nu - 1))
            ok_low &= low[n] == cl
            ok_up &= up[n] == cu
            ok_pos &= cl >= 0 and cu >= 0 and (n == 0 or cl > 0) and (n < 2 or cu > 0)
    out['R4_spread_minus_half_mean_is_a_positive_series'] = ok_low
    out['R4_mean_minus_spread_is_a_positive_series'] = ok_up
    out['R4_all_coefficients_non_negative'] = ok_pos
    # the count behind the coefficients: i successes among N+nu draws from N marked in 2N+2nu
    ok = True
    for nu in nus:
        for n in range(1, 30):
            ws = [comb(n, i)*comb(n + 2*nu, i + nu) for i in range(n + 1)]
            tot = sum(ws)
            m1 = F(sum(i*x for i, x in enumerate(ws)), tot)
            m2 = F(sum(i*i*x for i, x in enumerate(ws)), tot)
            ok &= tot == comb(2*n + 2*nu, n + nu) and m1 == F(n, 2)
            ok &= m2 - m1*m1 == F(n*(n + 2*nu), 4*(2*n + 2*nu - 1))
    out['R4_shell_spread_closed_form'] = ok
    # the constants are sharp: certified points where V/w is below 0.51 and above 0.99
    wl, wh = mean_interval(0, F(2500))
    out['R4_constant_one_half_is_sharp'] = (F(2500) - wl*wl)/wl < F(51, 100)
    wl, wh = mean_interval(0, F(1, 100))
    out['R4_constant_one_is_sharp'] = (F(1, 100) - wh*wh)/wh > F(99, 100)

    # bounded remainder:  V - w/2 < (2nu+1)/8  (nu >= 1),  < 1/4  (nu = 0)
    cap = lambda nu: F(1, 4) if nu == 0 else F(2*nu + 1, 8)
    out['R4_coefficients_below_the_cap'] = all(
        F(n*(2*nu + 1), 4*(2*n + 2*nu - 1)) <= cap(nu) and (nu == 0 or F(n*(2*nu + 1), 4*(2*n + 2*nu - 1)) < cap(nu))
        for nu in nus for n in range(1, 400))
    ok = True
    for nu in (0, 1, 2, 3, 5):
        for y in [F(1, 4), F(4), F(100), F(2500)]:
            wl, wh = mean_interval(nu, y)
            vh = y - wl*wl - nu*wl
            ok &= vh - wl/2 < cap(nu)
    out['R4_remainder_of_the_spread_is_bounded'] = ok
    wl, wh = mean_interval(1, F(2500))
    num['R4_remainder_at_nu1_y2500_against_three_eighths'] = round(float(F(2500) - wl*wl - wl - wl/2), 6)

    # R5 -------------------------------------------------------------------------------------------------
    grid_y = [F(1, 100), F(1, 4), F(1), F(4), F(25), F(100), F(2500)]
    ok, table = True, {}
    for nu in (0, 1, 2, 3, 5):
        for y in grid_y:
            wl, wh = mean_interval(nu, y)
            ok &= wh*wh + (nu + F(1, 2))*wh < y < wl*wl + (nu + 1)*wl
            vl, vh = y - wh*wh - nu*wh, y - wl*wl - nu*wl
            ok &= F(1, 2) < vl/wh and vh/wl < 1
            if nu in (0, 1):
                table[f'nu={nu},y={y}'] = round(float((vl/wh + vh/wl)/2), 6)
    out['R5_root_free_enclosure_on_the_grid'] = ok
    num['R5_spread_over_mean'] = table
    out['R5_spread_over_mean_runs_from_1_to_one_half'] = (
        table['nu=0,y=1/100'] > 0.99 and abs(table['nu=0,y=2500'] - 0.5) < 0.01)

    # R6 -------------------------------------------------------------------------------------------------
    ok_face = ok_centre = ok_t75 = ok_exact = ok_refined = True
    centre = {}
    for kappa in [F(1, 2), F(1), F(2), F(5), F(20), F(100)]:
        y = kappa*kappa/4
        for c2 in (0, 1, 2, 3):                                  # content c = c2/2,  nu = 2c + 1
            nu = c2 + 1
            wl, wh = mean_interval(nu, y)
            rl, rh = 2*wl/kappa, 2*wh/kappa
            ql, qh = kappa*(1 - rh*rh)/rh, kappa*(1 - rl*rl)/rl
            ok_face &= 2*c2 + 3 < ql and qh < 2*c2 + 4
            ok_t75 &= rh <= kappa/(2*(c2 + 2))                   # theorum/75-T1, recovered from R3
            if c2 == 0:
                sl, sh = 2*rl/(1 - rl*rl), 2*rh/(1 - rh*rh)      # sinh 2k_c with tanh k_c = r
                ok_centre &= kappa/2 < sl and sh < 2*kappa/3
                for wx in (wl, wh, F(1, 3)*kappa, F(2, 7)):       # an identity in (kappa, w), not only at the mean
                    rx, vx = 2*wx/kappa, y - wx*wx - wx
                    ok_exact &= 2*rx/(1 - rx*rx) == kappa/(1 + vx/wx)
                if kappa > 2:
                    ok_refined &= 2*kappa/3 - kappa/(3*(kappa - 2)) < sl
                centre[str(kappa)] = [round(float(kappa/2), 6), round(float((sl + sh)/2), 6),
                                      round(float(2*kappa/3), 6)]
    out['R6_face_ladder_two_sided_at_every_coupling'] = ok_face
    out['R6_ladder_law_of_theorum_75_recovered'] = ok_t75
    out['R6_centre_record_between_kappa_over_2_and_2kappa_over_3'] = ok_centre
    out['R6_centre_record_is_kappa_over_one_plus_spread_over_mean'] = ok_exact
    num['R6_centre_record_lower_value_upper'] = centre
    # refined lower bound at further couplings
    for kappa in [F(21, 10), F(3), F(7), F(50), F(400), F(3000)]:
        wl, wh = mean_interval(1, kappa*kappa/4)
        rl = 2*wl/kappa
        ok_refined &= 2*kappa/3 - kappa/(3*(kappa - 2)) < 2*rl/(1 - rl*rl)
        ok_refined &= kappa/2 - 1 < wl                           # w > sqrt(1 + y) - 1 > kappa/2 - 1
    out['R6_centre_record_within_a_bounded_distance_of_2kappa_over_3'] = ok_refined

    def sinh2k(kappa):                                           # enclosure of sinh 2k_c
        wl, wh = mean_interval(1, kappa*kappa/4)
        rl, rh = 2*wl/kappa, 2*wh/kappa
        return 2*rl/(1 - rl*rl), 2*rh/(1 - rh*rh)
    # CT1-C4's self-dual point sinh 2k_c = 1: below 1 at kappa = 3/2 and 1.886, above 1 at 1.887 and 2
    out['R6_self_dual_centre_record_bracketed'] = (
        sinh2k(F(3, 2))[1] < 1 and sinh2k(F(1886, 1000))[1] < 1 < sinh2k(F(1887, 1000))[0] and sinh2k(F(2))[0] > 1)

    # R7 -------------------------------------------------------------------------------------------------
    out['R7_sign_change_between_36_over_25_and_29_over_20'] = (
        alternating_sign(F(36, 25)) == 1 and alternating_sign(F(29, 20)) == -1)
    # no earlier zero: on (0, 36/25] the terms decrease from n = 1 (y < 4), so Z_0(-y) >= q(y), the partial sum
    # to n = 5; q' <= -1 + y/2 + y^3/144 < 0 there, so q >= q(36/25) > 0
    yy = F(36, 25)
    q = sum((-yy)**n*coeff(0, n) for n in range(6))
    out['R7_no_zero_before_36_over_25'] = q > 0 and -1 + yy/2 + yy**3/144 < 0
    wneg = [((-1)**m)*x for m, x in enumerate(mean_series(0, order))]     # the K < 0 branch, in |K| h
    out['R7_negative_K_branch_decreases'] = all(x < 0 for x in wneg[1:])

    # controls ----------------------------------------------------------------------------------------------
    one = [F(1, factorial(n)) for n in range(order + 1)]          # a single strand: weights y^n/n!
    w1 = sdiv(theta(one), one, order)
    rhs = sadd(shift([F(1)], order), smul(w1, w1, order), -1)
    out['control_single_strand_does_not_obey_the_equation'] = theta(w1) != rhs
    out['control_single_strand_has_spread_equal_to_mean'] = theta(w1) == w1
    bad = regular_branch(0, order)
    bad[2] += F(1, 1000)
    out['control_planted_coefficient_is_detected'] = bad != mean_series(0, order)

    out['pass'] = all(v for v in out.values() if isinstance(v, bool))
    return out, num


if __name__ == '__main__':
    out, num = run()
    with open(os.path.join(HERE, 'RW1_RESULT.json'), 'w') as fh:
        json.dump({'checks': out, 'numbers': num}, fh, indent=1, sort_keys=True)
        fh.write('\n')
    for k, v in out.items():
        print(f'{k}: {v}')
