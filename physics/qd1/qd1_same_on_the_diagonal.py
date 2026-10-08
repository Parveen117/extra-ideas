"""QD1: on the diagonal the reading rule and the fixed-observer rule are the same.

Owner's statement: on the diagonal quantum and classical are the same.
Framework side: a reading read by a cut has seen R and lost D with S = R + D, F = R - D (T24-6.1, IN1, SD1);
the observer is on the diagonal R = D (DO1).  "Classical" here is the rule with the observer fixed outside:
a mark set beforehand and read by sign.  Exact counting, sympy, mpmath.  Python 3.12."""
import itertools, json, math, os
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def fixed_mark_correlation(Q, k):
    """Q marks at half-odd places of a Q-mark clock; a cut at place j reads the sign of cos(mark - cut). Exact count."""
    def s(m, j):
        d = (2*(m - j) + 1) % (2*Q)                 # angle (m + 1/2 - j) 2 pi / Q in units of pi / Q
        return 1 if (d < Q//2 or d > 2*Q - Q//2) else -1
    return Fr(sum(s(m, 0)*s(m, k) for m in range(Q)), Q)


def run():
    out, num = {}, {}
    z = lambda e: sp.simplify(e) == 0
    th = sp.symbols('theta', real=True)

    # D1: one reading, one cut (IN1, SD1): seen cos^2(th/2), lost sin^2(th/2), flip cos th, spread 4 R D
    R, D = sp.cos(th/2)**2, sp.sin(th/2)**2
    out['D1_seen_plus_lost_is_one_flip_is_cos'] = z(R + D - 1) and z(R - D - sp.cos(th))
    out['D1_spread_is_four_seen_lost'] = z((1 - (R - D)**2 - 4*R*D).rewrite(sp.exp))
    out['D1_diagonal_is_a_fair_coin_with_largest_spread'] = R.subs(th, sp.pi/2) == sp.Rational(1, 2) and (4*R*D).subs(th, sp.pi/2) == 1 and z(sp.diff(4*R*D, th).subs(th, sp.pi/2))

    # D2: pairs.  Fixed-mark rule by exact count: 1 - 4k/Q; reading rule: cos(2 pi k / Q)
    ok_lin, ok_same, ok_sign = True, True, True
    for Q in (8, 16, 32, 64):
        for k in range(Q//2 + 1):
            e = fixed_mark_correlation(Q, k)
            ok_lin &= e == 1 - Fr(4*k, Q)
            diff = sp.cos(2*sp.pi*sp.Rational(k, Q)) - sp.Rational(e.numerator, e.denominator)
            if k in (0, Q//4, Q//2):
                ok_same &= sp.simplify(diff) == 0
            else:
                ok_sign &= (float(diff) > 0) == (k < Q//4)
    out['D2_fixed_mark_rule_is_linear_in_the_angle'] = ok_lin
    out['D2_rules_agree_on_the_cut_and_on_the_diagonal'] = ok_same
    out['D2_and_nowhere_else'] = ok_sign
    # any rule odd under seen <-> lost (theta -> pi - theta) is zero on the diagonal
    f = sp.Function('f')
    out['D2_exchange_odd_rules_vanish_on_the_diagonal'] = sp.solve(sp.Eq(f(sp.pi/2), -f(sp.pi - sp.pi/2)), f(sp.pi/2)) == [0]
    # sums over four pairs: every fixed assignment gives +-2 (enumeration); four-mark clock: both rules 2; eight-mark: 4 x 1/sqrt2 against 4 x 1/2
    out['D2_fixed_assignments_give_two'] = all(abs(a*b - a*b2 + a2*b + a2*b2) == 2 for a, a2, b, b2 in itertools.product((1, -1), repeat=4))
    lin = lambda t: 1 - 2*t/sp.pi
    four = lambda E: E(sp.pi/2) - E(sp.pi) + E(0) + E(sp.pi/2)
    eight = lambda E: 3*E(sp.pi/4) - E(3*sp.pi/4)
    out['D2_four_mark_clock_same'] = four(sp.cos) == 2 and four(lin) == 2
    out['D2_eight_mark_clock'] = z(eight(sp.cos) - 2*sp.sqrt(2)) and eight(lin) == 2 and sp.cos(sp.pi/4) == 1/sp.sqrt(2) and lin(sp.pi/4) == sp.Rational(1, 2)
    meas, sig = 2.82759, 0.00051                                    # arXiv:1506.01865
    num['pair_sum'] = dict(reading_rule=2*math.sqrt(2), fixed_observer=2, measured=meas, sigma=sig,
                           from_reading_rule_in_sigma=(meas - 2*math.sqrt(2))/sig, above_fixed_observer_in_sigma=(meas - 2)/sig)
    out['D2_measured_sum_is_the_reading_rule'] = abs(meas - 2*math.sqrt(2)) < 2*sig and (meas - 2)/sig > 1000

    # D3: counts.  n = 1/(e^u - eta): spread = n (1 + eta n).  With R = n, D = n^2: S = R + D, R, F = R - D
    u, x, n = sp.symbols('u x n', positive=True)
    ok = True
    for eta in (1, 0, -1):
        nn = 1/(sp.exp(u) - eta)
        ok &= z(-sp.diff(nn, u) - nn*(1 + eta*nn))
    out['D3_spread_of_a_count'] = ok
    sg = sp.symbols('s', positive=True)
    Gf = (1 - x)/(1 - x*sg)                                         # generating form of the repeating record p_k = (1 - x) x^k
    mean_g = sp.diff(Gf, sg).subs(sg, 1); var_g = sp.diff(Gf, sg, 2).subs(sg, 1) + mean_g - mean_g**2
    out['D3_repeating_record_is_S'] = z(mean_g - x/(1 - x)) and z(var_g - (mean_g + mean_g**2))
    out['D3_two_valued_record_is_F'] = z((n - n**2) - n*(1 - n)) and sp.solve(sp.diff(n - n**2, n), n) == [sp.Rational(1, 2)] and (n - n**2).subs(n, sp.Rational(1, 2)) == sp.Rational(1, 4)

    # D4: the diagonal mode R = D: n = 1, x = 1/2, u = ln 2
    out['D4_diagonal_mode'] = sp.solve(sp.Eq(n, n**2), n) == [1] and sp.solve(sp.Eq(x/(1 - x), 1), x) == [sp.Rational(1, 2)] and z(sp.log(1/sp.Rational(1, 2)) - sp.log(2))
    ent = (1 + n)*sp.log(1 + n) - n*sp.log(n)
    U_, S_, F_ = sp.log(2), ent.subs(n, 1), sp.log(1 - sp.Rational(1, 2))
    out['D4_one_bit_two_bits_minus_one_bit'] = z(S_ - 2*sp.log(2)) and z(F_ + sp.log(2)) and z(U_ - S_ - F_) and (1 - x).subs(x, sp.Rational(1, 2)) == sp.Rational(1, 2)
    # curvature of the entropy: count law -1/s'' = n, wave law n^2, the actual law their sum (equal parts at n = 1)
    s_count, s_wave = n - n*sp.log(n), sp.log(n)
    out['D4_entropy_curvatures_add'] = z(-1/sp.diff(s_count, n, 2) - n) and z(-1/sp.diff(s_wave, n, 2) - n**2) and z(-1/sp.diff(ent, n, 2) - (n + n**2))

    # D5: numbers.  The sky's thermal light (T = 2.7255 K, recalled): its diagonal mode; shares on the wave side (u < ln 2)
    mp.mp.dps = 25
    L = mp.log(2)
    cnt = mp.quad(lambda t: t**2/mp.expm1(t), [0, L]); en = mp.quad(lambda t: t**3/mp.expm1(t), [0, L])
    K = 2000
    cnt_s = mp.nsum(lambda j: (2 - mp.e**(-j*L)*(j**2*L**2 + 2*j*L + 2))/j**3, [1, K])
    en_s = mp.nsum(lambda j: (6 - mp.e**(-j*L)*(j**3*L**3 + 3*j**2*L**2 + 6*j*L + 6))/j**4, [1, K])
    out['D5_series_with_explicit_tail'] = 0 <= cnt - cnt_s < mp.mpf(1)/K**2 and 0 <= en - en_s < mp.mpf(2)/K**3
    num['share_of_count_on_the_wave_side'] = float(cnt/(2*mp.zeta(3)))
    num['share_of_energy_on_the_wave_side'] = float(en/(mp.pi**4/15))
    kB, h = 1.380649e-23, 6.62607015e-34
    num['diagonal_mode_of_the_sky_GHz'] = kB*2.7255*math.log(2)/h/1e9
    num['quantum_over_kT_on_the_diagonal'] = math.log(2)

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'QD1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k, v in o.items(): print(k, v)
    for k, v in n.items(): print(k, v)
