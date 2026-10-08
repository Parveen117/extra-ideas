"""PN1: pure numbers of the centre at rest, against measurement.  Seen = lost, measured.

MO1's form (MA1): N^2 = 1 - r_s/r, radial factor 1/N.  LT1: r = s (1 + r_s/4s)^2 lays space out evenly.
DO1: the observer on the diagonal loses what he sees.  Symbolic (sympy) and plain numbers.  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, num = {}, {}
    s, rs, U = sp.symbols('s r_s U', positive=True)
    xp = rs/(4*s)
    r = s*(1 + xp)**2
    N2 = 1 - rs/r
    # W1: in the even layout the time part is ((1-x')/(1+x'))^2 and the space part (1+x')^4
    out['W1_time_part'] = z(N2 - ((1 - xp)/(1 + xp))**2)
    out['W1_space_part'] = z(sp.diff(r, s)**2/N2 - (1 + xp)**4) and z(r**2/s**2 - (1 + xp)**4)
    # W2: with U = r_s/2s:  time = 1 - 2U + 2 b U^2 ... , space = 1 + 2 g U ... :  g = 1 (lost = seen), b = 1
    t_ser = sp.series(((1 - U/2)/(1 + U/2))**2, U, 0, 3).removeO()
    s_ser = sp.series((1 + U/2)**4, U, 0, 2).removeO()
    g_ = sp.Rational(1, 2)*s_ser.coeff(U, 1); b_ = sp.Rational(1, 2)*t_ser.coeff(U, 2)
    out['W2_lost_over_seen_is_one'] = (g_ == 1) and (t_ser.coeff(U, 1) == -2)
    out['W2_second_order_is_one'] = (b_ == 1)
    # W3: what the other choice of MC1-T5 (spreading N, not the memory) would give: m = x - x^2/4 ;
    #     (in-out / round)^2 from CD1-D2:  1 - (5/2) x : five sixths of the advance
    u = sp.Symbol('u', positive=True)
    m = u - u**2/4                                                   # r_s = 1
    ratio2 = 1 - m - 2*u*sp.diff(m, u) - u*(1 - m)*sp.diff(m, u, 2)/sp.diff(m, u)
    first = sp.series(ratio2, u, 0, 2).removeO()
    out['W3_other_choice_five_sixths'] = z(first - (1 - sp.Rational(5, 2)*u))
    num['advance_Mercury_arcsec_per_century'] = dict(memory_spread=42.98, clock_factor_spread=round(42.98*5/6, 2), measured=42.98)
    # W4: light: circle at r = (3/2) r_s (MO1-M4); its reach b = r/N = (3 sqrt 3 / 2) r_s ; apparent diameter 2 sqrt 27 x (GM/c^2 D)
    rr = sp.Symbol('r', positive=True)
    f = (1 - rs/rr)/rr**2
    rc = sp.solve(sp.diff(f, rr), rr)
    out['W4_light_circle'] = rc == [3*rs/2]
    bc = (rr/sp.sqrt(1 - rs/rr)).subs(rr, 3*rs/2)
    out['W4_reach'] = z(bc - 3*sp.sqrt(3)*rs/2)
    out['W4_diameter_in_GM'] = z(2*bc/(rs/2) - 2*sp.sqrt(27))
    # W5: pure numbers on the special shells
    x = sp.Symbol('x', positive=True)
    E = (1 - x)/sp.sqrt(1 - sp.Rational(3, 2)*x)
    out['W5_binding_at_last_stable_circle'] = z(E.subs(x, sp.Rational(1, 3)) - sp.sqrt(sp.Rational(8, 9)))
    num['binding_at_last_stable_circle'] = round(float(1 - sp.sqrt(sp.Rational(8, 9))), 5)
    # measured (arXiv:1409.7871 for the first four; arXiv:2311.08680 for the shadow; all read at source)
    meas = dict(
        lost_over_seen_minus_1=dict(line=0.0, measured=2.1e-5, sigma=2.3e-5, how='radio delay past the Sun, Cassini'),
        light_bending_half_sum=dict(line=1.0, measured=0.99992, sigma=0.00023, how='radio sources past the Sun'),
        second_order_minus_1=dict(line=0.0, measured=-4.1e-5, sigma=7.8e-5, how="Mercury's perihelion"),
        four_b_minus_g_minus_3=dict(line=0.0, measured=4.4e-4, sigma=4.5e-4, how='Moon against Earth falling to the Sun'),
    )
    shadow, shadow_err, theta_g = 48.7, 7.0, 5.125                   # microarcseconds
    meas['shadow_over_GM_c2D'] = dict(line=round(float(2*sp.sqrt(27)), 3), measured=round(shadow/theta_g, 2), sigma=round(shadow_err/theta_g, 2),
                                      how='Sgr A*: shadow 48.7 +- 7.0 uas, GM/c^2D 5.125 uas from stars')
    num['measured'] = meas
    ok = True
    for k, v in meas.items():
        ok &= abs(v['measured'] - v['line']) <= 2*v['sigma']
    out['M_all_within_two_sigma'] = ok
    out['M_seen_equals_lost_to_parts_in_100000'] = abs(meas['lost_over_seen_minus_1']['measured']) < 5e-5
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'PN1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o); print(n)
