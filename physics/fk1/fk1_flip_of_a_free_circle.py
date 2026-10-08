"""FK1: the flip of a free circle.

ONE_LAW: S = R + D, F = R - D.  For a free circle of a centre, seen R = (count factor)^2 and lost D = 1 - R.
AC1-A3: R = D at the last stable circle, in three cuts.  Here: for EVERY circle F is the squared ratio
(in-out rate / round rate), exactly in three cuts; what replaces it for a turning centre (SW2's form, QP1's rates).
sympy for the algebra, mpmath for numbers.  Python 3.12."""
import json, os, sys
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'qp1'))
import qp1_three_turns_and_the_rational_point as qp


def turning_count_factor():
    """(dtau/dt)^2 of the free circle in the plane of the equator on SW2's coframe (as built in QP1). Units r_s = 2."""
    w, a, th = sp.symbols('w a theta', positive=True)
    r, rs = w**2, 2
    rho, s = sp.sqrt(r**2 + a**2*sp.cos(th)**2), sp.sqrt(r**2 + a**2)
    b = sp.sqrt(rs*r)/rho
    TH = sp.Matrix([[1, 0, 0, 0], [b, rho/s, 0, -a*sp.sin(th)**2*b], [0, 0, rho, 0], [0, 0, 0, s*sp.sin(th)]])
    g = (TH.T*sp.diag(1, -1, -1, -1)*TH).applyfunc(sp.simplify)
    Om = 1/(w**3 + a)
    x0 = sp.Matrix([1, 0, 0, Om])
    return sp.simplify((x0.T*g*x0)[0].subs(th, sp.pi/2)), w, a


def run():
    out, num = {}, {}
    z = lambda e: sp.simplify(e) == 0

    # F1: centre at rest in d cuts.  MO1 form, fall^2 = m = c / r^(d-2).  Circle: Omega^2 = -(1/2r) d(1 - m)/dr... ; count factor^2 = 1 - m - r^2 Omega^2
    r, c = sp.symbols('r c', positive=True)
    ok_cf, ok_id = True, True
    for d in (3, 4, 5, 6, 7):
        m = c/r**(d - 2)
        f = 1 - m
        Om2 = sp.diff(f, r)/(2*r)
        cf2 = f - r**2*Om2
        kap2 = sp.Rational(1, 2)*f*sp.diff(f, r, 2) + 3*f*sp.diff(f, r)/(2*r) - sp.diff(f, r)**2
        ok_cf &= z(cf2 - (1 - sp.Rational(d, 2)*m)) and z(kap2/Om2 - ((4 - d) - d*m))          # CD1-D1 recovered
        ok_id &= z(kap2/Om2 - (2*cf2 - 1) - (3 - d))
    out['F1_count_factor_in_d_cuts'] = ok_cf
    out['F1_flip_is_rate_ratio_squared_up_to_3_minus_d'] = ok_id

    # F2: so in three cuts, with seen R = cf^2, lost D = 1 - R, q = in-out / round:  R - D = q^2 and R^2 = D^2 + q^2
    R, q = sp.symbols('R q', positive=True)
    D = 1 - R
    sol = sp.solve(sp.Eq(R - D, q**2), R)[0]
    out['F2_right_triangle'] = z(sol**2 - (1 - sol)**2 - q**2)
    ok_tri = True
    for (p_, q_) in [(1, 2), (1, 3), (2, 3), (1, 4), (3, 4), (2, 5), (3, 5)]:
        seen, lost = sp.Rational(q_**2 + p_**2, 2*q_**2), sp.Rational(q_**2 - p_**2, 2*q_**2)
        x = sp.Rational(q_**2 - p_**2, 3*q_**2)                                           # RC1: r_s / r at the rung p/q
        ok_tri &= (1 - sp.Rational(3, 2)*x == seen) and (q_**2 - p_**2)**2 + (2*p_*q_)**2 == (q_**2 + p_**2)**2 and seen - lost == sp.Rational(p_, q_)**2
    out['F2_rungs_are_pythagorean_triples'] = ok_tri
    out['F2_rung_one_third_is_3_4_5'] = (sp.Rational(5, 9), sp.Rational(4, 9), sp.Rational(1, 3)) == (sp.Rational(10, 18), sp.Rational(8, 18), sp.Rational(3, 9))

    # F3: turning centre
    cf2, w, a = turning_count_factor()
    y = a/w**3
    out['F3_count_factor_of_the_turning_centre'] = z(cf2 - (1 - 3/w**2 + 2*y)/(1 + y)**2)
    Om = 1/(w**3 + a)
    k2 = 1 - 6/w**2 + 8*a/w**3 - 3*a**2/w**4          # (in-out / round)^2, QP1-Q3
    n2 = 1 - 4*a/w**3 + 3*a**2/w**4                   # (up-down / round)^2, QP1-Q2
    out['F3_mean_of_the_two_cross_rates'] = z((k2 + n2)/2 - cf2/(1 - a*Om)**2)
    out['F3_reduces_to_F1_at_rest'] = z((k2 - (2*cf2 - 1)).subs(a, 0))
    # the last stable circle is no longer at seen = lost: exact value against the turn
    rr, aa = sp.symbols('r a', real=True)
    Gk = 1 - 6/rr + 8*aa/rr**sp.Rational(3, 2) - 3*aa**2/rr**2
    cfr = (1 - 3/rr + 2*aa/rr**sp.Rational(3, 2))/(1 + aa/rr**sp.Rational(3, 2))**2
    drda = (-sp.diff(Gk, aa)/sp.diff(Gk, rr)).subs({rr: 6, aa: 0})
    slope = sp.simplify((sp.diff(cfr, rr)*drda + sp.diff(cfr, aa)).subs({rr: 6, aa: 0}))
    num['slope_of_seen_at_the_last_stable_circle'] = str(slope)
    out['F3_slope'] = z(slope + sp.sqrt(6)/12)
    out['F3_against_the_turn_exact'] = z(Gk.subs({rr: 9, aa: -1})) and z(cfr.subs({rr: 9, aa: -1}) - sp.Rational(108, 169))
    tab = {}
    for av in (-1, -0.5, 0, 0.29, 0.5, 0.9, 0.99):
        rv = mp.findroot(lambda t: 1 - 6/t + 8*av/t**1.5 - 3*av**2/t**2, 6 - 3.2*av if av < 0.8 else (2.3 if av < 0.95 else 1.45))
        tab[str(av)] = dict(radius=float(rv), seen=float((1 - 3/rv + 2*av/rv**1.5)/(1 + av/rv**1.5)**2))
    num['last_stable_circle'] = tab
    out['F3_seen_is_half_only_at_rest'] = abs(tab['0']['seen'] - 0.5) < 1e-12 and all(abs(v['seen'] - 0.5) > 0.02 for k, v in tab.items() if k != '0')
    out['F3_seen_falls_with_the_turn'] = all(tab[str(x)]['seen'] > tab[str(y_)]['seen'] for x, y_ in zip((-1, -0.5, 0, 0.29, 0.5, 0.9), (-0.5, 0, 0.29, 0.5, 0.9, 0.99)))

    # F4: measured circle (GRO J1655-40: 441 +- 2, 298 +- 4, 17.3 +- 0.1 Hz, arXiv:1408.0884 Table 1)
    def seen_of(nphi, nper, nnod):
        M, av, rv = qp.fit3(nphi, nper, nnod)
        yv = av/rv**1.5
        direct = ((nphi - nper)**2 + (nphi - nnod)**2)/(2*nphi**2)          # (k2 + n2)/2, from the rates alone
        return [(1 - 3/rv + 2*yv)/(1 + yv)**2, direct, ((nphi - nper)/nphi)**2, 1 - 3/rv + 2*yv]
    base, err = qp.spread(seen_of, (441, 298, 17.3), (2, 4, 0.1))
    num['GRO_J1655_seen'] = dict(value=base[0], sigma=err[0])
    num['GRO_J1655_mean_cross_rate'] = dict(value=base[1], sigma=err[1])
    num['GRO_J1655_in_out_ratio_squared'] = dict(value=base[2], sigma=err[2])
    out['F4_rates_alone_give_the_mean_cross_rate'] = abs(base[1] - base[3]) < 1e-9

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'FK1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k, v in o.items(): print(k, v)
    for k, v in n.items(): print(k, v)
