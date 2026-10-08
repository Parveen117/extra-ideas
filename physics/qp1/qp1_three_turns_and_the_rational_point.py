"""QP1: the three rates of a free circle on the exact frame of SW2, and the rational point 1/3 against measured rates.

Count form of SW2 (memory Re(r_s/R)); free readings take the largest count (CL2, MO1).  A free circle in the plane of the
equator has a round rate, an in-out rate and an up-down rate.  Units in the symbolic part: r_s = 2, c = 1.
Symbolic (sympy) for the rates; mpmath for the numbers.  Python 3.12."""
import json, os
import sympy as sp
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))


def rates_from_the_form():
    w, a, th = sp.symbols('w a theta', positive=True)            # r = w^2
    r, rs = w**2, 2
    rho, s = sp.sqrt(r**2 + a**2*sp.cos(th)**2), sp.sqrt(r**2 + a**2)
    b = sp.sqrt(rs*r)/rho
    TH = sp.Matrix([[1, 0, 0, 0], [b, rho/s, 0, -a*sp.sin(th)**2*b], [0, 0, rho, 0], [0, 0, 0, s*sp.sin(th)]])   # SW2's coframe
    g = (TH.T*sp.diag(1, -1, -1, -1)*TH).applyfunc(sp.simplify)
    Om, om, rr = sp.symbols('Omega omega r', positive=True)
    x0 = sp.Matrix([1, 0, 0, Om])                                # the circle, parameter t
    F = (x0.T*g*x0)[0]
    Fr = F.subs(w, sp.sqrt(rr))
    out = {}
    # the circle: d_r ( g_ab x0^a x0^b ) = 0
    circle = sp.simplify(sp.diff(Fr, rr).subs(th, sp.pi/2).subs(rr, w**2))
    Omv = 1/(w**3 + a)
    out['Q1_round_rate'] = sp.simplify(circle.subs(Om, Omv)) == 0
    # up-down: the theta equation decouples by the mirror symmetry of the plane
    nu_th2 = sp.simplify((-sp.Rational(1, 2)*sp.diff(F, th, 2).subs(th, sp.pi/2)/g[2, 2].subs(th, sp.pi/2)).subs(Om, Omv))
    out['Q2_up_down_rate'] = sp.simplify(nu_th2/Omv**2 - (1 - 4*a/w**3 + 3*a**2/w**4)) == 0
    # in-out: the linearised equations on (t, r, phi)
    ge = g.subs(th, sp.pi/2)
    gr = sp.Matrix(4, 4, lambda i, j: sp.diff(g[i, j].subs(w, sp.sqrt(rr)), rr).subs(rr, w**2).subs(th, sp.pi/2))
    F2 = sp.diff(Fr, rr, 2).subs(rr, w**2).subs(th, sp.pi/2)
    idx = [0, 1, 3]
    Mm = sp.zeros(3)
    for i_, mu in enumerate(idx):
        for j_, lam in enumerate(idx):
            val = -om**2*ge[mu, lam]
            A_ml = (gr[mu, 0] + gr[mu, 3]*Om) if lam == 1 else 0
            A_lm = (gr[lam, 0] + gr[lam, 3]*Om) if mu == 1 else 0
            val += sp.I*om*(A_ml - A_lm)
            if mu == 1 and lam == 1:
                val += -sp.Rational(1, 2)*F2
            Mm[i_, j_] = val
    det = sp.simplify(Mm.det().subs(Om, Omv))
    kap2 = Omv**2*(1 - 6/w**2 + 8*a/w**3 - 3*a**2/w**4)
    out['Q3_in_out_rate'] = sp.simplify(det.subs(om**2, kap2).subs(om, sp.sqrt(kap2))) == 0 and sp.simplify(det/om**4).subs(om, 0) != 0
    # a = 0: MO1-M5 (ratio^2 = 1 - 3 r_s/r) and no turn of the plane
    out['Q4_centre_at_rest'] = sp.simplify((kap2/Omv**2).subs(a, 0) - (1 - 3*rs/r)) == 0 and sp.simplify(nu_th2.subs(a, 0) - Omv.subs(a, 0)**2) == 0
    # first order in a: the plane of the circle turns at SW1's swirl rate W = r_s a / r^3
    plane = sp.series(Omv - sp.sqrt(nu_th2), a, 0, 2).removeO()
    out['Q5_plane_turns_at_the_swirl'] = sp.simplify(plane - rs*a/r**3) == 0
    return out


BETA = 299792458.0**3/(2*mp.pi*1.32712440018e20)                  # c^3 / (2 pi G M_sun), Hz


def rates(M, a, r):
    """M in solar masses, a and r in units of GM/c^2 (= r_s/2).  Returns round rate, leftover-turn rate, plane-turn rate, ratio."""
    nphi = BETA/M/(r**1.5 + a)
    kr = mp.sqrt(1 - 6/r + 8*a/r**1.5 - 3*a**2/r**2)
    kt = mp.sqrt(1 - 4*a/r**1.5 + 3*a**2/r**2)
    return nphi, nphi*(1 - kr), nphi*(1 - kt), kr


def fit3(nphi, nper, nnod):
    return mp.findroot(lambda M, a, r: [rates(M, a, r)[0] - nphi, rates(M, a, r)[1] - nper, rates(M, a, r)[2] - nnod], (5.3, 0.29, 5.7))


def fit_closed(nphi=None, nnod=None, nper=None, ratio=mp.mpf(1)/3, start=(5.3, 0.29, 5.7)):
    """the circle at the rational point: in-out rate / round rate = ratio; two measured rates give mass, a and radius."""
    def f(M, a, r):
        v = rates(M, a, r)
        eqs = [v[3] - ratio, v[2] - nnod]
        eqs.append(v[0] - nphi if nphi is not None else v[1] - nper)
        return eqs
    return mp.findroot(f, start)


def spread(fun, vals, errs):
    """first-order propagation by finite differences."""
    base = fun(*vals)
    var = [mp.mpf(0)]*len(base)
    for i, e in enumerate(errs):
        up = list(vals); up[i] = vals[i] + e
        d = fun(*up)
        var = [var[k] + (d[k] - base[k])**2 for k in range(len(base))]
    return [float(b) for b in base], [float(mp.sqrt(v)) for v in var]


def run():
    out = rates_from_the_form()
    num = {}
    # GRO J1655-40: 17.3 +- 0.1, 298 +- 4, 441 +- 2 Hz (arXiv:1408.0884, Table 1)
    f3, e3 = spread(lambda p, q, n: list(fit3(p, q, n)), (441, 298, 17.3), (2, 4, 0.1))
    num['GRO_three_rates_M_a_r'] = dict(value=[round(v, 4) for v in f3], sigma=[round(v, 4) for v in e3])
    ratio = 1 - 298/441
    ratio_err = ((4/441)**2 + (298*2/441**2)**2)**0.5
    num['GRO_measured_ratio'] = dict(value=round(ratio, 4), sigma=round(ratio_err, 4), sigmas_from_one_third=round((1/3 - ratio)/ratio_err, 2))
    fc, ec = spread(lambda p, n: list(fit_closed(nphi=p, nnod=n)), (441, 17.3), (2, 0.1))
    num['GRO_closed_M_a_r'] = dict(value=[round(v, 4) for v in fc], sigma=[round(v, 4) for v in ec])
    num['GRO_closed_lower_rate'] = dict(predicted=round(441*2/3, 1), sigma=round(2*2/3, 1), measured=298, measured_sigma=4)
    num['GRO_mass_from_light_curve'] = dict(value=5.4, sigma=0.3)
    # XTE J1550-564: 13.08 +- 0.08 and 183 +- 5 Hz; mass from the companion's motion 9.10 +- 0.61 (same table)
    fx, ex = spread(lambda q, n: list(fit_closed(nper=q, nnod=n, start=(9.0, 0.34, 5.5))), (183, 13.08), (5, 0.08))
    num['XTE_closed_M_a_r'] = dict(value=[round(v, 4) for v in fx], sigma=[round(v, 4) for v in ex])
    num['XTE_mass_from_companion'] = dict(value=9.10, sigma=0.61)
    # H 1743-322: 165 (+9 -5) and 240 +- 3 Hz
    num['H1743_measured_ratio'] = dict(value=round(1 - 165/240, 4), plus=round(165/240 - (165 - 5)/243, 4), minus=round((165 + 9)/237 - 165/240, 4))
    # the ladder for a centre at rest
    lad = {}
    for p, q in ((1, 2), (1, 3), (1, 4), (2, 5)):
        x = (1 - mp.mpf(p)**2/q**2)/3
        lad['%d/%d' % (p, q)] = dict(r_over_r_s=float(1/x), rates_upper_to_lower='%d:%d' % (q, q - p), round_rate_times_mass_Hz=round(float(BETA/(2/x)**1.5), 1))
    num['ladder_centre_at_rest'] = lad
    ok = dict(out)
    ok['N1_three_rate_fit_matches_source'] = abs(f3[0] - 5.31) < 0.07 and abs(f3[1] - 0.2875) < 0.006 and abs(f3[2] - 5.68) < 0.04
    ok['N2_pair_at_one_third_within_two_sigma'] = abs(1/3 - ratio) < 2*ratio_err
    ok['N3_GRO_closed_mass_within_light_curve_mass'] = abs(fc[0] - 5.4) < (0.3**2 + ec[0]**2)**0.5
    ok['N4_XTE_closed_mass_within_companion_mass'] = abs(fx[0] - 9.10) < (0.61**2 + ex[0]**2)**0.5
    ok = {k: bool(v) for k, v in ok.items()}
    ok['pass'] = all(ok.values())
    json.dump(dict(checks=ok, numbers=num), open(os.path.join(HERE, 'QP1_RESULT.json'), 'w'), indent=1)
    return ok, num


if __name__ == '__main__':
    o, n = run()
    print(o)
    for k, v in n.items():
        print(k, v)
