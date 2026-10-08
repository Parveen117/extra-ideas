"""GM1: a charge on the displaced centre.  One displacement iota*a gives the turning of the centre and its magnetic moment.

Reads the contracted curvature of SW2's frame (sw2/sw2_contracted_curvature.json, written and checked by SW2):
frame e^1 = (rho/s) dr + b nu, b^2 = P(r)/rho^2, with P free.  Here the source is a light field.
Symbolic (sympy).  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))

r, x, a, rs, q2 = sp.symbols('r x a r_s q2', positive=True)
rho, s, sig, K = sp.symbols('rho s sigma K', positive=True)
P, P1, P2, P3 = sp.symbols('P P1 P2 P3')
SQ = [(rho, r**2 + a**2*x**2), (s, r**2 + a**2), (sig, 1 - x**2)]
ETA = sp.diag(1, -1, -1, -1)


def red(e):
    num, den = sp.fraction(sp.cancel(e))
    def lower(p):
        p = sp.expand(p)
        for sym, sq in SQ:
            pp = sp.Poly(p, sym)
            p = sp.expand(sum(co*sq**(k//2)*sym**(k % 2) for (k,), co in pp.terms()))
        return p
    return sp.cancel(lower(num)/lower(den))


def load():
    raw = json.load(open(os.path.join(HERE, '..', 'sw2', 'sw2_contracted_curvature.json')))
    ns = {str(v): v for v in (r, x, a, rho, s, sig, K, P, P1, P2, P3)}
    Ric = sp.zeros(4)
    for key, val in raw.items():
        j, l = int(key[0]), int(key[1])
        Ric[j, l] = Ric[l, j] = sp.sympify(val, locals=ns)
    return Ric


def run():
    out = {}
    z = lambda e: sp.simplify(e) == 0
    Ric = load()
    rho2, s2 = r**2 + a**2*x**2, r**2 + a**2
    # C1: the invariant part (trace) of the contracted curvature is -P''/rho^2 for every a
    trace = red(sum(ETA[i, i]*Ric[i, i] for i in range(4)))
    out['C1_invariant_part'] = red(trace + P2/rho2) == 0
    # so a source with no invariant part (a light field, EG1) leaves P = r_s r - q2 : two constants and no more
    ch = {P: rs*r - q2, P1: rs, P2: 0, P3: 0}
    Rq = Ric.subs(ch).applyfunc(red)
    u = q2/rho2**2
    v = a*sig/s
    # C2: components in the falling frame
    out['C2_components'] = (red(Rq[0, 0] - u*(1 + v**2)/(1 - v**2)) == 0 and red(Rq[3, 3] - u*(1 + v**2)/(1 - v**2)) == 0
                            and red(Rq[0, 3] + 2*u*v/(1 - v**2)) == 0 and red(Rq[1, 1] + u) == 0 and red(Rq[2, 2] - u) == 0
                            and all(Rq[j, l] == 0 for (j, l) in ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3))))
    # C3: in the frame that goes round the axis with the centre at speed v = a sin(theta)/sqrt(r^2+a^2) the pattern is EG1's
    g2 = 1/(1 - v**2)                                               # gamma^2
    L00 = g2*(Rq[0, 0] + 2*v*Rq[0, 3] + v**2*Rq[3, 3])
    L33 = g2*(Rq[3, 3] + 2*v*Rq[0, 3] + v**2*Rq[0, 0])
    L03 = g2*((1 + v**2)*Rq[0, 3] + v*(Rq[0, 0] + Rq[3, 3]))
    out['C3_pattern_of_EG1'] = red(L00 - u) == 0 and red(L33 - u) == 0 and red(L03) == 0
    out['C3_speed_below_one'] = red(1 - v**2 - rho2/s2) == 0
    # C4: u is the frame-independent energy of the light field F = -grad(q/R) of the displaced charge (EM1-T1: (1/2)|F.F|; SW2-F3)
    R = r - sp.I*a*x
    out['C4_modulus'] = z(sp.Abs(R**2)**2 - rho2**2) and z(sp.expand(R*sp.conjugate(R)) - rho2)
    out['C4_memory'] = z(sp.re(rs/R) - q2/(R*sp.conjugate(R)) - (rs*r - q2)/rho2)
    # C5: far away: q/R = q/d + iota q a cos/d^2 + ... : charge q and, in the turn part, a dipole of moment q a
    d, c_, e_ = sp.symbols('d c epsilon', positive=True)
    far = sp.series(1/sp.sqrt(1 - 2*c_*sp.I*e_ - e_**2), e_, 0, 2).removeO()
    out['C5_far_field'] = z(far - (1 + sp.I*e_*c_))
    # C5b: the field F = -grad(q/R) meets OB1-U3 with no source away from the ring: div F = 0 and curl F = 0, both parts
    X, Y, Z = sp.symbols('X Y Z', real=True)
    pot = 1/sp.sqrt(X**2 + Y**2 + (Z - sp.I*a)**2)
    F = [-sp.diff(pot, w_) for w_ in (X, Y, Z)]
    div = sum(sp.diff(F[i], w_) for i, w_ in enumerate((X, Y, Z)))
    curl = [sp.diff(F[2], Y) - sp.diff(F[1], Z), sp.diff(F[0], Z) - sp.diff(F[2], X), sp.diff(F[1], X) - sp.diff(F[0], Y)]
    out['C5b_no_source_away_from_the_ring'] = z(div) and all(z(v_) for v_ in curl)
    # on the axis far away the field is q/Z^2 + 2 iota q a / Z^3 : the field of a charge q and of a dipole q a along the axis
    on_axis = sp.series(F[2].subs({X: 0, Y: 0}).subs(Z, 1/e_), e_, 0, 4).removeO()
    out['C5c_axis_field'] = z(on_axis - (e_**2 + 2*sp.I*a*e_**3))
    # C6: the same displacement in the memory: r_s/R = r_s/d + iota r_s a cos/d^2 : swirl constant r_s a (SW2-L1, SW1)
    #     moment / charge = a = swirl constant / r_s .  With r_s = 2GM/c^2, swirl constant = 2GJ/c^3, moment = q a c:
    G, M, Jt, c0, q = sp.symbols('G M J c q', positive=True)
    a_len = (2*G*Jt/c0**3)/(2*G*M/c0**2)
    moment = q*a_len*c0
    gfac = moment/(q/(2*M)*Jt)
    out['C6_ratio_is_two'] = sp.simplify(gfac - 2) == 0
    # C7: a light-like shell needs r^2 - r_s r + a^2 + q2 = 0 to have a root: a^2 + q2 <= r_s^2/4
    disc = sp.discriminant(r**2 - rs*r + a**2 + q2, r)
    out['C7_shell_condition'] = z(disc - (rs**2 - 4*(a**2 + q2)))
    # numbers: the electron (J = hbar/2)
    hbar, c, Gn, me, e, k_e = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.602176634e-19, 8.9875517923e9
    a_e = (hbar/2)/(me*c); rs_e = 2*Gn*me/c**2; rq = (Gn*k_e*e**2)**0.5/c**2
    num = dict(a_m=a_e, r_s_m=rs_e, sqrt_q2_m=rq, a_over_half_r_s=a_e/(rs_e/2), g_model=2.0, g_measured_recalled=2.00231930436,
               shortfall=2.00231930436/2 - 1)
    out = {k: bool(val) for k, val in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'GM1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    print(run())
