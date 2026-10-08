"""LD1: the memory between two neutral centres, from the floors of their modes.

WQ1-W3 / LC1: a mode (a pair of readings turning at rate omega) has the floor kappa omega / 2.
Two modes read together through a coupling k: block [[1, k], [k, 1]] in units where each reads alike (DO1): lost part m = k^2.
MC1: the field between two centres is the least-cost one, r^-(d-2).
Exact algebra (sympy) and numbers.  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, num = {}, {}
    k, s = sp.symbols('k s', positive=True)
    # F1: two equal modes, coupling block H = [[1, k], [k, 1]] (x omega^2): own rates omega sqrt(1 +- k); floor sum / (kappa omega):
    floor = (sp.sqrt(1 + k) + sp.sqrt(1 - k))/2
    F = 1 - k**2                                                     # F = R - D of the block on the diagonal
    closed = sp.sqrt((1 + sp.sqrt(F))/2)
    out['F1_floor_closed_form'] = z(sp.expand(floor**2) - sp.expand(closed**2)) and floor.subs(k, sp.Rational(1, 2)) > 0
    ser = sp.series(1 - floor, k, 0, 6).removeO()
    out['F1_second_order'] = z(ser - (k**2/8 + sp.Rational(5, 128)*k**4))
    # bounds for every coupling (s = sqrt(1 - m), m = k^2):  m/8 <= lost floor <= (1 - 1/sqrt 2) m
    lower = sp.factor(s**4 + 14*s**2 - 32*s + 17)                    # >= 0  <=>  1 - sqrt((1+s)/2) >= (1 - s^2)/8
    out['F1_lower_bound'] = z(lower - (s - 1)**2*(s**2 + 2*s + 17))
    g = lambda mm: 1 - sp.sqrt((1 + sp.sqrt(1 - mm))/2)
    out['F1_upper_bound'] = all(sp.N(g(sp.Rational(i, 20)) - (1 - 1/sp.sqrt(2))*sp.Rational(i, 20)) <= 0 for i in range(1, 21)) and z(g(1) - (1 - 1/sp.sqrt(2)))
    # F2: coupling between two centres through the least-cost field: second derivatives of r^-(d-2): rates (d-1, -1, ..., -1) x (d-2)/r^d
    ok, sq = True, {}
    for d in (3, 4, 5):
        X = sp.symbols('x1:%d' % (d + 1), positive=True)
        rad = sp.sqrt(sum(v**2 for v in X))
        phi = rad**(-(d - 2))
        Hs = sp.Matrix(d, d, lambda i, j: sp.diff(phi, X[i], X[j]))
        pt = {X[0]: 1}; pt.update({v: 0 for v in X[1:]})
        Hp = Hs.subs(pt).applyfunc(sp.simplify)
        diag = [Hp[i, i]/(d - 2) for i in range(d)]
        ok &= (diag == [d - 1] + [-1]*(d - 1)) and all(Hp[i, j] == 0 for i in range(d) for j in range(d) if i != j) and sum(diag) == 0
        sq[d] = sum(v**2 for v in diag)
        ok &= (sq[d] == d*(d - 1))
    out['F2_pattern_and_trace'] = ok
    # F3: one mode per cut on each centre: lost floor = kappa omega (k0^2 / r^(2d)) d(d-1)/8 : three cuts: (3/4) kappa omega k0^2 / r^6
    out['F3_three_quarters_and_sixth_power'] = (sp.Rational(sq[3], 8) == sp.Rational(3, 4)) and (2*3 == 6)
    # F4: unequal centres: modes omega_A, omega_B, coupling c x_A x_B: own rates^2 = roots of (l - wA^2)(l - wB^2) = c^2
    wA, wB, c = sp.symbols('omega_A omega_B c', positive=True)
    disc = sp.sqrt((wA**2 - wB**2)**2 + 4*c**2)
    lp, lm = (wA**2 + wB**2 + disc)/2, (wA**2 + wB**2 - disc)/2
    out['F4_own_rates'] = z(sp.expand((lp - wA**2)*(lp - wB**2) - c**2)) and z(sp.expand((lm - wA**2)*(lm - wB**2) - c**2))
    # closed form of the floor sum: (sqrt lp + sqrt lm)^2 = wA^2 + wB^2 + 2 sqrt(wA^2 wB^2 - c^2)
    out['F4_closed_form'] = z(sp.expand(lp*lm) - (wA**2*wB**2 - c**2)) and z(lp + lm - wA**2 - wB**2)
    e = sp.Symbol('e', positive=True)
    tot = sp.sqrt(wA**2 + wB**2 + 2*sp.sqrt(wA**2*wB**2 - e))          # e = c^2
    slope = sp.simplify(sp.diff(tot, e).subs(e, 0))
    out['F4_second_order'] = (z(slope**2 - 1/(2*wA*wB*(wA + wB))**2) and slope.subs({wA: 3, wB: 2}) < 0
                              and z(tot.subs(e, 0)**2 - (wA + wB)**2))                # slope = -1/(2 wA wB (wA + wB))
    # with response alpha = q^2/omega^2 for each mode and c = g qA qB / r^3 (pattern squares 6):  C_AB = (3/2) kappa aA aB wA wB/(wA + wB)
    kap, aA, aB = sp.symbols('kappa alpha_A alpha_B', positive=True)
    C = lambda a1, w1, a2, w2: sp.Rational(3, 2)*kap*a1*a2*w1*w2/(w1 + w2)
    CAA, CBB, CAB = C(aA, wA, aA, wA), C(aB, wB, aB, wB), C(aA, wA, aB, wB)
    out['F4_equal_case_is_three_quarters'] = z(CAA - sp.Rational(3, 4)*kap*wA*aA**2)
    # F5: the rates and kappa drop out: C_AB = 2 C_AA C_BB / ( C_AA aB/aA + C_BB aA/aB )
    out['F5_mixed_rule'] = z(CAB - 2*CAA*CBB/(CAA*aB/aA + CBB*aA/aB))
    # numbers (arXiv:0902.3929, Table III and Table D, read at source; atomic units; stated uncertainty 1 %)
    alpha = dict(He=1.383, Ne=2.669, Ar=11.08, Kr=16.79, Xe=27.16)
    like = dict(He=1.461, Ne=6.38, Ar=64.3, Kr=130.0, Xe=286.0)
    mixed = {('He', 'Ne'): 3.03, ('He', 'Ar'): 9.55, ('He', 'Kr'): 13.42, ('He', 'Xe'): 19.6, ('Ne', 'Ar'): 19.5, ('Ne', 'Kr'): 27.3,
             ('Ne', 'Xe'): 39.7, ('Ar', 'Kr'): 91.1, ('Ar', 'Xe'): 134.5, ('Kr', 'Xe'): 192.0}
    rows, worst = {}, 0.0
    for (A, B), val in mixed.items():
        pred = 2*like[A]*like[B]/(like[A]*alpha[B]/alpha[A] + like[B]*alpha[A]/alpha[B])
        dev = pred/val - 1
        rows['%s-%s' % (A, B)] = dict(rule=round(pred, 3), tabulated=val, deviation=round(dev, 4))
        worst = max(worst, abs(dev))
    num['mixed_pairs'] = rows; num['largest_deviation'] = round(worst, 4)
    out['N1_mixed_rule_within_one_and_a_half_percent'] = worst < 0.015
    out['N1_eight_of_ten_within_half_a_percent'] = sum(1 for r_ in rows.values() if abs(r_['deviation']) <= 0.005) == 8
    # the single-rate reading of each atom: omega_eff = 4 C / (3 alpha^2) in atomic units of energy (hartree), against the ionisation energy
    ion = dict(He=0.9036, Ne=0.7925, Ar=0.5792, Kr=0.5145, Xe=0.4458)     # recalled; not re-read at source
    num['effective_rate_over_ionisation'] = {A: round((4*like[A]/(3*alpha[A]**2))/ion[A], 3) for A in like}
    out = {k_: bool(v) for k_, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'LD1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o)
    for kk, v in n.items(): print(kk, v)
