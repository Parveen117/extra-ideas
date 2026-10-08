"""ET2: (i) ET1's line against the usual condition for a charged turning centre; (ii) expansion as a third partner:
the line's two joining rules against the published mass formula, the correction needed, and the T = 0 surface.
Exact sympy."""
import json, os, sympy as sp
M, Q, a, S, J, h, R, F, v, x = sp.symbols('M Q a S J h R Phi0 v x', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
# ---- (i) usual variables at the T = 0 end: r_plus = M, a^2 + Q^2 = M^2, S = pi (r_plus^2 + a^2), J = a M
a2 = M**2 - Q**2
R2 = M**2 + a2                                              # R^2 = S/pi
Phi0sq = Q**2/R2; vsq = (a2*M**2/R2)/M**2                    # v = (J/R)/E with E = M
out['T1_line_is_the_usual_condition'] = z(Phi0sq + 2*vsq - 1)
# ---- (ii) expansion. Rest value by LC1/XU2: M0 = R/2 + Q^2/(2R) - h^2 R^3/2
M0 = R/2 + Q**2/(2*R) - h**2*R**3/2
naive = M0**2 + J**2/R**2                                    # SP1's rule unchanged
# published mass formula for a charged turning centre with a cosmological term (written with 1/l^2 -> -h^2, S = pi R^2)
Sx = sp.pi*R**2
known = Sx/(4*sp.pi) + sp.pi*(4*J**2 + Q**4)/(4*Sx) + Q**2/2 - J**2*h**2 - (Sx*h**2/(2*sp.pi))*(Q**2 + Sx/sp.pi - Sx**2*h**2/(2*sp.pi**2))
out['T2_rest_part_agrees'] = z(known.subs(J, 0) - M0**2)                 # charge + expansion: stored-cost rule holds
out['T2_naive_spin_rule_fails'] = not z(known - naive)
corrected = M0**2 + (J**2/R**2)*(1 - h**2*R**2)              # the rim's momentum weighed by the expansion's clock factor at R
out['T2_corrected_rule'] = z(known - corrected)
# ---- T3: heat reading with three partners. Phi0 = Q/R, x = h R, v^2 = (momentum term)/E^2
E2 = corrected
dE2 = sp.diff(E2, R)                                         # T = dE/dS = (dE2/dR)/(2E) / (2 pi R);  T0 = 1/(4 pi R)
subs = {Q: F*R, h: x/R}
M0s = sp.simplify(M0.subs(subs))                             # = (R/2)(1 + F^2 - x^2)
Es = M0s/sp.sqrt(1 - v**2)
Jsq = sp.solve(sp.Eq((J**2/R**2)*(1 - x**2), v**2*Es**2), J**2)[0]
ratio = sp.simplify((dE2.subs(subs).subs(J**2, Jsq)/(2*Es))*2)   # T/T0 = 2 dE/dR
target = ((1 - v**2)*(1 - F**2 - 3*x**2) - v**2*(1 + F**2 - x**2)/(1 - x**2))/sp.sqrt(1 - v**2)
pts = [(sp.Rational(1, 5), sp.Rational(1, 3), sp.Rational(1, 4)), (sp.Rational(1, 2), sp.Rational(1, 5), sp.Rational(1, 10)), (sp.Rational(1, 10), sp.Rational(3, 5), sp.Rational(2, 5))]
out['T3_heat_reading'] = all(abs(sp.N((ratio - target).subs({F: p[0], v: p[1], x: p[2]}))) < 1e-12 for p in pts)
out['T3_reduces_to_ET1'] = z(target.subs(x, 0) - (1 - F**2 - 2*v**2)/sp.sqrt(1 - v**2))
# ---- T4: the T = 0 surface and its three ends
surf = (1 - v**2)*(1 - x**2)*(1 - F**2 - 3*x**2) - v**2*(1 + F**2 - x**2)
out['T4_three_ends'] = surf.subs({v: 0, x: 0, F: 1}) == 0 and z(surf.subs({F: 0, x: 0, v: 1/sp.sqrt(2)})) and z(surf.subs({F: 0, v: 0, x: 1/sp.sqrt(3)}))
# the expansion end is XU1's state of unit 0: 27 h^2 M^2 = 1 at x = 1/sqrt 3
out['T4_expansion_end_is_unit_zero'] = z((27*h**2*(R*(1 - h**2*R**2)/2)**2).subs(h, 1/(sp.sqrt(3)*R)) - 1)
out = {k: bool(val) for k, val in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ET2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
