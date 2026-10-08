"""ET1: the heat reading of a turning (and charged) centre in the rim reader's coin; T = 0 as the quarter-turn clock.
Exact sympy."""
import json, os, sympy as sp
S, J, Q, x, v, chi, kap = sp.symbols('S J Q x v chi kappa', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
R = sp.sqrt(S/sp.pi)
E = sp.sqrt((R/2 + Q**2/(2*R))**2 + J**2/R**2)             # SP1-T5
T = sp.diff(E, S); T0 = 1/(4*sp.pi*R)                      # T0: neutral, non-turning centre of the same count
M0 = (R/2)*(1 + x)
sub = lambda e: sp.simplify(e.subs(Q, sp.sqrt(x)*R).subs(J, sp.sqrt(v**2*M0**2/(1 - v**2))*R))   # x = (Q/R)^2, v = rim speed
ratio = sub(T/T0)
target = (1 - x - 2*v**2)/sp.sqrt(1 - v**2)
out['T1_heat_reading'] = z(ratio**2 - target**2) and all(abs(sp.N((ratio - target).subs({x: a, v: b}))) < 1e-12 for a, b in ((sp.Rational(1, 5), sp.Rational(1, 3)), (sp.Rational(1, 2), sp.Rational(3, 5)), (sp.Rational(9, 10), sp.Rational(1, 10)), (sp.Rational(1, 10), sp.Rational(4, 5))))
out['T1_v_is_rim_speed'] = z(sub((J/R)/E) - v)
# T2: the line T = 0 in the two partners: Phi0^2 + 2 v^2 = 1  (Phi0 = Q/R)
out['T2_zero_line'] = sp.solve(sp.Eq(1 - x - 2*v**2, 0), x) == [1 - 2*v**2] and (1 - x - 2*v**2).subs({x: 1, v: 0}) == 0 and z((1 - x - 2*v**2).subs({x: 0, v: 1/sp.sqrt(2)}))
# T3: spin only, coin of the rim reader (N, beta) = (cos chi, sin chi) (UN2, RT1): T/T0 = cos 2chi / cos chi
r_spin = ratio.subs(x, 0).subs(v, sp.sin(chi))
unit = lambda a: sp.cos(a)                                  # HB1: unit of a clock with tick phase a
out['T3_ratio_of_units'] = all(abs(sp.N(r_spin.subs(chi, t) - unit(2*t)/unit(t))) < 1e-12 for t in (0.1, 0.4, 0.7, 0.78))
out['T3_zero_at_eighth_turn'] = unit(2*sp.pi/4) == 0 and sp.simplify(sp.sin(sp.pi/4) - 1/sp.sqrt(2)) == 0
# the tick and its square at the zero: zeta = exp(i pi/4) (eight marks), zeta^2 = i (four marks, HB1 unit 0)
zeta = sp.exp(sp.I*sp.pi/4)
out['T3_squared_tick_is_quarter_turn'] = sp.simplify(zeta**2 - sp.I) == 0 and sp.simplify(zeta**8 - 1) == 0
# T4: charge only, Phi0 = tan kappa: T/T0 = 1 - Phi0^2 = cos 2kappa / cos^2 kappa, zero at kappa = pi/4 as well
r_ch = ratio.subs(v, 0).subs(x, sp.tan(kap)**2)
out['T4_charge'] = z(r_ch - sp.cos(2*kap)/sp.cos(kap)**2) and r_ch.subs(kap, sp.pi/4) == 0
# T5: both: T/T0 = (cos 2chi - tan^2 kappa) / cos chi
both = ratio.subs({x: sp.tan(kap)**2, v: sp.sin(chi)})
out['T5_both'] = all(abs(sp.N((both - (sp.cos(2*chi) - sp.tan(kap)**2)/sp.cos(chi)).subs({chi: a, kap: b}))) < 1e-12 for a, b in ((0.2, 0.3), (0.5, 0.1), (0.1, 0.6)))
out = {k: bool(val) for k, val in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ET1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
