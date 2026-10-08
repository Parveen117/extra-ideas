"""BH1: the turn-part of a charged horizon. Centre potential M(S,Q) put in; exact sympy."""
import json, os, sympy as sp
S, Q, x = sp.symbols('S Q x', positive=True)
M = sp.sqrt(S/sp.pi)/2*(1 + sp.pi*Q**2/S)           # put in (units G = c = 1, S = area/4)
T, Phi = sp.diff(M, S), sp.diff(M, Q)
DS = lambda g: S*sp.diff(g, S); DQ = lambda g: Q*sp.diff(g, Q)
D = lambda g: 2*DS(g) + DQ(g)
z = lambda e: sp.simplify(e) == 0
X = sp.pi*Q**2/S                                     # pure number: 0 neutral, 1 extremal
out = {}
out['T1_self_similar_2_1'] = z(D(M) - M) and z(D(X))
# charge corner: ratio constant -> commutes exactly
epsF = DS(Phi)/DQ(Phi)
out['T2_charge_corner_commutes'] = z(epsF + sp.Rational(1, 2))
# heat corner: one defect
epsT = DS(T)/DQ(T)
out['T3_eps_heat'] = z(epsT - (1 - 3*X)/(4*X))
w = DS(1/epsT)*T*S/2                                 # DW1-T2
wU = sp.simplify(w/M)
closed = -x*(1 - x)/((1 + x)*(1 - 3*x)**2)
out['T4_turn'] = z(wU - closed.subs(x, X))
out['T4_closed_form'] = str(closed)
# T5: zeros at x=0 (neutral) and x=1 (extremal, T=0); double pole at x=1/3; one sign throughout
out['T5_zeros'] = closed.subs(x, 0) == 0 and closed.subs(x, 1) == 0 and z(T.subs(Q, sp.sqrt(S/sp.pi)))
out['T5_pole'] = sp.limit(closed*(1 - 3*x)**2, x, sp.Rational(1, 3)) == sp.Rational(-1, 6)
out['T5_one_sign'] = all(closed.subs(x, v) < 0 for v in (sp.Rational(1, 10), sp.Rational(3, 10), sp.Rational(2, 5), sp.Rational(9, 10)))
out['T5_pole_is_Q2_over_M2'] = str(sp.simplify((Q**2/M**2).subs(Q, sp.sqrt(S/(3*sp.pi)))))   # 3/4
# T6: BL1-T2 reading: first cycle multiplies exp(I0) by 1+2w^2
g = sp.simplify(1 + 2*closed**2)
out['T6_first_cycle_x_1_2'] = str(g.subs(x, sp.Rational(1, 2)))
out['T6_rest_at_ends'] = g.subs(x, 0) == 1 and g.subs(x, 1) == 1
out['pass'] = all(v for v in out.values() if isinstance(v, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BH1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
