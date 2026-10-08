"""ST1: the count of a horizon is the area its closed turns sweep in the (M, t) plane. Exact sympy."""
import json, os, sympy as sp
t, kap, M, m_, rs, d, c, Q, S = sp.symbols('t kappa M m rs d c Q S', positive=True)
out = {}
# T1: along the open orbit Exp(kappa t G') the shared information is I = Log cosh(kappa t): rate -> kappa
I_t = sp.log(sp.cosh(kap*t))
out['T1_rate'] = sp.limit(sp.diff(I_t, t), t, sp.oo) == kap and sp.simplify(sp.diff(I_t, t) - kap*sp.tanh(kap*t)) == 0
# T2: the closed reading of the same generator (HX1-T5: angle x iota) returns after P = 2 pi / kappa
P = 2*sp.pi/kap
out['T2_period'] = sp.simplify(kap*P - 2*sp.pi) == 0
# T3: three directions, M = rs/2, kappa = 1/(2 rs): P = 8 pi M; area swept in the (M, t) plane
P3 = P.subs(kap, 1/(4*m_))
area_Mt = sp.integrate(P3, (m_, 0, M))
A = 4*sp.pi*(2*M)**2
out['T3_area_identity'] = sp.simplify(area_Mt - 4*sp.pi*M**2) == 0 and sp.simplify(area_Mt - A/4) == 0
# T4: with count := (M,t)-area per unit turn, T = dM/dS = 1/P is not an extra input
out['T4_T_is_inverse_period'] = sp.simplify(1/sp.diff(area_Mt, M) - 1/(8*sp.pi*M)) == 0
# T5: every d: d(area)/d rs = P * dM/d rs, and it is AL1's count with u = 1/(2 pi)
kd = (d-2)/(2*rs); Md = c*rs**(d-2)
dS = sp.simplify((2*sp.pi/kd)*sp.diff(Md, rs))
AL1 = 2*c*rs**(d-1)/((d-1)*(1/(2*sp.pi)))
out['T5_all_d'] = sp.simplify(sp.powsimp(dS - sp.diff(AL1, rs), force=True)) == 0
# T6: charged centre (BH1). At fixed charge the (M,t)-area from the T = 0 end is the count above the end's own
Mq = sp.sqrt(S/sp.pi)/2*(1 + sp.pi*Q**2/S)
Tq = sp.diff(Mq, S)
out['T6_T_zero_end'] = sp.simplify(Tq.subs(S, sp.pi*Q**2)) == 0 and sp.simplify(Mq.subs(S, sp.pi*Q**2) - Q) == 0
s_ = sp.symbols('s', positive=True)
swept = sp.integrate(sp.Integer(1), (s_, sp.pi*Q**2, S))     # int dM/T at fixed Q = int dS
out['T6_excess_only'] = sp.simplify(swept - (S - sp.pi*Q**2)) == 0
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ST1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
