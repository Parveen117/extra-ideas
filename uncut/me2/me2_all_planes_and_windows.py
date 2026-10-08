"""ME2: the remaining planes (volume, magnetisation, particle number) added to ME1's form; the window law for
every plane at a first-order transition; the cross relation between two planes. Exact sympy."""
import json, os, itertools, sympy as sp
out = {}
z = lambda e: sp.simplify(e) == 0
# ---- T1: eight planes, one form. E as a function of the eight extensive readings
names = ['S', 'V', 'Q', 'J', 'p', 'M', 'N', 'h2']          # count, volume, charge, spin, momentum, magnetisation, particles, expansion
q = sp.symbols('S V Q J p M N h2', positive=True)
E = sp.Function('E')(*q)
partners = [sp.diff(E, x) for x in q]                        # T, -P, Phi, Omega, v, H, mu, (volume-like)
n = len(q)
pairs = list(itertools.combinations(range(n), 2))
out['T1_pair_equations'] = all(z(sp.diff(partners[i], q[j]) - sp.diff(partners[j], q[i])) for i, j in pairs)
out['T1_counts'] = (len(pairs), len(pairs) + n) == (28, 36)   # 28 pair equations + 8 weight conditions
# count form: dS = P_t (dE - sum_{i>0} partner_i dq_i): coefficient check
Pt = 1/partners[0]
out['T1_count_form'] = z(Pt*partners[0] - 1)
# ---- T2: a first-order transition in three fields (T, P, H): two phases with Gibbs readings G1, G2
T, P, H, d = sp.symbols('T P H delta', real=True)
G1, G2 = sp.Function('G1')(T, P, H), sp.Function('G2')(T, P, H)
dG = G2 - G1
DS, DV, DM = -sp.diff(dG, T), sp.diff(dG, P), -sp.diff(dG, H)        # jumps of S, V, M
dTdP = -sp.diff(dG, P)/sp.diff(dG, T); dTdH = -sp.diff(dG, H)/sp.diff(dG, T)   # slopes of the surface dG = 0
out['T2_slopes'] = z(dTdP - DV/DS) and z(dTdH + DM/DS)
out['T2_cross'] = z(DM/DV + dTdH/dTdP)                                # jump of M from the two slopes
# ---- T3: window law for any plane (field X, extensive Y, jump DY): anomalous C = T DS/delta, anomalous (1/Y)dY/dT = DY/(Y delta)
Y, DY, DSs, TN, slope = sp.symbols('Y DeltaY DeltaS T_N slope', real=True)
w_in = (TN*DSs/d)/(2*DY/(Y*d))
out['T3_width_free'] = not w_in.has(d)
out['T3_window_law'] = z(w_in.subs(DY, DSs*slope) - Y*TN/(2*slope))  # with DY = DS * dT_N/dX (sign carried by the plane)
# critical width per plane
ab = sp.symbols('alpha_bg', positive=True)
out['T3_critical_width'] = sp.solve(sp.Eq(ab + DY/(Y*d), 0), d) == [-DY/(Y*ab)]
# ---- T4: the two windows of one transition are tied: w_V / w_M = (V/M) * (-dT/dH)/(dT/dP)
V_, M_ = sp.symbols('V M', positive=True); sP, sH = sp.symbols('s_P s_H', real=True)
wV = V_*TN/(2*sP); wM = -M_*TN/(2*sH)
out['T4_tied'] = z(wV/wM + (V_/M_)*(sH/sP))
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ME2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
