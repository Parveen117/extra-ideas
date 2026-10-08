"""FD1: the frame defect of a substance as a measurable pure number. Exact sympy."""
import json, os, sympy as sp
T, V, S, a, b, Nk, cv, B0 = sp.symbols('T V S a b Nk c_v B0', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
R_ = sp.Rational
def responses(P, CV):
    """from P(T,V) and C_V(T,V): alpha, C_P, and Hh = C_P/(alpha T)"""
    PT, PV = sp.diff(P, T), sp.diff(P, V)
    KT = -V*PV
    alpha = PT/KT
    CP = CV + T*V*alpha**2*KT
    return alpha, CP, sp.simplify(CP/(alpha*T))
D_of = lambda Hh, CV: sp.simplify((T/CV)*sp.diff(Hh, T))          # D = (d Hh / d S) at fixed V
# T1: DW1's second turn-part in measurable terms. With P(S,V): V P_V / P_S = -C_P/(alpha T)  [entropy units]
Pid = Nk*T/V; CVid = cv*Nk
al, CP, Hh = responses(Pid, CVid)
out['T1_ideal'] = z(Hh - CP) and z(D_of(Hh, CVid))                 # Hh = C_P, constant: D = 0
# so w2 = (1/2) S d_S(eps^-1) T S = (1/2) T Hh = C_P T/2 = H/2 ;  w2/U = gamma/2 = (1/2) exp(2 I_cut) (CI1)
gam = (cv + 1)/cv
out['T1_turn_is_half_enthalpy'] = z(R_(1, 2)*T*Hh - (cv + 1)*Nk*T/2) and z((R_(1, 2)*T*Hh)/(cv*Nk*T) - gam/2)
# T2: pure powers (no count unit): radiation U ~ S^(4/3) V^(-1/3), horizon M ~ S^(1/2): eps constant, both defects zero
U_rad = S**R_(4, 3)*V**R_(-1, 3); P_rad = -sp.diff(U_rad, V)
eps_rad = sp.simplify(S*sp.diff(P_rad, S)/(V*sp.diff(P_rad, V)))
out['T2_radiation_commutes'] = eps_rad == -1 and z(sp.diff(eps_rad, S))
# ideal gas in (S,V): P ~ V^-gamma exp(S/C_V): eps = -S/(gamma C_V): runs with S. The scale is the count unit C_V ~ Nk
g_, C_ = sp.symbols('gamma C', positive=True)
P_ig = V**(-g_)*sp.exp(S/C_)
eps_ig = sp.simplify(S*sp.diff(P_ig, S)/(V*sp.diff(P_ig, V)))
out['T2_ideal_gas_scale_is_count_unit'] = z(eps_ig + S/(g_*C_)) and z(S*sp.diff(eps_ig, S) + C_*sp.diff(eps_ig, C_))   # SS1 with b0 = C
# T3: hard cores alone give no frame defect
al, CP, Hh = responses(Nk*T/(V - b), cv*Nk)
out['T3_hard_core_no_defect'] = z(D_of(Hh, cv*Nk))
# T4: attraction is the source. van der Waals: D = 2 a (V - b) / (Nk T V^2)
al, CP, Hh = responses(Nk*T/(V - b) - a/V**2, cv*Nk)
Dv = D_of(Hh, cv*Nk)
out['T4_vdW'] = z(Dv - 2*a*(V - b)/(Nk*T*V**2))
Tr, vr = sp.symbols('T_r v_r', positive=True)
red = Dv.subs({V: 3*b*vr, T: Tr*8*a/(27*Nk*b)})
out['T4_reduced'] = z(red - R_(9, 4)*(1 - 1/(3*vr))/(Tr*vr))       # one function for every such substance (SS1)
# T5: low density, any substance with second virial B(T), C_V = c_v Nk + correction: leading D in 1/V
Bf = sp.Function('B')(T)
Pv = Nk*T/V*(1 + Bf/V)
CVv = cv*Nk - (Nk/V)*(2*T*sp.diff(Bf, T) + T**2*sp.diff(Bf, T, 2))  # from d C_V/d V = T P_TT
out['T5_CV_consistent'] = z(sp.diff(CVv, V) - T*sp.diff(Pv, T, 2))
al, CP, Hh = responses(Pv, CVv)
eps_ = sp.symbols('epsilon', positive=True)
Dl = sp.series(D_of(Hh, CVv).subs(V, 1/eps_), eps_, 0, 2).removeO()
lead = sp.simplify(Dl/eps_)
out['T5_leading'] = str(sp.simplify(lead))
out['T5_matches_vdW'] = z(lead.subs(Bf, b - a/(Nk*T)).doit() - 2*a/(Nk*T))
# T6: D' = (T/C_V) d[C_P/(alpha T) - C_V]/dT: zero for an ideal gas with ANY C_V(T); same as D for van der Waals;
#     low density: D' V -> [2 T B' - (c_v - 1) T^2 B''] / c_v
cT = sp.Function('c')(T)
al6, CP6, Hh6 = responses(Nk*T/V, cT*Nk)
out['T6_ideal_any_heat_capacity'] = z((T/(cT*Nk))*sp.diff(Hh6 - cT*Nk, T))
Dp = sp.simplify((T/CVv)*sp.diff(Hh - CVv, T))
lead6 = sp.simplify(sp.series(Dp.subs(V, 1/eps_), eps_, 0, 2).removeO()/eps_)
out['T6_low_density'] = z(lead6 - (2*T*sp.diff(Bf, T) - (cv - 1)*T**2*sp.diff(Bf, T, 2))/cv)
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FD1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
