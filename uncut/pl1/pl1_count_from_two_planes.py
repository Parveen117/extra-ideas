"""PL1: the count as area in two planes - (M, t) and (charge or spin, its angle). Exact sympy."""
import json, os, sympy as sp
M, Q, J, q, j_ = sp.symbols('M Q J q j', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
# ---- charged centre, in the variables (M, Q)
rp = M + sp.sqrt(M**2 - Q**2)
Sc = sp.pi*rp**2
Pc = sp.simplify(sp.diff(Sc, M))                    # period of the closed turn = dS/dM
angQ = sp.simplify(-sp.diff(Sc, Q))                 # angle swept in the charge direction over one period
Phi = Q/rp
out['T1_gauge_angle_is_Phi_P'] = z(angQ - Phi*Pc)
out['T1_closed_form'] = z(sp.diff(Pc, Q) + sp.diff(angQ, M))        # one equation of the Maxwell type
# path: Q = 0 up to M, then across at fixed M
a1 = sp.integrate(Pc.subs(Q, 0), (M, 0, M))
a2 = sp.integrate(sp.simplify(angQ.subs(Q, q)), (q, 0, Q))
out['T2_two_areas'] = z(a1 - 4*sp.pi*M**2) and z(a1 - a2 - Sc)
ext = sp.limit(a2, Q, M)
out['T2_extremal_charge_plane'] = z(ext - 3*sp.pi*M**2) and z(a1 - ext - sp.pi*M**2)
# ---- rotating centre, in the variables (M, J)
rk = M + sp.sqrt(M**2 - J**2/M**2)
Sk = 2*sp.pi*M*rk
Pk = sp.simplify(sp.diff(Sk, M)); angJ = sp.simplify(-sp.diff(Sk, J))
Om = (J/M)/(2*M*rk)
out['T3_rotation_angle_is_Omega_P'] = z(angJ - Om*Pk)
out['T3_closed_form'] = z(sp.diff(Pk, J) + sp.diff(angJ, M))
b2 = sp.integrate(sp.simplify(angJ.subs(J, j_)), (j_, 0, J))
out['T4_two_areas'] = z(4*sp.pi*M**2 - b2 - Sk)
extk = sp.limit(b2, J, M**2)
out['T4_extremal_spin_plane'] = z(extk - 2*sp.pi*M**2) and z(4*sp.pi*M**2 - extk - 2*sp.pi*M**2)
# T5: at the extremal end the angle swept per period: charge -> unbounded with P; ratio angle/P = Phi -> 1, Omega -> 1/(2M)
out['T5_rates_at_end'] = sp.limit(Phi, Q, M) == 1 and z(sp.limit(Om, J, M**2) - 1/(2*M))
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PL1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
