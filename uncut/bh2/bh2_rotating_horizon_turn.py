"""BH2: the turn-parts of a rotating horizon. Centre potential M(S,J) put in; exact sympy."""
import json, os, sympy as sp
S, J, y = sp.symbols('S J y', positive=True)
M = sp.sqrt(S/(4*sp.pi) + sp.pi*J**2/S)             # put in (G = c = 1, S = area/4)
T, Om = sp.diff(M, S), sp.diff(M, J)
DS = lambda g: S*sp.diff(g, S); DJ = lambda g: J*sp.diff(g, J)
sub = lambda e: sp.simplify(e.subs(J, sp.sqrt(y)*S/(2*sp.pi)))   # y = 4 pi^2 J^2 / S^2 ; 1 = extremal
z = lambda e: sp.simplify(e) == 0
R = sp.Rational
out = {}
out['T1_self_similar_2_2'] = z(2*DS(M) + 2*DJ(M) - M) and z(sub(T).subs(y, 1))
epsO, epsT = DS(Om)/DJ(Om), DS(T)/DJ(T)
out['T2_eps'] = z(sub(epsO) + (y + 3)/2) and z(sub(epsT) + (3*y**2 + 6*y - 1)/(2*y*(y + 3)))
w_rot = -y**2/(2*(y + 1))
w_heat = y*(y - 1)*(3*y**2 + 2*y + 3)/((y + 1)*(3*y**2 + 6*y - 1)**2)
out['T3_w_rot'] = z(sub(DJ(epsO)*Om*J/(2*M)) - w_rot)
out['T4_w_heat'] = z(sub(DS(1/epsT)*T*S/(2*M)) - w_heat)
y0 = 2/sp.sqrt(3) - 1
out['T5_heat_shape'] = (w_heat.subs(y, 0) == 0 and w_heat.subs(y, 1) == 0 and z((3*y**2 + 6*y - 1).subs(y, y0))
                        and all(w_heat.subs(y, v) < 0 for v in (R(1, 10), R(1, 7), R(1, 5), R(1, 2), R(9, 10))))
out['T5_pole_J2_over_M4'] = str(sp.nsimplify(sp.simplify(4*y0/(1 + y0)**2)))      # 2 sqrt(3) - 3
out['T6_rot_at_extremal'] = w_rot.subs(y, 1) == R(-1, 4) and sp.diff(w_rot, y).subs(y, R(1, 2)) < 0
out['T6_first_cycle_extremal'] = str(1 + 2*w_rot.subs(y, 1)**2)                    # 9/8
out = {k: (bool(v) if not isinstance(v, str) else v) for k, v in out.items()}
out['pass'] = all(v for v in out.values() if isinstance(v, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BH2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
