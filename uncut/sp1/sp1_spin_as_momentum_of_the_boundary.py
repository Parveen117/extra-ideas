"""SP1: the rotating centre from the momentum plane. Spin J carried by the boundary of radius R = sqrt(S/pi) is the
momentum p = J/R of MT1; charge enters as LC1's stored cost. Exact sympy."""
import json, os, sympy as sp
S, J, Q, y = sp.symbols('S J Q y', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
R = sp.sqrt(S/sp.pi)                            # radius of the count (ST1: S = pi R^2)
M0 = R/2                                        # neutral, non-turning centre (AL1/WD1-T4)
E = sp.sqrt(M0**2 + (J/R)**2)                   # MT1: energy read = sqrt(rest^2 + momentum^2)
# T1: this is BH2's centre potential
out['T1_BH2_centre_derived'] = z(E**2 - (S/(4*sp.pi) + sp.pi*J**2/S))
# T2: the partner of spin is the rim speed over the radius: Omega = v / R with v = p / E (MT1-T5)
Om = sp.diff(E, J); v = (J/R)/E
out['T2_partner_is_rim_speed'] = z(Om - v/R)
# T3: BH2's pure number y = 4 pi^2 J^2 / S^2 is (gamma v)^2
yv = 4*sp.pi**2*J**2/S**2
out['T3_y_is_gamma_v_squared'] = z(yv - v**2/(1 - v**2))
# T4: where T = 0 (y = 1): rim speed 1/sqrt 2, and the co-turning reader's unit sqrt(1 - v^2) = 1/sqrt 2 (RT1, UN2)
vy = sp.sqrt(y/(1 + y))
out['T4_extremal_rim_speed'] = z(vy.subs(y, 1) - 1/sp.sqrt(2)) and z(sp.sqrt(1 - vy**2).subs(y, 1) - 1/sp.sqrt(2))
T = sp.diff(E, S)
out['T4_T_zero_there'] = z(T.subs(J, S/(2*sp.pi)))
# T5: with charge: rest value = R/2 + Q^2/(2R) (LC1), then the same momentum rule: ME1's combined potential
M0q = R/2 + Q**2/(2*R)
Eq = sp.sqrt(M0q**2 + (J/R)**2)
ME1 = (S/(4*sp.pi))*(1 + sp.pi*Q**2/S)**2 + sp.pi*J**2/S
out['T5_combined_centre_derived'] = z(Eq**2 - ME1)
# T6: BH2's turn-parts in the rim speed: w_rot/M = -y^2/(2(1+y)) = -(gamma v)^2 v^2 / 2
w_rot = -y**2/(2*(1 + y)); vv = sp.symbols('v', positive=True)
out['T6_w_rot_in_v'] = z(w_rot.subs(y, vv**2/(1 - vv**2)) + (vv**2/(1 - vv**2))*vv**2/2)
out = {k: bool(val) for k, val in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SP1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
