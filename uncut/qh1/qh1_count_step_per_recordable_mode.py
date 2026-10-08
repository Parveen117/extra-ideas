"""QH1: is the count of a horizon a whole number? The line's own WQ1 and BC1 against today's ST1/CK1.
Exact sympy."""
import json, os, sympy as sp
out = {}
z = lambda e: sp.simplify(e) == 0
ku, w, kh, n, N, wl, M = sp.symbols('kappa_u omega kappa_h n N omega_loc M', positive=True)
# WQ1-W2: a recordable mode of rate omega has energy n * kappa_u * omega (kappa_u: the unit; its value is free, SC1/DC1)
dM = ku*w                                         # one count of the mode given to the centre
# ST1/SM1 with the unit written out: T = kappa_u * kappa_h / (2 pi)  (unit x horizon rate per turn); dS = dM / T
T = ku*kh/(2*sp.pi)
dS = dM/T
out['T1_step'] = z(dS - 2*sp.pi*w/kh)
out['T1_unit_drops_out'] = z(sp.diff(dS, ku))
# T2: in turns, the step is the number of turns the mode makes during one closing turn of the horizon (CK1)
P = 2*sp.pi/kh
out['T2_turns_of_the_mode_per_closed_turn'] = z(dS/(2*sp.pi) - w*P/(2*sp.pi))
# T3: whole exactly when omega/kappa_h is whole; a mode at the horizon's own rate adds one turn: dS = 2 pi
out['T3_own_rate'] = z(dS.subs(w, kh) - 2*sp.pi)
# three directions, S = area/4: the boundary grows by 8 pi (in units of the unit area) per count of such a mode
out['T3_area_step'] = z(4*dS.subs(w, kh) - 8*sp.pi)
# T4: BC1-B3: a mode of local rate omega_loc at clock factor N is read far at omega_loc * N: the step closes up at the boundary
dS_loc = dS.subs(w, wl*N)
out['T4_steps_close_up'] = sp.limit(dS_loc, N, 0) == 0 and z(dS_loc - 2*sp.pi*wl*N/kh)
# T5: neutral centre, three directions: kappa_h = 1/(4M): step = 8 pi M omega - any real number; nothing makes it whole
out['T5_not_forced_whole'] = z(dS.subs(kh, 1/(4*M)) - 8*sp.pi*M*w) and (8*sp.pi*M*w).subs({M: sp.Rational(1, 3), w: sp.Rational(1, 7)})/(2*sp.pi) == sp.Rational(4, 21)
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QH1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
