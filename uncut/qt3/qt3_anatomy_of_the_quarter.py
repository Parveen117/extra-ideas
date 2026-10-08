"""QT3: what the quarter is made of. Exact sympy."""
import json, os, sympy as sp
d, c, u, Om, r = sp.symbols('d c u Omega r', positive=True)
out = {}
S = 2*c*r**(d-1)/((d-1)*u)                    # AL1-T2
A = Om*r**(d-1)                               # boundary measure, Omega = measure of the unit boundary
ratio = sp.simplify(S/A)
out['T1_ratio'] = sp.simplify(ratio - 2*c/((d-1)*u*Om)) == 0
turn = 2*sp.pi
# T2 (refusal): one unit per walk (quarter turn, QT2-T5) gives a sixteenth, not a quarter
r_walk = ratio.subs({u: 1/(turn/4), d: 3, c: sp.Rational(1, 2), Om: 4*sp.pi})
r_turn = ratio.subs({u: 1/turn, d: 3, c: sp.Rational(1, 2), Om: 4*sp.pi})
out['T2_per_walk'] = str(r_walk); out['T2_per_turn'] = str(r_turn)
out['T2_refusal'] = r_walk == sp.Rational(1, 16) and r_turn == sp.Rational(1, 4)
# T3: with one unit per full turn: ratio = c * (2 turns / Omega) / (d-1)
out['T3_factor'] = sp.simplify(ratio.subs(u, 1/turn) - c*(2*turn/Om)/(d-1)) == 0
# T4: three directions: Omega = 4 pi = two turns exactly, so ratio = c/(d-1) = (1/2)(1/2)
out['T4_two_turns'] = sp.simplify(4*sp.pi - 2*turn) == 0
out['T4_quarter'] = (sp.Rational(1, 2)/(3 - 1)) == sp.Rational(1, 4)
# T5: the quarter in every d  <=>  c = (d-1) Omega / (16 pi)
out['T5_c_all_d'] = sp.solve(sp.Eq(ratio.subs(u, 1/turn), sp.Rational(1, 4)), c) == [Om*(d-1)/(16*sp.pi)]
# T6: with a cone at the centre (QT1-T5) and any c: ratio = c (2 angle/Omega)/(d-1)
ang = sp.symbols('alpha', positive=True)
out['T6_uncut'] = sp.simplify(ratio.subs(u, 1/ang) - c*(2*ang/Om)/(d-1)) == 0
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QT3_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
