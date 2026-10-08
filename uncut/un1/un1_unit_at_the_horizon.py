"""UN1: the unit of the horizon frame and the count. Exact sympy."""
import json, os, sympy as sp
R_ = sp.Rational
A, cf, Q, beta = sp.symbols('A c_f Q beta', positive=True)
out = {}
unit = lambda q: sp.cos(2*sp.pi/q)                       # HB1 ladder: even part of a tick on q marks
out['T1_ladder'] = unit(2) == -1 and unit(4) == 0 and sp.simplify(unit(8) - 1/sp.sqrt(2)) == 0 and sp.limit(unit(Q), Q, sp.oo) == 1
# T2 (refusal): one tick is the turn 2 pi/Q; it is a full turn only for Q = 1 (no clock at all)
out['T2_tick_not_full_turn'] = sp.solve(sp.Eq(2*sp.pi/Q, 2*sp.pi), Q) == [1]
# T3: the minimum of count x response is (1/2) c_f (HB1-T2/T3). One such unit per closed turn:
u = cf/(2*sp.pi)
S = A/(8*sp.pi*u)                                        # AL1-T4
out['T3_count'] = sp.simplify(S - A/(4*cf)) == 0
# T4: HB1 section 3: unit = cone speed = 1 - curvature/2. MO1: the frame's speed at the horizon is 1.
N2 = 1 - beta**2
out['T4_horizon_is_lightlike'] = sp.solve(sp.Eq(N2, 0), beta) == [1]
out['T4_quarter'] = sp.simplify(S.subs(cf, 1) - A/4) == 0
# T5: what other clocks would give
out['T5_eight_mark'] = str(sp.simplify(S.subs(cf, unit(8))/A))        # sqrt(2)/4
out['T5_quarter_turn_clock_no_count'] = sp.limit(1/S, cf, 0) == 0     # unit 0: count unbounded, no minimum
# T6: uncut form with QT1-T5 cone angle alpha and unit c_f together
al = sp.symbols('alpha', positive=True)
out['T6_general'] = sp.simplify(A/(8*sp.pi*(cf/al)) - (A/4)*(al/(2*sp.pi))/cf) == 0
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'UN1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
