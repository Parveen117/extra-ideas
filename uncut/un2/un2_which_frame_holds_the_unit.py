"""UN2: which frame holds the unit. GR1's dictionary against UN1-T4. Exact sympy."""
import json, os, sympy as sp
N, n, r, rs, th = sp.symbols('N n r rs theta', positive=True)
out = {}
p = (1 + N)/2                                   # GR1: share
j = (2*p - 1)*n                                 # PR4: current
sig2 = n**2*4*p*(1 - p)                         # PR4: proper density squared
out['T1_legs'] = sp.simplify(j/n - N) == 0 and sp.simplify(sig2/n**2 - (1 - N**2)) == 0
# so the walk's coin at a place is (c, s) = (N, beta): the speed-like leg is N, not the fall speed
Nr = sp.sqrt(1 - rs/r)
out['T2_far_frame_full_unit'] = sp.limit(Nr, r, sp.oo) == 1          # HB1: unit = c = 1 for the far frame
out['T2_horizon_static_unit_zero'] = Nr.subs(r, rs) == 0              # quarter-turn clock: no minimum there
# T3: UN1-T4 used the other leg (fall speed beta -> 1). The two legs are one coin, a quarter turn apart
out['T3_quarter_turn_apart'] = sp.simplify(sp.cos(sp.pi/2 - th) - sp.sin(th)) == 0 and sp.simplify(sp.sin(sp.pi/2 - th) - sp.cos(th)) == 0
# T4: the readings M and T of AL1 are far readings: kappa is the rate against the far clock
kap = sp.simplify(sp.diff(Nr**2, r).subs(r, rs)/2)
out['T4_far_rate'] = sp.simplify(kap - 1/(2*rs)) == 0
A = 4*sp.pi*rs**2; cf = sp.limit(Nr, r, sp.oo)
S = sp.integrate(sp.diff(rs/2, rs)/(cf*kap/(2*sp.pi)), (rs, 0, rs))
out['T4_quarter_with_far_unit'] = sp.simplify(S - A/4) == 0
# T5: a static frame at finite r has unit N < 1: read with its own unit the same count is area/(4 N)
out['T5_local'] = str(sp.simplify((A/(4*N))/(A/4)))                   # 1/N = n, GR1's blueshift
out = {k_: (bool(x) if not isinstance(x, str) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'UN2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
