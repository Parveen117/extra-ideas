"""AL1: the count of a one-scale horizon. Exact sympy, any dimension d."""
import json, os, sympy as sp
r, c, u, d, S = sp.symbols('r c u d S', positive=True)
out = {}
# put in: MC1 potential r^-(d-2) -> fall speed beta^2 = (rs/r)^(d-2) (MO1 at d=3); GR1 N^2 = 1 - beta^2
rs = sp.symbols('rs', positive=True)
N2 = 1 - (rs/r)**(d-2)
kappa = sp.simplify(sp.diff(N2, r).subs(r, rs)/2)            # clock-rate gradient at the horizon
out['T1_rate'] = sp.simplify(kappa - (d-2)/(2*rs)) == 0
M = c*rs**(d-2)                                               # one scale: centre value fixed by rs
T = u*kappa                                                   # put in: reading = unit x horizon rate
Sc = 2*c*rs**(d-1)/((d-1)*u)                                  # ~ rs^(d-1): the boundary measure
out['T2_count'] = sp.simplify(sp.powsimp(sp.diff(Sc, rs)*T/sp.diff(M, rs), force=True)) == 1   # dS = dM/T, zero with no horizon
# T3: degree of the centre potential in the count, and its heat ratio
k = sp.simplify(sp.diff(sp.log(M), rs)/sp.diff(sp.expand_log(sp.log(Sc), force=True), rs))
out['T3_degree'] = sp.simplify(k - (d-2)/(d-1)) == 0
out['T3_heat_ratio_is_ID1_ratio'] = sp.simplify((k - 1) + 1/(d-1)) == 0   # S T_S / T = -1/(d-1)
# T4: d = 3, c = 1/2 : count = area / (8 pi u); quarter <=> u = 1/(2 pi)
S3 = Sc.subs({d: 3, c: sp.Rational(1, 2)}); A = 4*sp.pi*rs**2
out['T4_d3'] = sp.simplify(S3 - A/(8*sp.pi*u)) == 0
out['T4_quarter_iff'] = sp.solve(sp.Eq(S3, A/4), u) == [1/(2*sp.pi)]
# T5: with that unit BH1's centre is recovered: M = sqrt(S/pi)/2
out['T5_BH1_centre'] = sp.simplify((rs/2) - sp.sqrt(S3.subs(u, 1/(2*sp.pi))/sp.pi)/2) == 0
# T6: for any u the turn laws of BH1/BH2 are unchanged in shape: w/M depends on S only through pure numbers
out['T6_unit_free_degree'] = not k.has(u) and not k.has(c)
out = {k_: (bool(v) if not isinstance(v, str) else v) for k_, v in out.items()}
out['pass'] = all(v for v in out.values() if isinstance(v, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'AL1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
