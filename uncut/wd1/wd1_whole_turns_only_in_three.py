"""WD1: the boundary of d directions is a whole number of turns only for d = 2 and d = 3. Exact sympy."""
import json, os, sympy as sp
R_ = sp.Rational
out = {}
turn = 2*sp.pi
Om = lambda d: 2*sp.pi**R_(d, 2)/sp.gamma(R_(d, 2))          # measure of the unit boundary in d directions
ratio = {d: sp.simplify(Om(d)/turn) for d in range(2, 41)}
out['T1_d2_d3'] = ratio[2] == 1 and ratio[3] == 2
out['T1_only'] = [d for d in ratio if ratio[d].is_rational] == [2, 3]
# closed forms: even d = 2m: pi^(m-1)/(m-1)! ; odd d = 2m+1: 2^m pi^(m-1)/(2m-1)!!  -> power of pi is m-1
m = sp.symbols('m', integer=True, positive=True)
out['T2_even_form'] = all(sp.simplify(ratio[2*k] - sp.pi**(k-1)/sp.factorial(k-1)) == 0 for k in range(1, 20))
out['T2_odd_form'] = all(sp.simplify(ratio[2*k+1] - 2**k*sp.pi**(k-1)/sp.factorial2(2*k-1)) == 0 for k in range(1, 20))
# T3: d = 2 has no scale (MC1: log r; AL1 degree 0): the centre value does not depend on the scale
d, rs, c = sp.symbols('d rs c', positive=True)
out['T3_d2_no_scale'] = sp.diff((c*rs**(d-2)).subs(d, 2), rs) == 0 and ((d-2)/(d-1)).subs(d, 2) == 0
# T4: d = 3: c from the line's own period law (MO1: (dphi/dt)^2 = rs/(2 r^3) = M/r^3) -> M = rs/2
r = sp.symbols('r', positive=True)
N2 = 1 - rs/r
w2 = sp.simplify(sp.diff(N2, r)/(2*r))
out['T4_c_half'] = sp.simplify(w2 - (rs/2)/r**3) == 0
# T5: horizon count / boundary with one unit per full turn (QT3-T3): rational only at d = 3 (c rational)
cnt = lambda dd: sp.simplify(R_(1, 2)*(2*turn/Om(dd))/(dd-1))
out['T5_count_share_d3'] = cnt(3) == R_(1, 4)
out['T5_rational_only_d3'] = [dd for dd in range(3, 41) if cnt(dd).is_rational] == [3]
out = {k_: (bool(x) if not isinstance(x, (str, list)) else x) for k_, x in out.items()}
out['pass'] = all(x for x in out.values() if isinstance(x, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'WD1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
