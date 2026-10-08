"""XU1: in an expanding space no frame has the full unit; the count read by the best frame depends on the path.
Exact sympy. Clock factor N^2 = 1 - rs/r - h^2 r^2 (MO1 with CO1's expansion rate h)."""
import json, os, sympy as sp
r, M, h, rb = sp.symbols('r M h r_b', positive=True)
out = {}
N2 = 1 - 2*M/r - h**2*r**2
rstar = (M/h**2)**sp.Rational(1, 3)                        # where the clock factor is largest
out['T1_best_frame'] = sp.simplify(sp.diff(N2, r).subs(r, rstar)) == 0 and sp.simplify(sp.diff(N2, r, 2).subs(r, rstar)) < 0
Nmax2 = sp.simplify(N2.subs(r, rstar))
out['T1_unit'] = sp.simplify(Nmax2 - (1 - 3*(h*M)**sp.Rational(2, 3))) == 0
# T2: flat limit: unit -> 1 (UN2). Largest centre: unit 0 at 27 h^2 M^2 = 1 (no frame can read a count)
out['T2_limits'] = Nmax2.subs(h, 0) == 1 and sp.simplify(Nmax2.subs(M, 1/(sp.sqrt(27)*h))) == 0
# T3: variables (r_b, h): r_b = horizon of the centre; M = r_b (1 - h^2 r_b^2)/2 ; flat count S = pi r_b^2
Mb = rb*(1 - h**2*rb**2)/2
out['T3_horizon'] = sp.simplify(N2.subs({M: Mb, r: rb})) == 0
k = 1/sp.sqrt(Nmax2.subs(M, Mb))                           # UN1-T6 / SM1-T5: count read with unit c_f = N_max
S = sp.pi*rb**2
curl = sp.simplify(sp.diff(k, rb)*sp.diff(S, h) - sp.diff(k, h)*sp.diff(S, rb))   # d(k dS) = dk ^ dS on (r_b, h)
out['T4_defect_form'] = sp.simplify(curl + sp.diff(k, h)*2*sp.pi*rb) == 0
# T5: sign: dk/dh > 0 wherever a best frame exists (h r_b < 1/sqrt 3): the defect has one sign on the whole region
dk = sp.diff(k, h)
pts = [(sp.Rational(1, 2), sp.Rational(1, 5)), (1, sp.Rational(1, 3)), (2, sp.Rational(1, 5)), (3, sp.Rational(1, 10)), (sp.Rational(1, 10), 5)]
out['T5_one_sign'] = all(sp.N(dk.subs({rb: a, h: b})) > 0 for a, b in pts)
x = sp.symbols('x', positive=True)                         # x = h r_b
g = sp.simplify((h*Mb).subs(rb, x/h))                      # h M = x (1 - x^2)/2 : increasing for x < 1/sqrt 3
out['T5_monotone'] = sp.simplify(g - x*(1 - x**2)/2) == 0 and sp.solve(sp.diff(g, x), x) == [sp.sqrt(3)/3]
# T6: one circuit: r_b 1 -> 2 at h = 1/10, h -> 1/5, r_b 2 -> 1, h -> 1/10
import mpmath as mp
fc = sp.lambdify((rb, h), curl, 'mpmath')
val = mp.quad(lambda a_: mp.quad(lambda b_: fc(a_, b_), [mp.mpf(1)/10, mp.mpf(1)/5]), [1, 2])
out['T6_loop_value'] = mp.nstr(val, 8); out['T6_nonzero'] = bool(abs(val) > 0)
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'XU1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
