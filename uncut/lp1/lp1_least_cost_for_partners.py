"""LP1: MC1's least-cost rule applied to the partners. Commuting plane: 1/r, and its stored cost is a source in
the clock plane - which gives BH1's centre potential instead of putting it in. Non-commuting plane: the
stationarity condition carries the quadratic term. Exact sympy."""
import json, os, sympy as sp
from sympy.calculus.euler import euler_equations
r, Q, M, S, rp = sp.symbols('r Q M S r_plus', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
# ---- T1: one commuting partner Phi(r), three directions, cost (1/2) int r^2 Phi'^2 dr  (MC1 with w = r^(d-2), continuous form)
Phi = sp.Function('Phi')
eq = euler_equations(sp.Rational(1, 2)*r**2*sp.diff(Phi(r), r)**2, Phi(r), r)[0]
sol = Q/r
out['T1_least_cost_is_1_over_r'] = z(eq.lhs.subs(Phi(r), sol).doit()) and z(r**2*sp.diff(sol, r) + Q)      # constant flux Q
# ---- T2: cost stored outside the shell r
x = sp.symbols('x', positive=True)
stored = sp.integrate(sp.Rational(1, 2)*x**2*sp.diff(Q/x, x)**2, (x, r, sp.oo))
out['T2_stored_outside'] = z(stored - Q**2/(2*r))
# ---- T3: the centre value read at r is the far value less what is stored outside; memory m = 2 M(r)/r (GR1/MC1)
Mr = M - stored
N2 = sp.expand(1 - 2*Mr/r)
out['T3_clock_factor'] = z(N2 - (1 - 2*M/r + Q**2/r**2))
# ---- T4: horizon N^2 = 0 at r_plus; count = pi r_plus^2 (ST1/AL1): the centre potential of BH1 follows
Mh = sp.solve(sp.Eq(N2.subs(r, rp), 0), M)[0]
out['T4_BH1_centre_derived'] = z(Mh.subs(rp, sp.sqrt(S/sp.pi)) - sp.sqrt(S/sp.pi)/2*(1 + sp.pi*Q**2/S))
# and the partner at the horizon is the Phi of BH1
out['T4_partner_at_horizon'] = z(sp.diff(sp.sqrt(S/sp.pi)/2*(1 + sp.pi*Q**2/S), Q) - (Q/rp).subs(rp, sp.sqrt(S/sp.pi)))
# ---- T5: non-commuting plane, two base variables: cost density (1/2) sum_c (f^c)^2, f^c = d a^c - s eps a^a a^b (NA1-T2)
X, Y = sp.symbols('X Y', real=True); s = sp.symbols('s', real=True)
a = [[sp.Function(f'a{c}{i}')(X, Y) for i in range(2)] for c in range(3)]
f = [sp.diff(a[c][1], X) - sp.diff(a[c][0], Y) - s*sum(sp.LeviCivita(c, p, q)*a[p][0]*a[q][1] for p in range(3) for q in range(3)) for c in range(3)]
Lag = sp.Rational(1, 2)*sum(fc**2 for fc in f)
eqs = euler_equations(Lag, [a[c][i] for c in range(3) for i in range(2)], [X, Y])
# expected: for partner c, component x:  d_Y f^c - s eps_{cpq} a^p_Y f^q = 0  (the other partners enter as sources)
exp_x = [sp.diff(f[c], Y) - s*sum(sp.LeviCivita(c, p, q)*a[p][1]*f[q] for p in range(3) for q in range(3)) for c in range(3)]
got_x = [eqs[2*c].lhs for c in range(3)]
out['T5_stationarity_has_quadratic_term'] = all(sp.simplify(sp.expand(got_x[c] - exp_x[c])) == 0 or sp.simplify(sp.expand(got_x[c] + exp_x[c])) == 0 for c in range(3))
# one cut active: reduces to d f = 0, MC1's constant flux
one = {a[1][0]: 0, a[1][1]: 0, a[2][0]: 0, a[2][1]: 0}
out['T5_one_cut_is_constant_flux'] = sp.simplify(exp_x[0].subs(one).doit() - sp.diff(sp.diff(a[0][1], X) - sp.diff(a[0][0], Y), Y)) == 0
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LP1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
