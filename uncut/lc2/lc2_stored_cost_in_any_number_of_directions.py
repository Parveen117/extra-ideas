"""LC2: LC1's stored-cost rule in d = 3..6 directions, and the heat reading of a charged centre. Exact sympy."""
import json, os, sympy as sp
from sympy.calculus.euler import euler_equations
out = {}
z = lambda e: sp.simplify(e) == 0
r, x, q, c, M, rp = sp.symbols('r x q c M r_plus', positive=True)
ok1 = ok2 = ok3 = True
for d in (3, 4, 5, 6):
    Phi = sp.Function('Phi')
    eq = euler_equations(sp.Rational(1, 2)*r**(d-1)*sp.diff(Phi(r), r)**2, Phi(r), r)[0]     # MC1's weight in d directions
    sol = q/r**(d-2)
    ok1 = ok1 and z(eq.lhs.subs(Phi(r), sol).doit())
    stored = sp.integrate(sp.Rational(1, 2)*x**(d-1)*sp.diff(q/x**(d-2), x)**2, (x, r, sp.oo))
    ok2 = ok2 and z(stored - (d-2)*q**2/(2*r**(d-2)))
    Mr = M - stored                                           # centre value read at r
    N2 = 1 - Mr/(c*r**(d-2))                                  # memory normalised so that the neutral horizon is M = c R^(d-2)
    Mh = sp.solve(sp.Eq(N2.subs(r, rp), 0), M)[0]
    F2 = (d-2)*q**2/(2*c*rp**(2*(d-2)))                       # pure number Phi_hat^2
    s = 4*sp.pi*c/(d-1)                                       # count S = s r^(d-1), fixed by the neutral T0 = (d-2)/(4 pi r)
    T = sp.diff(Mh, rp)/sp.diff(s*rp**(d-1), rp)
    ok3 = ok3 and z(Mh - c*rp**(d-2)*(1 + F2)) and z(4*sp.pi*rp*T - (d-2)*(1 - F2))
out['T1_least_cost_partner'] = ok1
out['T2_stored_cost'] = ok2
out['T3_heat_reading_charged'] = ok3
# T4: three directions, charge and turning together at the same horizon radius: additive (from SP1/LC1's centre)
S, J, Q = sp.symbols('S J Q', positive=True)
R = sp.sqrt(S/sp.pi)
E = sp.sqrt((R/2 + Q**2/(2*R))**2 + J**2/R**2)
T3 = sp.diff(E, S); v2 = (J/R)**2/E**2; rplus = R*sp.sqrt(1 - v2)
pts = [(sp.Rational(3), sp.Rational(1, 3), sp.Rational(1, 2)), (sp.Rational(5), sp.Rational(1), sp.Rational(1, 4)), (sp.Rational(2), sp.Rational(1, 10), sp.Rational(1, 3))]
out['T4_additive_d3'] = all(abs(sp.N((4*sp.pi*rplus*T3 - (1 - 2*v2 - Q**2/R**2)).subs({S: p[0], J: p[1], Q: p[2]}))) < 1e-12 for p in pts)
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LC2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
