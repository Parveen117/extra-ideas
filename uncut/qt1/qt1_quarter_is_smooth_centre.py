"""QT1: the quarter of the horizon count is the smooth-centre condition F4. Exact sympy."""
import json, os, sympy as sp
r, rs, d, rho, kap, t, th, Pd = sp.symbols('r rs d rho kappa t theta P', positive=True)
R = sp.Matrix([[0, -1], [1, 0]]); K = sp.Matrix([[1, 0], [0, -1]]); Sx = R*K; I2 = sp.eye(2)
Exp = lambda g, a: (g*a).exp().applyfunc(sp.simplify)
out = {}
# T1: near the horizon the clock rate is N = kappa * (proper distance), every d
N2 = 1 - (rs/r)**(d-2); k0 = (d-2)/(2*rs)
lead = sp.diff(N2, r).subs(r, rs)                       # N^2 ~ 2 kappa (r - rs); rho = int dr/N -> N = kappa rho
out['T1_rate'] = sp.simplify(lead - 2*k0) == 0 and sp.simplify(sp.sqrt(2*kap*(kap*rho**2/2)) - kap*rho) == 0
# T2: the frame's orbit at fixed rho is a boost: open, never returns
B = Exp(Sx, kap*t)
out['T2_boost_open'] = sp.simplify(B - sp.Matrix([[sp.cosh(kap*t), sp.sinh(kap*t)], [sp.sinh(kap*t), sp.cosh(kap*t)]])) == sp.zeros(2) \
    and Sx**2 == I2 and sp.solve(sp.Eq(sp.cosh(th), 1), th) == []
# T3: the same orbit read by the turn (iota^2 = -1) closes, first at angle 2 pi; half-way it is -1
Tn = Exp(R, th)
out['T3_turn_closes'] = R**2 == -I2 and Tn.subs(th, 2*sp.pi) == I2 and Tn.subs(th, sp.pi) == -I2 \
    and all(Tn.subs(th, sp.pi*q) != I2 for q in (sp.Rational(1, 2), 1, sp.Rational(3, 2)))
# T4: a cycle of period P in t sweeps angle kappa P; circle of radius rho has length kappa P rho.
# length / radius = 2 pi (no cone at the centre: F4)  <=>  P = 2 pi / kappa
circ = sp.integrate(kap*rho, (t, 0, Pd))
out['T4_smooth_iff'] = sp.solve(sp.Eq(circ/rho, 2*sp.pi), Pd) == [2*sp.pi/kap]
# T5: one unit per closed cycle: T = 1/P = u kappa  ->  u = 1/(kappa P); AL1: S = area/(8 pi u)
A = sp.symbols('A', positive=True)
u = 1/(kap*Pd); Scount = A/(8*sp.pi*u)
out['T5_quarter'] = sp.simplify(Scount.subs(Pd, 2*sp.pi/kap) - A/4) == 0
out['T5_general'] = sp.simplify(Scount - (A/4)*(kap*Pd/(2*sp.pi))) == 0   # with a cone: quarter x (angle / 2 pi)
# T6: the three kinds of generator u K + v S + w iota: square = u^2+v^2-w^2 ; closes only when negative
uu, vv, ww = sp.symbols('u_ v_ w_', real=True)
g = uu*K + vv*Sx + ww*R
out['T6_square'] = sp.simplify(g**2 - (uu**2 + vv**2 - ww**2)*I2) == sp.zeros(2)
out = {k_: (bool(v) if not isinstance(v, str) else v) for k_, v in out.items()}
out['pass'] = all(v for v in out.values() if isinstance(v, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QT1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
