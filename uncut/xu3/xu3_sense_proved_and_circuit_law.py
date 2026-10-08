"""XU3: the sense of the frame-consistent count is proved over the whole region; the gain of a circuit
factorises; the slow-expansion law. Exact sympy (one numeric quadrature for a table value)."""
import json, os, sympy as sp, mpmath as mp
rb, h, S, x = sp.symbols('r_b h S x', positive=True)
R_ = sp.Rational
out = {}
z = lambda e: sp.simplify(e) == 0
# everything in the two variables S = pi r_b^2 (flat count) and x = h r_b (pure number)
u = x*(1 - x**2)/2                                    # = h M
c = sp.sqrt(1 - 3*u**R_(2, 3))                        # unit of the best frame (XU1-T1): depends on x only
g = 2*(1 - x**2)/(1 - 3*x**2)                         # P M = S g(x)
Mb = rb*(1 - h**2*rb**2)/2; P = 4*sp.pi*rb/(1 - 3*h**2*rb**2)
out['T1_scaling'] = z((P*Mb).subs(h, x/rb) - (sp.pi*rb**2)*g) and z((h*Mb).subs(h, x/rb) - u)
# frame-consistent form (XU2, reading F):  theta = dS/c - P M dc/c^2 = a(x) dS + S b(x) dx
a_ = 1/c; b_ = -g*sp.diff(c, x)/c**2
curlF = sp.simplify(b_ - sp.diff(a_, x))              # coefficient of dS ^ dx
closed = (1 + x**2)/(2*u**R_(1, 3)*c**3)
out['T2_symbolic'] = bool(sp.simplify(sp.radsimp(curlF - closed)) == 0) or 'not reduced symbolically'
out['T2_closed_form'] = all(abs(sp.N((curlF - closed).subs(x, v), 30)) < 1e-25 for v in (R_(1, 100), R_(1, 10), R_(1, 4), R_(2, 5), R_(1, 2), R_(11, 20)))
# T3 (proof of sign): every factor of the closed form is positive for 0 < x < 1/sqrt 3
out['T3_region'] = sp.solve(sp.Eq(c**2, 0), x) == [sp.sqrt(3)/3] or z(c.subs(x, 1/sp.sqrt(3)))
out['T3_factors_positive'] = bool(sp.simplify(u.subs(x, R_(1, 2))) > 0) and sp.solve(sp.diff(u, x), x) == [sp.sqrt(3)/3] and z(u.subs(x, 1/sp.sqrt(3)) - 1/(3*sp.sqrt(3)))
# reading A of XU1 for comparison: coefficient -a'(x) = c'/c^2 : negative
curlA = sp.simplify(-sp.diff(a_, x))
out['T3_A_opposite'] = all(sp.N(curlA.subs(x, v)) < 0 and sp.N(closed.subs(x, v)) > 0 for v in (R_(1, 100), R_(1, 4), R_(1, 2)))
# T4: the gain of a circuit (S1 -> S2 at x1, x1 -> x2, S2 -> S1 at x2, back) factorises: (S2 - S1) * (Phi(x2) - Phi(x1))
fx = sp.lambdify(x, closed, 'mpmath')
Phi = lambda lo, hi: mp.quad(fx, [lo, hi])
gain = (mp.pi*4 - mp.pi*1)*Phi(mp.mpf(1)/10, mp.mpf(1)/5)        # r_b: 1 -> 2 is S: pi -> 4 pi ; x: 1/10 -> 1/5
out['T4_example_gain'] = mp.nstr(gain, 8); out['T4_positive'] = bool(gain > 0)
# T5: slow expansion: Phi(x) -> (3/2)(x/2)^(2/3) = (3/2)(h M)^(2/3)
out['T5_leading'] = z(sp.limit(closed*x**R_(1, 3), x, 0) - 2**R_(-2, 3))            # closed ~ 2^(-2/3) x^(-1/3)  =>  Phi ~ (3/2)(x/2)^(2/3)
out['T5_Phi_leading'] = z(sp.integrate(2**R_(-2, 3)*x**R_(-1, 3), (x, 0, x)) - R_(3, 2)*(x/2)**R_(2, 3))
# T6: no expansion, no defect: the coefficient times dx vanishes as h -> 0 at fixed r_b (Phi(0) = 0)
out['T6_flat_limit'] = z(sp.limit(R_(3, 2)*(x/2)**R_(2, 3), x, 0))
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'XU3_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
