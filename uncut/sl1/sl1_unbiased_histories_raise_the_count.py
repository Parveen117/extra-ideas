"""SL1: unbiased (symmetric) histories raise the count. GE2-T3 steps +-v applied to the horizon count. Exact sympy + sampled signs."""
import json, os, sympy as sp, mpmath as mp
eta, r, M, d_, k = sp.symbols('eta r M delta k', positive=True); R_ = sp.Rational
out = {}
z = lambda e: sp.simplify(e.rewrite(sp.exp)) == 0
# T1: one symmetric step +-eta in Log r (weights 1/2, 1/2) on the flat count S = pi r^2
S = lambda rr: sp.pi*rr**2
mean = (S(r*sp.exp(eta)) + S(r*sp.exp(-eta)))/2
out['T1_mean_count'] = z(mean/S(r) - sp.cosh(2*eta))
w = sp.sinh(eta)
out['T1_is_BL1_first_cycle'] = z(sp.cosh(2*eta) - (1 + 2*w**2))          # BL1-T2: 1 + 2 w^2, sinh eta = w
# T2: ledger of the step (GE2 (10)): h = (1/2) sum w |v|^2 = eta^2/2 ; small steps: d Log(mean S)/d tau = 4 = (degree)^2
h = eta**2/2
out['T2_rate'] = sp.limit((sp.cosh(2*eta) - 1)/h, eta, 0) == 4
deg = sp.symbols('n', positive=True)
out['T2_degree_squared'] = sp.limit((sp.cosh(deg*eta) - 1)/h, eta, 0) == deg**2
# T3: the gap. Centre value M ~ r (degree 1): mean count / count of the mean centre = 1 + tanh^2 eta
gap = sp.cosh(2*eta)/sp.cosh(eta)**2
out['T3_gap'] = z(gap - (1 + sp.tanh(eta)**2))
# T4: two centres exchanging at fixed total (flat count 4 pi M^2 each): any unbiased exchange raises the count
M1, M2 = sp.symbols('M1 M2', positive=True)
Sc = lambda a, b: 4*sp.pi*(a**2 + b**2)
gain = sp.simplify((Sc(M1 + d_, M2 - d_) + Sc(M1 - d_, M2 + d_))/2 - Sc(M1, M2))
out['T4_exchange_gain'] = sp.simplify(gain - 8*sp.pi*d_**2) == 0
# equal centres are the least count at fixed total, not the greatest
t_ = sp.symbols('t', real=True)
out['T4_equal_is_minimum'] = sp.solve(sp.diff(Sc(M/2 + t_, M/2 - t_), t_), t_) == [0] and sp.diff(Sc(M/2 + t_, M/2 - t_), t_, 2) > 0
# T5: with expansion, frame-consistent count along fixed h: K'(r) = 2 pi r [1/c + g x (-c')/(2 c^2)], x = h r
x = sp.symbols('x', positive=True)
u = x*(1 - x**2)/2; c = sp.sqrt(1 - 3*u**R_(2, 3)); g = 2*(1 - x**2)/(1 - 3*x**2)
dens = 1/c + g*x*(-sp.diff(c, x))/(2*c**2)                 # K'(r) / (2 pi r)
# drift under symmetric steps in Log r:  (r d/dr)^2 K = 2 pi r^2 [2 dens + x dens'(x)]
drift = sp.lambdify(x, 2*dens + x*sp.diff(dens, x), 'mpmath')
pts = [mp.mpf(k_)/100 for k_ in (1, 5, 10, 20, 30, 40, 50, 55, 57)]
out['T5_drift_positive'] = all(drift(p) > 0 for p in pts)
out['T5_flat_limit'] = abs(drift(mp.mpf(10)**-12) - 2) < mp.mpf(10)**-6     # -> 2, i.e. (r d/dr)^2 S = 4 S
# T6: in two directions (Log r_b, Log h) only the symmetric part of the count form enters the mean; XU3's circuit part averages out
a_, b_, c_ = sp.symbols('a b c', real=True)
v1, v2 = sp.symbols('v1 v2', real=True)
T = sp.Matrix([[a_, b_ + c_], [b_ - c_, d_]])               # gradient of a one-form: symmetric part b, circuit part c
vv = sp.Matrix([v1, v2])
out['T6_only_symmetric_part'] = sp.simplify((vv.T*T*vv)[0] - (a_*v1**2 + 2*b_*v1*v2 + d_*v2**2)) == 0
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SL1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
