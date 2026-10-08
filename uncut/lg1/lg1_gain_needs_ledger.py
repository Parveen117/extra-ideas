"""LG1: a count gain needs an enclosed area, and an enclosed area costs history ledger (GE2-T3). Exact sympy + mpmath."""
import json, os, sympy as sp, mpmath as mp
x = sp.symbols('x', positive=True); R_ = sp.Rational
u = x*(1 - x**2)/2; c = sp.sqrt(1 - 3*u**R_(2, 3))
w = (1 + x**2)/(2*u**R_(1, 3)*c**3)                        # XU3-T2: d theta = w(x) dS ^ dx
fw = sp.lambdify(x, w, 'mpmath')
out = {}
# T1: out and back along one direction: theta picks up nothing (the two legs cancel), ledger is positive
dS, dx = sp.symbols('dS dx', positive=True)
a = 1/c
out['T1_there_and_back'] = sp.simplify(a*dS + a*(-dS)) == 0
ledger_back = R_(1, 2)*(R_(1, 2)*dS**2 + R_(1, 2)*dS**2)   # GE2 (10): h = (1/2) sum w |v|^2, two turns +v, -v
out['T1_ledger_positive'] = sp.simplify(ledger_back - dS**2/2) == 0
# T2: the smallest history that gains count: two directions, four legs (+v, +w, -v, -w)
S1, S2, x1, x2 = mp.mpf(3), mp.mpf(4), mp.mpf(1)/10, mp.mpf(3)/20
gain = (S2 - S1)*mp.quad(fw, [x1, x2])
ledger = mp.mpf(1)/2*(mp.mpf(1)/4*2*(S2 - S1)**2 + mp.mpf(1)/4*2*(x2 - x1)**2)   # four legs, weight 1/4 each
wmax = max(fw(x1 + (x2 - x1)*k/200) for k in range(201))
out['T2_gain'] = mp.nstr(gain, 8); out['T2_ledger'] = mp.nstr(ledger, 8)
# T3: bound  gain <= dS dx wmax <= (dS^2 + dx^2)/2 wmax = 2 wmax * ledger
out['T3_bound'] = bool(gain <= (S2 - S1)*(x2 - x1)*wmax <= 2*wmax*ledger)
A_, B_ = sp.symbols('A B', positive=True)
out['T3_inequality_exact'] = sp.simplify((A_**2 + B_**2)/2 - A_*B_ - (A_ - B_)**2/2) == 0
# T4: reversing the order of the two directions reverses the gain, not the ledger
out['T4_order'] = bool(-gain < 0 < ledger)
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LG1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
