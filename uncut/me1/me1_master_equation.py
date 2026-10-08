"""ME1: one form for all planes, and its defect when the flatness conditions are dropped. Exact sympy."""
import json, os, itertools, sympy as sp
S, Q, J, p = sp.symbols('S Q J p', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
# a centre with charge and spin (BH1 + BH2 combined), then read by a moving reader (MT1)
M2 = (S/(4*sp.pi))*(1 + sp.pi*Q**2/S)**2 + sp.pi*J**2/S
M = sp.sqrt(M2)
out['T0_reduces'] = z(M.subs(J, 0) - sp.sqrt(S/sp.pi)/2*(1 + sp.pi*Q**2/S)) and z(M.subs(Q, 0)**2 - (S/(4*sp.pi) + sp.pi*J**2/S))
E = sp.sqrt(M2 + p**2)                         # energy read by the mover; four variables, four planes
T, Phi, Om, v = [sp.diff(E, x) for x in (S, Q, J, p)]
T0, Phi0, Om0 = [sp.diff(M, x) for x in (S, Q, J)]
gam = E/M
# T1: every partner is lowered by the same factor; the velocity is the fourth partner
out['T1_partners'] = z(T - T0/gam) and z(Phi - Phi0/gam) and z(Om - Om0/gam) and z(v - p/E)
# T2: the angles swept per closed turn do not depend on the reader
P = 1/T; P0 = 1/T0
out['T2_angles_reader_free'] = z(Phi*P - Phi0*P0) and z(Om*P - Om0*P0)
# T3: the one form:  dS = P (dE - Phi dQ - Omega dJ - v dp) ; all six pair equations hold (one potential)
vars_ = (S, Q, J, p); parts = (T, Phi, Om, v)
out['T3_six_pair_equations'] = all(z(sp.diff(parts[i], vars_[j]) - sp.diff(parts[j], vars_[i])) for i, j in itertools.combinations(range(4), 2))
# T4: drop the flatness conditions. theta = k * theta0, theta0 = sum a_i dq_i (no potential assumed)
q = sp.symbols('q1:5'); k = sp.Function('k')(*q); a = [sp.Function(f'a{i}')(*q) for i in range(1, 5)]
ok = True
for i, j in itertools.combinations(range(4), 2):
    full = sp.diff(k*a[j], q[i]) - sp.diff(k*a[i], q[j])
    weight = sp.diff(k, q[i])*a[j] - sp.diff(k, q[j])*a[i]
    frame = k*(sp.diff(a[j], q[i]) - sp.diff(a[i], q[j]))
    ok = ok and z(full - weight - frame)
out['T4_master_equation'] = ok                 # d(theta) = dk ^ theta0 + k d(theta0), component by component
# T5: count of flatness conditions for n planes: n(n-1)/2 pair equations + n conditions on k
n = sp.symbols('n', positive=True, integer=True)
cnt = n*(n - 1)/2 + n
out['T5_counts'] = [int(cnt.subs(n, m)) for m in (2, 3, 4, 5)] == [3, 6, 10, 15]
out['T5_two_planes_is_UP5'] = int(cnt.subs(n, 2)) == 3      # UP5-T5: (alpha, beta, F)
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ME1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
