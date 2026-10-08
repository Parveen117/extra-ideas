"""FD1 against the line's existing response curvature (QC1-T2, in QC2's entropy representation), van der Waals.
Exact sympy."""
import json, os, sympy as sp
u, v, a, b, cv, T, eps = sp.symbols('u v a b c_v T epsilon', positive=True)
out = {}
def curvature(s, X):
    H = -sp.hessian(s, X); Hi = H.inv()
    A = [sp.diff(H, x) for x in X]; B = [sp.diff(Hi, x) for x in X]
    return sp.simplify(sp.Rational(1, 4)*Hi*(A[0]*B[1] - A[1]*B[0])*H)     # QC1-T2
X = (u, v)
F = curvature(cv*sp.log(u + a/v) + sp.log(v - b), X)                        # van der Waals, Nk = 1
f2 = sp.simplify(F.det().subs(u, cv*T - a/v))                               # F = f J_H, det F = f^2
Dp = 2*a*(v - b)/(T*v**2)                                                   # FD1-T4
out['T1_same_source'] = sp.simplify(f2.subs(a, 0)) == 0 and sp.simplify(Dp.subs(a, 0)) == 0      # both vanish without attraction; hard cores alone are flat
phi = sp.Function('phi')
out['T2_ideal_any_heat_capacity_flat'] = curvature(phi(u) + sp.log(v), X) == sp.zeros(2)         # QC2 section 2; FD1-T6
# low density: f du dv = [D'/(2 sqrt(c_v))] dlnT dlnv : f * (c_v T) * v -> D'/(2 sqrt c_v)
ratio2 = sp.simplify(f2*(cv*T*v)**2/(Dp**2/(4*cv)))
out['T3_low_density_relation'] = sp.limit(ratio2.subs(v, 1/eps), eps, 0) == 1
out['T3_not_identical_at_finite_density'] = sp.simplify(ratio2 - 1) != 0
out = {k: bool(val) for k, val in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FD1_RELATION_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
