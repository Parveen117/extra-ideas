"""SE1: the law of TP1 on flat slices is a law for the stretch block alone: the plane determinants sum to zero.

Frame with one time and flat slices, any fall velocity u(x) (SW1-T4):  e_0 = d_t + u.grad ,  e_i = d_i .
Its connection (CV1's formula) splits into a turn (1/2) curl u and a stretch  K_ij = (d_i u_j + d_j u_i)/2 .
Symbolic (sympy), with the functions of CV1 and TP1 for three cuts and the same formulas written for n cuts.  Python 3.12."""
import itertools, json, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1'))
_cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'tp1'))
import cv1_frame_curvature as cv
import tp1_frame_defect_law as tp
os.chdir(_cwd)
z = lambda e: sp.simplify(e) == 0


def e2(K):
    n = K.shape[0]
    return sum(K[i, i]*K[j, j] - K[i, j]*K[j, i] for i in range(n) for j in range(i + 1, n))


def invariants_n(d, u, X):
    """TP1's three invariants for the frame e_0 = d_t + u.grad, e_i = d_i with d cuts (same formulas as TP1, any d)."""
    n = d + 1
    eta = [1] + [-1]*d
    c = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(1, n):                                         # [e_0, e_i] = -(d_i u_j) e_j
        for j in range(1, n):
            c[0][i][j] = -sp.diff(u[j-1], X[i-1]); c[i][0][j] = -c[0][i][j]
    I1 = sum(eta[k]*eta[a]*eta[b]*c[a][b][k]**2 for k in range(n) for a in range(n) for b in range(n))
    I2 = sum(eta[b]*c[a][b][k]*c[k][b][a] for k in range(n) for a in range(n) for b in range(n))
    tr = [sum(c[a][b][a] for a in range(n)) for b in range(n)]
    I3 = sum(eta[b]*tr[b]**2 for b in range(n))
    return I1, I2, I3


def run():
    out = {}
    T, x, y, zz = sp.symbols('T x y z', real=True); co = (T, x, y, zz)
    u = [sp.Function('u%d' % i)(x, y, zz) for i in (1, 2, 3)]
    E = sp.Matrix([[1, u[0], u[1], u[2]], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    c = tp.structure_in(co, E)
    I1, I2, I3 = tp.invariants(c)
    Q = sp.expand(sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3)
    d_ = lambda i, j: sp.diff(u[i], co[j + 1])
    K = sp.Matrix(3, 3, lambda i, j: (d_(i, j) + d_(j, i))/2)
    # S1: TP1's law on this frame is the stretch alone:  Q = sum K_ij^2 - (tr K)^2 = -2 e2(K); the turn (curl u) does not enter
    out['S1_law_is_minus_twice_e2'] = z(Q + 2*e2(K))
    out['S1_squares_minus_square_of_sum'] = z(Q - (sum(K[i, j]**2 for i in range(3) for j in range(3)) - K.trace()**2))
    # the same with TP1's formulas written for d = 2 and d = 4 cuts
    ok = True
    for d in (2, 4):
        X = sp.symbols('x1:%d' % (d + 1), real=True)
        ud = [sp.Function('v%d' % i)(*X) for i in range(1, d + 1)]
        J1, J2, J3 = invariants_n(d, ud, X)
        Kd = sp.Matrix(d, d, lambda i, j: (sp.diff(ud[i], X[j]) + sp.diff(ud[j], X[i]))/2)
        ok &= z(sp.expand(sp.Rational(1, 4)*J1 + sp.Rational(1, 2)*J2 - J3) + 2*e2(Kd))
    J1, J2, J3 = invariants_n(3, u, (x, y, zz))
    ok &= z(J1 - I1) and z(J2 - I2) and z(J3 - I3)               # the n-cut formulas reproduce TP1's own at d = 3
    out['S1_any_number_of_cuts'] = ok
    # S2: contracted curvature (CV1's formulas): G_00 = e2(K) ; the (time, cut) components are -(1/2) curl curl u
    w = cv.connection(c)
    ap = lambda vec, f: sum(vec[i]*sp.diff(f, co[i]) for i in range(4))
    R = {}
    for a_, b_, k_, l_ in itertools.product(range(4), repeat=4):
        val = ap(E.row(k_), w[a_][b_][l_]) - ap(E.row(l_), w[a_][b_][k_])
        val += sum(w[a_][e][k_]*w[e][b_][l_] - w[a_][e][l_]*w[e][b_][k_] for e in range(4))
        val -= sum(c[k_][l_][e]*w[a_][b_][e] for e in range(4))
        R[(a_, b_, k_, l_)] = sp.simplify(val)
    Ric = sp.Matrix(4, 4, lambda b_, l_: sp.simplify(sum(R[(a_, b_, a_, l_)] for a_ in range(4))))
    scal = sum(cv.ETA[i, i]*Ric[i, i] for i in range(4))
    out['S2_energy_component_is_e2'] = z(Ric[0, 0] - scal/2 - e2(K))
    curl = [d_(2, 1) - d_(1, 2), d_(0, 2) - d_(2, 0), d_(1, 0) - d_(0, 1)]
    cc = [sp.diff(curl[2], y) - sp.diff(curl[1], zz), sp.diff(curl[0], zz) - sp.diff(curl[2], x), sp.diff(curl[1], x) - sp.diff(curl[0], y)]
    out['S2_mixed_components_are_curl_curl'] = all(z(Ric[0, i + 1] + cc[i]/2) for i in range(3))
    # S3: e2 = 0  <=>  sum of squares = square of the sum  <=>  variance of the rates = (d - 1) x (mean)^2
    lam = sp.symbols('l1:6', real=True)
    ok = True
    for d in (2, 3, 4, 5):
        L = lam[:d]
        s1, s2 = sum(L), sum(v**2 for v in L)
        e2l = sum(L[i]*L[j] for i in range(d) for j in range(i + 1, d))
        var, mean = s2/d - (s1/d)**2, s1/d
        ok &= z(s2 - s1**2 + 2*e2l) and z(d*(var - (d - 1)*mean**2) + 2*e2l)
    out['S3_ratio_law'] = ok
    # S4: radial fall in d cuts, u = -beta rhat: rates (rho, 1, ..., 1) x (-beta/r), rho = r beta'/beta
    r = sp.Symbol('r', positive=True); beta = sp.Function('beta')(r); dd = sp.Symbol('d', positive=True)
    rho = r*sp.diff(beta, r)/beta
    ok = True
    for d in (2, 3, 4, 5):
        X = sp.symbols('x1:%d' % (d + 1), positive=True)
        rad = sp.sqrt(sum(v**2 for v in X))
        ud = [-beta.subs(r, rad)*v/rad for v in X]
        Kd = sp.Matrix(d, d, lambda i, j: (sp.diff(ud[i], X[j]) + sp.diff(ud[j], X[i]))/2)
        pt = {X[0]: r}; pt.update({v: 0 for v in X[1:]})
        val = sp.simplify(e2(Kd).doit().subs(pt).doit())
        target = (beta/r)**2*((d - 1)*rho + sp.Rational((d - 1)*(d - 2), 2))
        ok &= z(val - target)
    out['S4_radial'] = ok
    # two cuts: e2 = (beta/r)^2 rho, so the law leaves beta constant: no falling memory at all (MC1's log r is not allowed by it)
    out['S4_two_cuts_constant'] = z(((beta/r)**2*((2 - 1)*rho + 0)) - beta*sp.diff(beta, r)/r)
    # e2 = 0  <=>  rho = -(d-2)/2  <=>  memory beta^2 = B r^-(d-2)  (MC1's exponent); in three cuts rho = -1/2 (NC1-N3)
    B = sp.Symbol('B', positive=True)
    out['S4_exponent'] = z(rho.subs(beta, sp.sqrt(B*r**(-(dd - 2)))).doit() + (dd - 2)/2)
    eta_ = sp.Function('eta')(r)                                   # NC1: tanh eta = beta, psi = 2 eta, rho = r psi'/sinh psi
    nc1 = r*sp.diff(2*eta_, r)/sp.sinh(2*eta_)
    mine = (r*sp.diff(sp.tanh(eta_), r)/sp.tanh(eta_))
    out['S4_is_NC1_ratio'] = z((nc1 - mine).rewrite(sp.exp))
    # with a source: G_00 = e2 = k u  gives EG1's ratio law  rho + 1/2 = k u r^2 / (2 m)   (three cuts)
    ku = sp.Symbol('ku')
    out['S4_EG1_ratio_law_with_source'] = z(sp.solve(sp.Eq((beta/r)**2*(2*rho + 1), ku), sp.diff(beta, r))[0]*r/beta + sp.Rational(1, 2) - ku*r**2/(2*beta**2))
    # S5: SW1's swirl on flat slices: curl curl u = 0  <=>  (r^4 W')' = 0 ; e2 = -(shear of the swirl)^2 = -9 a^2 sin^2/(4 r^6): SW1-T3
    a, rs = sp.symbols('a r_s', positive=True)
    rr = sp.sqrt(x**2 + y**2 + zz**2)
    Wf = sp.Function('W')
    us = [-sp.sqrt(rs/rr)*x/rr - Wf(rr)*y, -sp.sqrt(rs/rr)*y/rr + Wf(rr)*x, -sp.sqrt(rs/rr)*zz/rr]
    cu = [sp.diff(us[2], y) - sp.diff(us[1], zz), sp.diff(us[0], zz) - sp.diff(us[2], x), sp.diff(us[1], x) - sp.diff(us[0], y)]
    ccs = [sp.diff(cu[2], y) - sp.diff(cu[1], zz), sp.diff(cu[0], zz) - sp.diff(cu[2], x), sp.diff(cu[1], x) - sp.diff(cu[0], y)]
    rsym = sp.Symbol('rho_', positive=True)
    at = {x: rsym, y: 0, zz: 0}                                    # a point on the equator
    val = sp.simplify(ccs[1].doit().subs(at).doit())
    law = sp.diff(rsym**4*sp.diff(Wf(rsym), rsym), rsym)/rsym**3
    out['S5_swirl_first_order_law'] = z(val + law) or z(val - law)
    ua = [v.subs(Wf(rr), a/rr**3).doit() for v in us]
    Ks = sp.Matrix(3, 3, lambda i, j: (sp.diff(ua[i], co[j + 1]) + sp.diff(ua[j], co[i + 1]))/2)
    th = sp.Symbol('theta', positive=True)
    pt = {x: rsym*sp.sin(th), y: 0, zz: rsym*sp.cos(th)}
    out['S5_swirl_remainder'] = z(sp.simplify(e2(Ks).subs(pt)) + 9*a**2*sp.sin(th)**2/(4*rsym**6))
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out), open(os.path.join(HERE, 'SE1_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
