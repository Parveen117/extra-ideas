"""CD1: closure and the number of space dimensions.

MC1: the least-cost memory in d dimensions is B r^-(d-2) (log r for d = 2).  MO1: the form with any memory m(r):
    dtau^2 = dt^2 - (dr + beta dt)^2 - r^2 dphi^2 ,  beta^2 = m ,   static form N^2 = 1 - m (MO1-M1), radial law (MO1-M3).
RC1: a bound history closes iff (in-out rate)/(round rate) is rational.
Symbolic (sympy).  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


def ricci_static(d, m, r):
    """coordinate route: static form with memory m(r) in d space dimensions; returns the contracted curvature (diagonal)."""
    t = sp.Symbol('t')
    ang = sp.symbols('th1:%d' % d, positive=True)
    X = (t, r) + tuple(ang)
    N2 = 1 - m
    diag = [N2, -1/N2]
    f = r**2
    for th in ang:
        diag.append(-f)
        f = f*sp.sin(th)**2
    n = len(X)
    g = sp.diag(*diag); gi = sp.diag(*[1/v for v in diag])
    G = [[[sp.simplify(sum(gi[i, l]*(sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(n))/2)
           for k in range(n)] for j in range(n)] for i in range(n)]
    def ric(j, k):
        v = sum(sp.diff(G[i][j][k], X[i]) - sp.diff(G[i][j][i], X[k]) for i in range(n))
        v += sum(G[i][i][l]*G[l][j][k] - G[i][k][l]*G[l][j][i] for i in range(n) for l in range(n))
        return sp.simplify(v)
    return [sp.simplify(ric(j, j)/g[j, j]) for j in range(n)]


def run():
    out = {}
    z = lambda e: sp.simplify(e) == 0
    u, L, mu, n_ = sp.symbols('u L mu n', positive=True)
    mf = sp.Function('m')
    # D1: orbit equation for any memory m(u), u = 1/r, from MO1-M3:  (u')^2 = ( E^2 - (1 - m)(1 + L^2 u^2) ) / L^2
    E = sp.Symbol('E', positive=True)
    rhs2 = (E**2 - (1 - mf(u))*(1 + L**2*u**2))/L**2
    force = sp.diff(rhs2, u)/2                                     # u'' = force(u)
    out['D1_orbit_equation'] = z(force - (-u + sp.diff(mf(u), u)*(1 + L**2*u**2)/(2*L**2) + mf(u)*u))
    # circle: force = 0 fixes L; (in-out / round)^2 = - d force/du there
    Lsol = sp.solve(sp.Eq(force, 0), L**2)[0]
    ratio2 = sp.simplify((-sp.diff(force, u)).subs(L**2, Lsol))
    m0, m1, m2 = mf(u), sp.diff(mf(u), u), sp.diff(mf(u), u, 2)
    general = 1 - m0 - 2*u*m1 - u*(1 - m0)*m2/m1
    out['D2_ratio_for_any_memory'] = z(ratio2 - general)
    # D3: least-cost memory of MC1 in d dimensions, m = mu u^(d-2):  ratio^2 = (4 - d) - d m
    for d in (3, 4, 5, 6):
        mm = mu*u**(d - 2)
        val = sp.simplify(general.subs(mf(u), mm).doit())
        out['D3_power_law_d%d' % d] = z(val - ((4 - d) - d*mm))
    dd = sp.Symbol('d', positive=True)
    mm = mu*u**(dd - 2)
    out['D3_power_law_any_d'] = z(sp.simplify(general.subs(mf(u), mm).doit()) - ((4 - dd) - dd*mm))
    # d = 2: m = c + mu log u :  ratio^2 = 2 - 2 m - 2 mu
    c0 = sp.Symbol('c0', positive=True)
    ml = c0 + mu*sp.log(u)
    out['D3_log_d2'] = z(sp.simplify(general.subs(mf(u), ml).doit()) - (2 - 2*ml - 2*mu))
    # D4: first order (m -> 0): ratio = sqrt(4 - d): 1 in d = 3 (closes), sqrt 2 in d = 2 (never), 0 in d = 4, none beyond
    top = {d: sp.sqrt(4 - d) for d in (2, 3, 4, 5, 6)}
    out['D4_only_three_closes'] = [d for d in top if top[d].is_rational and top[d] > 0] == [3]
    out['D4_two_is_irrational'] = top[2].is_rational is False
    out['D4_four_and_more_have_no_stable_circle'] = all(((4 - d) - d*sp.Rational(1, 100)) < 0 for d in (4, 5, 6))
    # D5: carried direction / count on a free circle: (dtau/dt)^2 = 1 - m - r^2 Omega^2 , r^2 Omega^2 = u m'/2  ->  1 - (d/2) m
    count2 = 1 - mm - u*sp.diff(mm, u)/2
    out['D5_count_factor'] = z(sp.simplify(count2) - (1 - dd*mm/2))
    # special shells in d dimensions: last stable circle m = (4-d)/d, light circle m = 2/d, clock stops m = 1;
    # in d = 3 they are 1/3, 2/3, 1
    shells = lambda d: (sp.Rational(4 - d, d), sp.Rational(2, d), sp.Integer(1))
    out['D5_thirds_in_three'] = shells(3) == (sp.Rational(1, 3), sp.Rational(2, 3), 1)
    # in every d >= 3 the three shells are equally spaced in memory, step (d-2)/d; the first lies at positive memory only for d = 3
    out['D5_equal_steps_in_every_d'] = all(shells(d)[1] - shells(d)[0] == shells(d)[2] - shells(d)[1] == sp.Rational(d - 2, d) for d in range(3, 12))
    out['D5_first_shell_exists_only_in_three'] = [d for d in range(3, 12) if shells(d)[0] > 0] == [3]
    # D6: the same exponent from the frame law: static form in d = 3, 4, 5: contracted curvature zero iff (r^(d-2) m)' = 0
    r = sp.Symbol('r', positive=True)
    B = sp.Symbol('B', positive=True)
    ok = True
    for d in (3, 4, 5):
        ric = ricci_static(d, B*r**(-(d - 2)), r)
        ok &= all(v == 0 for v in ric)
        ric_other = ricci_static(d, B*r**(-(d - 1)), r)
        ok &= any(v != 0 for v in ric_other)
    out['D6_frame_law_gives_the_same_exponent'] = ok
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out), open(os.path.join(HERE, 'CD1_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
