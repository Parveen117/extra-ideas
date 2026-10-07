"""CV1: the curvature of the free-fall frame field, as the order defect of its directions.

Frame of the thesis (MO1, GB1), c = 1:   e0 = d_t - beta(r) d_r ,  e1 = d_r ,  e2 = (1/r) d_theta ,
                                          e3 = (1/(r sin theta)) d_phi ;   invariant form eta = diag(+,-,-,-) (IN1).
Native route (EMK-C1 form):
    order defect of the directions     [e_a, e_b] = c_ab^c e_c
    the connection that keeps eta and has no torsion, from the c's alone
    curvature = order defect of the connected directions:
        R^a_b(c,d) = e_c w^a_bd - e_d w^a_bc + w^a_ec w^e_bd - w^a_ed w^e_bc - c_cd^e w^a_be
Cross-check: the same curvature from the interval dtau^2 = dt^2 - (dr + beta dt)^2 - r^2 dOmega^2 by the
coordinate route.  With m = beta^2 (the memory):
T1  the two routes agree, and every curvature component depends on beta only through m, m', m''.
T2  the contracted curvature vanishes  iff  (r m)' = 0  iff  m = r_s / r.
T3  for m = r_s / r the uncontracted curvature is (r_s / 2 r^3) x (2, -1, -1) on the planes (time, direction):
    it is not zero, and no change of frame at a point removes it (its invariant 12 r_s^2 / r^6).
T4  (r m)' = L r^2 gives contracted curvature = constant x eta: the constant-curvature case, m = r_s/r + L r^2/3.
T5  the time-time contracted curvature is exactly -(1/2) x MC1's variance operator (1/r^2)(r^2 m')';
    least variance alone allows m = A/r + B; the remaining components (r m)'/r^2 remove the constant B.
Symbolic algebra with sympy (not standard library).  Python 3.12.
"""
import itertools
import json

import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi', positive=True)
X = (t, r, th, ph)
beta = sp.Function('beta')(r)
ETA = sp.diag(1, -1, -1, -1)


def frame(b):
    """rows: components of e_a in the coordinate basis."""
    return sp.Matrix([[1, -b, 0, 0], [0, 1, 0, 0], [0, 0, 1/r, 0], [0, 0, 0, 1/(r*sp.sin(th))]])


def apply(vec, f):
    return sum(vec[i]*sp.diff(f, X[i]) for i in range(4))


def structure(E):
    Einv = E.inv()
    c = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            comm = [apply(E.row(a), E[b, i]) - apply(E.row(b), E[a, i]) for i in range(4)]
            for k in range(4):
                c[a][b][k] = sp.simplify(sum(comm[i]*Einv[i, k] for i in range(4)))
    return c


def connection(c):
    """w[a][b][k] = w^a_{b k}: connected derivative along e_k of e_b has component a."""
    low = lambda a, b, k: ETA[k, k]*c[a][b][k]          # c_{ab k} with last index lowered
    w = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for k in range(4):
                # Koszul for an eta-orthonormal frame: Gamma_{a b k} = <nabla_k e_b, e_a>
                val = sp.Rational(1, 2)*(low(k, b, a) - low(b, a, k) + low(a, k, b))
                w[a][b][k] = sp.simplify(ETA[a, a]*val)
    return w


def curvature(E, c, w):
    R = {}
    for a, b, k, l in itertools.product(range(4), repeat=4):
        val = apply(E.row(k), w[a][b][l]) - apply(E.row(l), w[a][b][k])
        val += sum(w[a][e][k]*w[e][b][l] - w[a][e][l]*w[e][b][k] for e in range(4))
        val -= sum(c[k][l][e]*w[a][b][e] for e in range(4))
        R[(a, b, k, l)] = sp.simplify(val)
    return R


def contracted(R):
    return sp.Matrix(4, 4, lambda b, l: sp.simplify(sum(R[(a, b, a, l)] for a in range(4))))


def coordinate_route(b):
    g = sp.Matrix([[1 - b**2, -b, 0, 0], [-b, -1, 0, 0], [0, 0, -r**2, 0], [0, 0, 0, -r**2*sp.sin(th)**2]])
    gi = g.inv()
    G = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, bb], X[cc]) + sp.diff(g[d, cc], X[bb]) - sp.diff(g[bb, cc], X[d])) for d in range(4))/2)
           for cc in range(4)] for bb in range(4)] for a in range(4)]
    def riem(a, bb, cc, d):
        v = sp.diff(G[a][bb][d], X[cc]) - sp.diff(G[a][bb][cc], X[d])
        v += sum(G[a][e][cc]*G[e][bb][d] - G[a][e][d]*G[e][bb][cc] for e in range(4))
        return v
    ric = sp.Matrix(4, 4, lambda bb, d: sp.simplify(sum(riem(a, bb, a, d) for a in range(4))))
    return g, ric


def run():
    m = sp.Function('m')(r)
    E = frame(beta)
    c = structure(E)
    w = connection(c)
    R = curvature(E, c, w)
    Ric = contracted(R)
    # T1: express through m = beta^2
    def through_m(expr):
        e = sp.simplify(expr)
        e = e.subs(sp.diff(beta, r, 2), (sp.diff(m, r, 2)/2 - sp.diff(beta, r)**2)/beta)
        e = e.subs(sp.diff(beta, r), sp.diff(m, r)/(2*beta))
        e = sp.simplify(e.subs(beta**2, m))
        return sp.simplify(e.subs(beta, sp.sqrt(m)))
    Ric_m = Ric.applyfunc(through_m)
    if any(x.has(beta) for x in Ric_m):
        raise ValueError('contracted curvature must depend on beta through m only')
    # coordinate route cross-check: frame components of the coordinate Ricci tensor
    g, ric_coord = coordinate_route(beta)
    ric_frame = sp.simplify(E*ric_coord*E.T)
    diff = (ric_frame - Ric).applyfunc(sp.simplify)
    if diff != sp.zeros(4, 4):
        raise ValueError('the two routes disagree')
    rm1 = sp.diff(r*m, r)
    rm2 = sp.diff(r*m, r, 2)
    expected = sp.diag(rm2/(2*r), -rm2/(2*r), -rm1/r**2, -rm1/r**2)
    signs = None
    for s in (1, -1):
        if (Ric_m - s*expected).applyfunc(sp.simplify) == sp.zeros(4, 4):
            signs = s
    if signs is None:
        raise ValueError('contracted curvature is not of the form ((r m)\'\'/2r, (r m)\'/r^2)')
    # T5a: the time-time contracted curvature is MC1's variance operator, -(1/2) (1/r^2)(r^2 m')'
    lap = sp.diff(r**2*sp.diff(m, r), r)/r**2
    if sp.simplify(Ric_m[0, 0] + lap/2) != 0:
        raise ValueError('time-time contracted curvature is not minus half the radial Laplacian of the memory')
    # T2
    rs, L, A, B = sp.symbols('r_s L A B', positive=True)
    vac = Ric_m.subs(m, rs/r).doit().applyfunc(sp.simplify)
    if vac != sp.zeros(4, 4):
        raise ValueError('m = r_s/r must have zero contracted curvature')
    const = Ric_m.subs(m, B).doit().applyfunc(sp.simplify)
    if const == sp.zeros(4, 4):
        raise ValueError('a constant memory must not be curvature-free')
    # T3: uncontracted curvature for m = r_s/r, frame components R_{0i0i}
    Rv = {k: sp.simplify(v.subs(sp.diff(beta, r, 2), sp.diff(sp.sqrt(rs/r), r, 2)).subs(sp.diff(beta, r), sp.diff(sp.sqrt(rs/r), r)).subs(beta, sp.sqrt(rs/r))) for k, v in R.items()}
    tidal = [sp.simplify(Rv[(0, i, 0, i)]) for i in (1, 2, 3)]
    low = lambda a, b, k, l: ETA[a, a]*Rv[(a, b, k, l)]
    up = lambda a, b, k, l: ETA[b, b]*ETA[k, k]*ETA[l, l]*Rv[(a, b, k, l)]
    square = sp.simplify(sum(low(*idx)*up(*idx) for idx in itertools.product(range(4), repeat=4)))
    if sp.simplify(sum(tidal)) != 0 or sp.simplify(square - 12*rs**2/r**6) != 0:
        raise ValueError('vacuum curvature pattern failed')
    # T4
    lam = Ric_m.subs(m, rs/r + L*r**2/3).doit().applyfunc(sp.simplify)
    ratio = [sp.simplify(lam[i, i]/ETA[i, i]) for i in range(4)]
    if len(set(ratio)) != 1 or ratio[0] == 0:
        raise ValueError('constant-curvature case failed')
    # T5
    mc1 = sp.simplify(sp.diff(r**2*sp.diff(A/r + B, r), r))
    return dict(contracted=[str(Ric_m[i, i]) for i in range(4)], overall_sign=signs,
                vacuum_tidal=[str(x) for x in tidal], curvature_square=str(square),
                constant_curvature_ratio=str(ratio[0]), constant_memory_contracted=[str(const[i, i]) for i in range(4)],
                mc1_operator_on_A_over_r_plus_B=str(mc1),
                structure_nonzero={f'{a}{b}{k}': str(c[a][b][k]) for a in range(4) for b in range(a+1, 4) for k in range(4) if c[a][b][k] != 0})


if __name__ == '__main__':
    res = run()
    json.dump(res, open('CV1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
