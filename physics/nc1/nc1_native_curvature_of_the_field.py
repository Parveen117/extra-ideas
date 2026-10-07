"""NC1: the framework's OWN curvature (NT-3) of the gravity field, and the laws it offers.

Native definition (Publications, native-thermodynamic-curvature, NT-2/NT-3):
    response element H (self-dagger),  X_i = H^-1 d_i H,  connection A = q X,
    F^{qX}_ij = q (q - 1) [X_i, X_j] ;   q = 1 is pure gauge (flat);  the strain choice q = 1/2 gives
    F^H_ij = -(1/4) [X_i, X_j] .
Field of the thesis (GB1): the self-dagger unit block  exp(eta(r) n.C),  tanh(eta) = beta,  m = beta^2.
Its response element (pairing) is G = exp(psi n.C) with psi = 2 eta.
N1  NT-3 re-certified on this field: F^{qX} = q(q-1)[X_i, X_j] for q = 1/2, and q = 1 is flat.
N2  The native curvature has two invariant sizes (sc of F^2, frame-independent):
        radial-angular   -(psi' sinh psi)^2 / (4 r^2)        angular-angular   -(sinh psi)^4 / (4 r^4)
    and their ratio is the pure number   rho = r psi' / sinh psi .
N3  RATIO LAW.  rho = -1/2 at every radius   <=>   tanh^2(eta) = r_s / r   exactly: the field of CV1/MO1.
N4  VARIATIONAL LAW.  Stationary total invariant curvature-square gives  w'' = w (w^2 - 1) / r^2 , w = cosh psi.
    Linearised: w - 1 = A/r + B r^2  (the 1/r field and the constant-curvature term).
    Beyond first order: w = 1 + A/r + (3/4) A^2/r^2 + ... ; the field of N3 does NOT solve it.
N5  Consequence for the orbit advance (MO1's radial law, m = r_s/r + kappa r_s^2/r^2 gives (3 + 2 kappa) pi r_s/p):
        ratio law N3:        kappa = 0      ->  42.98" per century for Mercury
        variational law N4:  kappa = +1/2   ->  57.3"      (response element psi = 2 eta)
                             kappa = -3/8   ->  32.2"      (if the unit block itself, psi = eta, is used)
Computation in the 2x2 representation of three anticommuting self-dagger cuts (FR1) with sympy.  Python 3.12.
"""
import json

import sympy as sp

x, y = sp.symbols('x y', real=True)
z = sp.symbols('z', positive=True)
r = sp.symbols('r', positive=True)
XS = (x, y, z)
C = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
ONE = sp.eye(2)
psi = sp.Function('psi')


def sc(M):
    return sp.simplify(M.trace()/2)


def field(q):
    """raw (unsimplified) X_i and F_ij of the connection q X; simplify only after evaluation on the axis."""
    rad = sp.sqrt(x**2 + y**2 + z**2)
    n = (x/rad)*C[0] + (y/rad)*C[1] + (z/rad)*C[2]
    G = sp.cosh(psi(rad))*ONE + sp.sinh(psi(rad))*n
    Ginv = sp.cosh(psi(rad))*ONE - sp.sinh(psi(rad))*n
    X = [Ginv*sp.diff(G, XS[i]) for i in range(3)]
    A = [q*Xi for Xi in X]
    F = {}
    for i in range(3):
        for j in range(i+1, 3):
            F[(i, j)] = sp.diff(A[j], XS[i]) - sp.diff(A[i], XS[j]) + A[i]*A[j] - A[j]*A[i]
    return X, F


def at_pole(M):
    """evaluate on the z axis at distance r (n = C_3)."""
    out = M.subs({x: 0, y: 0}).doit()
    return sp.simplify(out.subs(z, r).doit())


def run():
    half = sp.Rational(1, 2)
    X, F_half = field(half)
    _, F_one = field(1)
    # N1
    comm = {(i, j): X[i]*X[j] - X[j]*X[i] for i in range(3) for j in range(i+1, 3)}
    n1_half = all(at_pole(F_half[k] + comm[k]/4) == sp.zeros(2, 2) for k in comm)
    n1_one = all(at_pole(F_one[k]) == sp.zeros(2, 2) for k in comm)
    if not (n1_half and n1_one):
        raise ValueError('NT-3 on the field failed')
    # N2
    p = psi(r)
    dp = sp.diff(p, r)
    Fzx = at_pole(F_half[(0, 2)])
    Fzy = at_pole(F_half[(1, 2)])
    Fxy = at_pole(F_half[(0, 1)])
    s_ra = sp.simplify(sc(Fzx*Fzx))
    s_ra2 = sp.simplify(sc(Fzy*Fzy))
    s_aa = sp.simplify(sc(Fxy*Fxy))
    want_ra = -(dp*sp.sinh(p))**2/(4*r**2)
    want_aa = -sp.sinh(p)**4/(4*r**4)
    if sp.simplify(s_ra - want_ra) != 0 or sp.simplify(s_ra2 - want_ra) != 0 or sp.simplify(s_aa - want_aa) != 0:
        raise ValueError('invariant sizes of the native curvature failed')
    # N2b: this is RMG1's formula f = (1/4) sinh(l/2) dl ^ dphi with eigenvalue ratio e^l, l = 2 psi,
    #      written per unit coordinate area r dr dphi
    ell = 2*p
    rmg1 = sp.Rational(1, 4)*sp.sinh(ell/2)*sp.diff(ell, r)/r
    if sp.simplify(rmg1**2 - (-s_ra)) != 0:
        raise ValueError('native size is not the pull-back of the hyperbolic area form (RMG1)')
    # N3: ratio law
    rs = sp.Symbol('r_s', positive=True)
    tt = sp.sqrt(rs/r)                                   # tanh(eta) of the field of CV1
    deta = sp.diff(tt, r)/(1 - tt**2)                    # eta' from tanh(eta) = tt
    sinh2 = 2*tt/(1 - tt**2)
    ratio_gr = sp.simplify(r*2*deta/sinh2)
    if sp.simplify(ratio_gr + half) != 0:
        raise ValueError('the ratio law does not give the field of CV1')
    eta_f = sp.Function('eta')(r)
    ode = sp.dsolve(sp.Eq(r*sp.diff(eta_f, r), -sp.Rational(1, 4)*sp.sinh(2*eta_f)), eta_f) if False else None
    # converse: rho = -1/2 means d ln tanh(eta) / d ln r = -1/2
    converse = sp.simplify(r*sp.diff(sp.log(sp.tanh(eta_f)), r) - (r*sp.diff(2*eta_f, r)/sp.sinh(2*eta_f)))
    if converse != 0:
        raise ValueError('ratio identity failed')
    # N4: variational law
    w = sp.Function('w')(r)
    lagr = sp.Rational(1, 2)*sp.diff(w, r)**2 + sp.Rational(1, 4)*(w**2 - 1)**2/r**2     # -(r^2/1) x total invariant square, in w = cosh psi
    check = sp.simplify((2*(-s_ra) + (-s_aa))*r**2 - (sp.Rational(1, 2)*dp**2*sp.sinh(p)**2 + sp.Rational(1, 4)*sp.sinh(p)**4/r**2))
    if check != 0:
        raise ValueError('reduced curvature-square failed')
    from sympy.calculus.euler import euler_equations
    el = sp.simplify(euler_equations(lagr, [w], [r])[0].lhs)
    law = sp.simplify(el - (-(sp.diff(w, r, 2)) + w*(w**2 - 1)/r**2))
    if law != 0:
        raise ValueError('variational law is not w\'\' = w (w^2 - 1) / r^2')
    A, B = sp.symbols('A B')
    lin = lambda f: sp.simplify(sp.diff(f, r, 2) - 2*f/r**2)
    if lin(A/r + B*r**2) != 0:
        raise ValueError('linearised solutions failed')
    series = 1 + A/r + sp.Rational(3, 4)*A**2/r**2
    resid = sp.series(sp.simplify(sp.diff(series, r, 2) - series*(series**2 - 1)/r**2), r, sp.oo, 5).removeO()
    if sp.simplify(resid) != 0:
        raise ValueError('second-order coefficient 3/4 failed')
    residuals = {}
    for name, wgr in (('psi = 2 eta', (1 + rs/r)/(1 - rs/r)), ('psi = eta', 1/sp.sqrt(1 - rs/r))):
        res = sp.simplify(sp.diff(wgr, r, 2) - wgr*(wgr**2 - 1)/r**2)
        residuals[name] = str(sp.simplify(sp.series(res, r, sp.oo, 5).removeO()))
        if sp.simplify(res) == 0:
            raise ValueError('the field of N3 must not solve the variational law')
    # N5: kappa of the variational law in each reading
    eps = sp.Symbol('epsilon', positive=True)       # eps = r_s / r
    def kappa(reading):
        if reading == 'psi = 2 eta':                 # w = cosh 2 eta = (1 + m)/(1 - m)
            a = 2*rs                                  # first order: w - 1 = 2 m = 2 r_s / r
            wser = 1 + a/r + sp.Rational(3, 4)*a**2/r**2
            m = sp.series(((wser - 1)/(wser + 1)).subs(r, rs/eps), eps, 0, 3).removeO()
        else:                                         # w = cosh eta = 1/sqrt(1 - m)
            a = rs/2
            wser = 1 + a/r + sp.Rational(3, 4)*a**2/r**2
            m = sp.series((1 - 1/wser**2).subs(r, rs/eps), eps, 0, 3).removeO()
        return sp.nsimplify(sp.expand(m).coeff(eps, 2)), sp.expand(m).coeff(eps, 1)
    kap = {}
    for reading in ('psi = 2 eta', 'psi = eta'):
        k2, k1 = kappa(reading)
        if k1 != 1:
            raise ValueError('first order must be r_s / r')
        kap[reading] = k2
    mercury = 42.98
    advance = {name: float(mercury*(3 + 2*k)/3) for name, k in kap.items()}
    return dict(nt3_on_field=True, equals_rmg1_area_form=True, radial_angular=str(s_ra), angular_angular=str(s_aa),
                ratio_of_field=str(ratio_gr), variational_residual_of_ratio_field=residuals,
                second_order_coefficient='3/4', kappa={k: str(v) for k, v in kap.items()},
                mercury_arcsec_per_century=dict(ratio_law=mercury, **advance))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('NC1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, ':', v)
