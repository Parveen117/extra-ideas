"""UP2: the two-point potential - one function on both halves of the diagram; the cut surface is its rest set,
and its layers are the lambda level and every higher level.

No measurement, no physical constant.  Sources read (unchanged):
  Recognition-Kernel-Framework  theorum/42 Theorem 3.1 (typed tower L_{n+1} = nabla L_n), section 7 (higher layers repair)
  Publications research branch  coherence-first-thermodynamics CF-3 (7) (C_V, C_P, K_S, K_T from the Hessian)
                                thermo-compass-foundations RESULTS_INDEX (C_P/C_V = K_S/K_T)
  response-geometry series      RMG2 T1 (number sectors; boundary rho = 1 is the dual sector), RMG9 T5 (Schur ratio)
  extra-ideas                   uncut/up1 (degree k; silent depth 1/k), physics/lt1 (single null reading)
Symbolic algebra with sympy.  Python 3.12.
"""
import json

import sympy as sp

k = sp.Symbol('k', positive=True)


def grad(phi, xs):
    return [sp.diff(phi, q) for q in xs]


def two_point(phi, xs, base):
    """G(x; x0) = phi(x) - phi(x0) - grad phi(x0).(x - x0)"""
    sub = dict(zip(xs, base))
    g0 = [q.subs(sub, simultaneous=True) for q in grad(phi, xs)]
    return phi - phi.subs(sub, simultaneous=True) - sum(g*(q - b) for g, q, b in zip(g0, xs, base))


def zero(e):
    return sp.simplify(e) == 0


def run():
    res = {}
    x, x0, y, s = sp.symbols('x x0 y s', positive=True)
    ks = k + 1                                    # work with degree K = k + 1 > 1 so that the conjugate degree exists
    kc = ks/(ks - 1)

    # T1 pure degree, one pair: the two-point potential in the readings (x, y)
    phi = x**ks/ks
    dual = y**kc/kc
    G = phi + dual - x*y
    surface = {y: x**(ks - 1)}
    if not zero(two_point(phi, (x,), (x0,)) - G.subs(y, x0**(ks - 1))):
        raise ValueError('two forms of the two-point potential differ')
    if not (zero(G.subs(surface)) and zero(sp.diff(G, x).subs(surface)) and zero(sp.diff(G, y).subs(surface))):
        raise ValueError('cut surface is not the rest set')
    Hd = sp.hessian(G, (x, y)).subs(surface)
    h = (ks - 1)*x**(ks - 2)
    if sp.simplify(Hd - sp.Matrix([[h, -1], [-1, 1/h]])) != sp.zeros(2, 2) or not zero(Hd.det()):
        raise ValueError('doubled response')
    m = Hd.trace()/2
    u, v = (Hd[0, 0] - Hd[1, 1])/(2*m), Hd[0, 1]/m
    if not zero(u**2 + v**2 - 1):
        raise ValueError('not on the boundary rho = 1')
    if not zero(G.subs({x: s*x, y: s**(ks - 1)*y}, simultaneous=True) - s**ks*G):
        raise ValueError('weights')
    if not zero(1/ks + 1/kc - 1):
        raise ValueError('conjugate depths')
    res['one_pair'] = {'rest_set': 'y = x^(K-1): value 0, gradient 0',
                       'doubled_response': '[[h, -1], [-1, 1/h]], determinant 0, rho = 1',
                       'weights': 'x: 1, y: K-1, G: K',
                       'depths': '1/K + 1/K* = 1'}

    # T2 the layers of G at a base point are the tower layers from the second on; the first two are absent
    c = sp.symbols('c0:7', real=True)
    xi = sp.Symbol('xi', real=True)
    poly = sum(c[i]*x**i for i in range(7))
    Gp = two_point(poly, (x,), (x0,)).subs(x, x0 + xi)
    layers = sum(sp.diff(poly, x, n).subs(x, x0)*xi**n/sp.factorial(n) for n in range(2, 7))
    if not zero(sp.expand(Gp - layers)):
        raise ValueError('layers')
    a0, a1 = sp.symbols('a0 a1', real=True)
    if not zero(two_point(poly + a0 + a1*x, (x,), (x0,)) - two_point(poly, (x,), (x0,))):
        raise ValueError('affine part should be absent')
    x1 = sp.Symbol('x1', positive=True)
    gp = lambda p: sp.diff(poly, x).subs(x, p)
    three = two_point(poly, (x,), (x1,)) + two_point(poly, (x,), (x0,)).subs(x, x1) + (gp(x1) - gp(x0))*(x - x1)
    if not zero(sp.expand(two_point(poly, (x,), (x0,)) - three)):
        raise ValueError('change of base point')
    res['layers'] = 'G(x0 + xi; x0) = sum over n >= 2 of L_n(x0) xi^n / n!  (value and first layer absent)'

    # T3 two pairs, coupled, pure degree: phi = r^K / K
    p, q, p0, q0 = sp.symbols('p q p0 q0', positive=True)
    yp, yq = sp.symbols('y_p y_q', positive=True)
    r = sp.sqrt(p*p + q*q)
    phi2 = r**ks/ks
    dual2 = sp.sqrt(yp*yp + yq*yq)**kc/kc
    G2 = phi2 + dual2 - p*yp - q*yq
    surf2 = {yp: r**(ks - 2)*p, yq: r**(ks - 2)*q}
    for expr in (G2, sp.diff(G2, p), sp.diff(G2, q), sp.diff(G2, yp), sp.diff(G2, yq)):
        if not zero(sp.powsimp(sp.simplify(expr.subs(surf2, simultaneous=True)), force=True)):
            raise ValueError('two pairs: rest set')
    num = {p: sp.Rational(3), q: sp.Rational(4), k: sp.Rational(2)}        # K = 3, K* = 3/2, r = 5
    H4 = sp.hessian(G2, (p, q, yp, yq)).subs(surf2, simultaneous=True).subs(num)
    H4 = H4.applyfunc(sp.nsimplify)
    if H4.rank() != 2:
        raise ValueError('doubled response should have rank 2 of 4')
    res['two_pairs'] = 'rest set of dimension 2 in 4; doubled response of rank 2'

    # T4 cutting one pair of G by its rest condition: the response becomes the Schur complement
    A, B, C = sp.symbols('A B C', positive=True)
    quad = (A*p*p + 2*B*p*q + C*q*q)/2
    Gq = two_point(quad, (p, q), (p0, q0))
    pstar = sp.solve(sp.diff(Gq, p), p)[0]
    reduced = sp.simplify(Gq.subs(p, pstar))
    if not zero(reduced - (A*C - B*B)/(2*A)*(q - q0)**2):
        raise ValueError('Schur complement')
    T_, V_ = sp.symbols('T V', positive=True)
    det = A*C - B*B
    CV, CP, KS, KT = T_/A, T_*C/det, V_*C, V_*det/A                      # CF-3 (7)
    if not zero(CP/CV - KS/KT) or not zero((CP/CV)*(KT/KS) - 1):
        raise ValueError('closure')
    res['one_pair_cut'] = 'response C -> (AC - B^2)/A ; C_P/C_V = K_S/K_T is this ratio read on the two pairs'
    return res


def dual_check(K, xv, x0v):
    """exact check of G_phi(x; x0) = G_dual(y0; y) for phi = x^K / K at rational points"""
    xv, x0v, K = sp.nsimplify(xv), sp.nsimplify(x0v), sp.nsimplify(K)
    Kc = K/(K - 1)
    f = lambda t: t**K/K
    fd = lambda t: t**Kc/Kc
    yv, y0v = xv**(K - 1), x0v**(K - 1)
    left = f(xv) - f(x0v) - x0v**(K - 1)*(xv - x0v)
    right = fd(y0v) - fd(yv) - yv**(Kc - 1)*(y0v - yv)
    return sp.simplify(left - right) == 0


if __name__ == '__main__':
    out = run()
    out['duality_at_rational_points'] = all(dual_check(*t) for t in ((3, 4, 9), (2, 5, 7), (sp.Rational(3, 2), 4, 9)))
    with open('UP2_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
