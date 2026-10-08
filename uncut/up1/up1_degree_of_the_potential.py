"""UP1: the potential before any cut, its degree, and what the degree does to the diagram and to the tower.

No measurement, no physical constant.  Sources read (unchanged):
  Publications research branch  coherence-first-thermodynamics CF-1 (the centre is a boundary of the cut chart),
                                CF-3 (W = kappa r^nu + eps s^2/2, 1 < nu < 2), CF-4 (density ~ r^(nu-2))
                                uncut-cut-measurement MANUSCRIPT section 1 (the uncut is not a set of states)
                                cut-first-equivalence (the classical laws are the flat, memoryless sector)
  Recognition-Kernel-Framework  theorum/42 Theorem 3.1 (typed tower L_{n+1} = nabla L_n)
  response-geometry series      RMG2 T2, T3 (tower T_lambda(H) = H + lambda H^2; kappa, fixed sets)
  extra-ideas physics           LT1 (one generation at infinite lambda is squaring)
Symbolic algebra with sympy.  Python 3.12.
"""
import itertools
import json

import sympy as sp

k, s, lam, tt = sp.symbols('k s lambda t', positive=True)


def generic(n):
    """a family of potentials of degree k in n positive variables: free coefficients, free exponents, and a radial term"""
    xs = sp.symbols(f'x1:{n + 1}', positive=True)
    c1, c2, c3 = sp.symbols('c1 c2 c3', real=True)
    a, b = sp.symbols('a b', real=True)
    rad = sp.sqrt(sum(x*x for x in xs))**k
    if n == 2:
        return xs, c1*xs[0]**(k - a)*xs[1]**a + c2*xs[0]**(k - b)*xs[1]**b + c3*rad
    return xs, c1*xs[0]**(k - a - b)*xs[1]**a*xs[2]**b + c2*xs[1]**(k - a)*xs[2]**a + c3*rad


def euler(expr, xs):
    return sum(x*sp.diff(expr, x) for x in xs)


def zero(e):
    return sp.simplify(e) == 0


def cut(phi, xs, subset):
    """the reading of phi after cutting the pairs in subset: phi - sum x_i y_i"""
    return phi - sum(xs[i]*sp.diff(phi, xs[i]) for i in subset)


def partial_cut(phi, xs, depth):
    return phi - sum(d*x*sp.diff(phi, x) for d, x in zip(depth, xs))


def tower(H, l):
    return H + l*H*H


def run():
    res = {}
    for n in (2, 3):
        xs, phi = generic(n)
        # T1 the dilation generator reads the degree; each layer of the tower drops it by one
        if not zero(euler(phi, xs) - k*phi):
            raise ValueError('degree')
        for i in range(n):
            yi = sp.diff(phi, xs[i])
            if not zero(euler(yi, xs) - (k - 1)*yi):
                raise ValueError('first layer')
            for j in range(n):
                hij = sp.diff(yi, xs[j])
                if not zero(euler(hij, xs) - (k - 2)*hij):
                    raise ValueError('second layer')
        if n == 2:
            third = sp.diff(phi, xs[0], 2, xs[1])
            if not zero(euler(third, xs) - (k - 3)*third):
                raise ValueError('third layer')
        # T2 the corners: every face closes, the mean of all corners is (1 - k/2) phi, the far corner (1 - k) phi
        corners = {A: cut(phi, xs, A) for r in range(n + 1) for A in itertools.combinations(range(n), r)}
        for i, j in itertools.combinations(range(n), 2):
            rest = [q for q in range(n) if q not in (i, j)]
            for r in range(len(rest) + 1):
                for B in itertools.combinations(rest, r):
                    a = tuple(sorted(B))
                    f = corners[a] + corners[tuple(sorted(B + (i, j)))] - corners[tuple(sorted(B + (i,)))] - corners[tuple(sorted(B + (j,)))]
                    if not zero(f):
                        raise ValueError('face relation')
        mean = sum(corners.values())/len(corners)
        if not zero(mean - (1 - k/2)*phi) or not zero(corners[tuple(range(n))] - (1 - k)*phi):
            raise ValueError('centre or far corner')
        # T3 along the diagonal the reading is (1 - t k) phi: silent at depth 1/k
        diag = partial_cut(phi, xs, [tt]*n)
        if not zero(diag - (1 - tt*k)*phi) or not zero(diag.subs(tt, 1/k)):
            raise ValueError('diagonal')
    res['corners'] = 'mean of all corners = (1 - k/2) phi ; far corner = (1 - k) phi ; diagonal depth t reads (1 - t k) phi'

    # T4 a potential of mixed degree: depth 1/k removes exactly its degree-k part (CF-3's family)
    sv, vv, kap, eps, nu, T0, P0 = sp.symbols('s_ v_ kappa epsilon nu T_0 P_0', positive=True)
    r = sp.sqrt(sv**2 + vv**2)
    parts = {1: T0*sv - P0*vv, nu: kap*r**nu, 2: eps*sv**2/2}
    total = sum(parts.values())
    for d, piece in parts.items():
        got = partial_cut(total, (sv, vv), [1/sp.sympify(d)]*2)
        want = sum((1 - dd/sp.sympify(d))*pp for dd, pp in parts.items())
        marks = {1: (T0, P0), nu: (kap,), 2: (eps,)}[d]
        if not zero(got - want) or any(sp.simplify(got).has(q) for q in marks):
            raise ValueError('degree filter')
    centre = sp.simplify(partial_cut(total, (sv, vv), [sp.Rational(1, 2)]*2))
    if not zero(centre - ((T0*sv - P0*vv)/2 + (1 - nu/2)*kap*r**nu)) or centre.has(eps):
        raise ValueError('centre of CF-3')
    far = sp.simplify(partial_cut(total, (sv, vv), [1, 1]))
    if not zero(far - ((1 - nu)*kap*r**nu - eps*sv**2/2)) or far.has(T0) or far.has(P0):
        raise ValueError('far corner of CF-3')
    res['cf3'] = {'centre': 'half the linear part + (1 - nu/2) kappa r^nu, no epsilon',
                  'far_corner': '(1 - nu) kappa r^nu - eps s^2/2, no T_0, P_0',
                  'window_1_lt_nu_lt_2': 'formation part positive at the centre, negative at the far corner'}

    # T5 the tower and the place are one parameter: T_lambda at s x = s^(k-2) T_{lambda s^(k-2)} at x
    xs, phi = generic(2)
    H = sp.hessian(phi, xs)
    Hs = sp.simplify(sp.hessian(phi.subs({xs[0]: s*xs[0], xs[1]: s*xs[1]}, simultaneous=True), xs)/s**2)
    if sp.simplify(Hs - s**(k - 2)*H) != sp.zeros(2, 2):
        raise ValueError('response element degree')
    lhs = tower(s**(k - 2)*H, lam)
    rhs = s**(k - 2)*tower(H, lam*s**(k - 2))
    if sp.simplify(lhs - rhs) != sp.zeros(2, 2):
        raise ValueError('tower covariance')
    res['tower'] = 'lambda and place enter only through lambda s^(k-2)'

    # T6 RMG2's kappa at the two ends of lambda m: identity far away, squaring at the centre (k < 2)
    m, rho = sp.symbols('m rho', positive=True)
    kappa_t = (1 + 2*lam*m)/(1 + lam*m*(1 + rho**2))
    if sp.limit(kappa_t, m, 0) != 1 or not zero(sp.limit(kappa_t, m, sp.oo) - 2/(1 + rho**2)):
        raise ValueError('ends of the tower')
    if not zero(sp.limit(kappa_t, m, sp.oo)*rho - 2*rho/(1 + rho**2)):
        raise ValueError('squaring map')
    res['ends'] = {'lambda m -> 0': 'identity', 'lambda m -> infinity': 'rho -> 2 rho/(1 + rho^2), no lambda left'}

    # T7 which layer is at rest under dilation: layer n has degree k - n
    res['layer_at_rest'] = 'n = k ; none when k is not a whole number'
    return res


if __name__ == '__main__':
    out = run()
    with open('UP1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
