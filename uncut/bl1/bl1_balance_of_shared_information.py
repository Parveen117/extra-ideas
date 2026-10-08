"""BL1: the balance of the shared information I - what it is exchanged for, and what one cycle adds.

Sources read (unchanged): Publications emk-ugd-algebra EMK-1 (det = seen + lost); extra-ideas uncut/ci1 (I), tt1 (cycle = boost);
physics/tl1 L2 (u_A N_A = u_B N_B), ln1 N5 (seen part = lost part for the unit block), rmg10 T6 via tt1 (tanh eta = t/rho).
No measurement.  sympy, exact.  Python 3.12.
"""
import json

import sympy as sp

I2 = sp.eye(2)
R = sp.Matrix([[0, -1], [1, 0]])
K = sp.Matrix([[1, 0], [0, -1]])
S = R*K


def parts(M):
    return ((M[0, 0] + M[1, 1])/2, (M[0, 0] - M[1, 1])/2, (M[0, 1] + M[1, 0])/2, (M[1, 0] - M[0, 1])/2)


def channels(M, across=False):
    """(seen, lost) for the cut along K (or across it)"""
    a, b, c, d = parts(M)
    return (a*a - c*c, d*d - b*b) if across else (a*a - b*b, d*d - c*c)


def zero(e):
    return sp.simplify(sp.expand(sp.sympify(e).rewrite(sp.exp))) == 0


def run():
    res = {}
    a, b, c, d = sp.symbols('a b c d', real=True)
    M = a*I2 + b*K + c*S + d*R

    # T1 seen + lost = det, in every frame; so exp(2 I_cut) - 1 = -lost/det : I rises exactly as the lost channel grows
    seen, lost = channels(M)
    if not zero(seen + lost - M.det()):
        raise ValueError('channels')
    if not zero(seen/M.det() - 1 + lost/M.det()):
        raise ValueError('balance')
    psi = sp.Symbol('psi', real=True)
    G = sp.cosh(psi)*I2 + sp.sinh(psi)*S                   # unit block across the cut
    sg, lg = channels(G)
    if not (zero(sg - sp.cosh(psi)**2) and zero(lg + sp.sinh(psi)**2) and zero(sg + lg - 1)):
        raise ValueError('unit block channels')
    res['balance'] = 'seen + lost = det ; for a unit block seen = cosh^2, lost = -sinh^2 : what is seen more is what is lost more'

    # T2 k cycles of an open diagram from rest: exp(I_0) = cosh(2 k eta); steps grow toward 2 eta = Log((rho + t)/(rho - t))
    wv = sp.Rational(3, 4)
    kap = sp.Rational(5, 4)*K + wv*R
    cyc = (R*kap)**2
    vals, Gk = [], I2
    for k in range(1, 6):
        Gk = cyc*Gk
        vals.append(parts(Gk)[0])                          # = cosh(2 k eta), exact
    eta2 = sp.log(sp.Rational(17, 8) + sp.Rational(15, 8)) # exp(2 eta) = cosh + sinh = 4
    if sp.simplify(eta2 - sp.log(4)) != 0:
        raise ValueError('rapidity of one cycle')
    for k, v in enumerate(vals, 1):
        if v != (sp.Integer(4)**k + sp.Rational(1, 4)**k)/2:
            raise ValueError('cosh(2 k eta)')
    ratios = [vals[i + 1]/vals[i] for i in range(len(vals) - 1)]
    if not (vals[0] == 1 + 2*wv*wv and all(ratios[i] < ratios[i + 1] < 4 for i in range(len(ratios) - 1)) and ratios[0] > vals[0]):
        raise ValueError('steps')
    yv = sp.Symbol('y', positive=True)
    if sp.simplify(((1 + sp.tanh(yv))/(1 - sp.tanh(yv)) - sp.exp(2*yv)).rewrite(sp.exp)) != 0:
        raise ValueError('limit step')
    res['cycles'] = {'first': 'exp(I_0) = 1 + 2 w^2', 'exp_I0_after_k_cycles_w_3_4': [str(v) for v in vals],
                     'step_factors': [str(r) for r in ratios], 'limit_step': 'exp(Delta I) -> exp(2 eta) = (rho + t)/(rho - t)  (here 4)'}

    # T3 the reverse cycle undoes it: the order is the order of repeating one sense
    if (cyc.inv()*cyc) != I2 or parts(cyc.inv())[0] != parts(cyc)[0]:
        raise ValueError('reverse')
    res['reverse'] = 'the reverse cycle returns to rest: nothing here forbids it'

    # T4 between two places (TL1-L2): u N is kept, N = exp(-I) : Log(u) - I is the same at both places
    NA, NB, uA = sp.Rational(8, 17), sp.Rational(3, 5), sp.Rational(2)
    uB = uA*NA/NB
    if sp.simplify((sp.log(uA) + sp.log(NA)) - (sp.log(uB) + sp.log(NB))) != 0 or not (uB < uA and NB > NA):
        raise ValueError('units and clock')
    res['exchange'] = 'Log(local unit) - I is kept between places that exchange: where I is larger the local unit is larger by exp(Delta I)'
    return res


if __name__ == '__main__':
    out = run()
    with open(__file__.replace('bl1_balance_of_shared_information.py', 'BL1_RESULT.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
