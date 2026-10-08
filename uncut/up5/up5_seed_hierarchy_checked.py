"""UP5: the owner's seed document (thermodynamic hierarchy L_0 -> L_1 -> L_2, higher-order lambda thermodynamics)
checked against UP3 and UP4.

Seed definitions (layer 1):  lambda_p = -ST/C_p,  lambda_v = -ST/C_v,  lambda_s = V (dP/dV)_S,  z_s = -PV/lambda_s,
                             lambda_t = -P (dV/dP)_T,  z_t = -PV/lambda_t.
No measurement, no physical constant.  Sources read (unchanged): the seed document; uncut/up3, up4; RMG2 T5.
Symbolic algebra with sympy.  Python 3.12.
"""
import json

import sympy as sp

a, b = sp.symbols('a b', real=True)


def br(f, g):
    return sp.diff(f, a)*sp.diff(g, b) - sp.diff(f, b)*sp.diff(g, a)


def at_fixed(y, x, c):
    """(dy/dx) at fixed c"""
    return br(y, c)/br(x, c)


def zero(e):
    return sp.simplify(e) == 0


def seed_layer(T, V, S, P):
    lp = -S*T/(T*at_fixed(S, T, P))
    lv = -S*T/(T*at_fixed(S, T, V))
    ls = V*at_fixed(P, V, S)
    lt = -P*at_fixed(V, P, T)
    return dict(lp=lp, lv=lv, ls=ls, lt=lt, zs=-P*V/ls, zt=-P*V/lt)


def run():
    res = {}
    T, V, S, P = [sp.Function(n)(a, b) for n in 'TVSP']

    # T1 the four Maxwell-type relations of any layer are one equation between the two area forms
    t, v, s, p = T, V, S, P
    rel = {'F': at_fixed(s, v, t) - at_fixed(p, t, v),
           'G': at_fixed(s, p, t) + at_fixed(v, t, p),
           'H': at_fixed(t, p, s) - at_fixed(v, s, p),
           'U': at_fixed(t, v, s) + at_fixed(p, s, v)}
    gap = br(t, s) - br(p, v)
    den = {'F': -br(v, t), 'G': -br(p, t), 'H': br(p, s), 'U': br(v, s)}
    for name in rel:
        if not zero(rel[name] - gap/den[name]):
            raise ValueError('relation ' + name)
    res['four_relations'] = 'each equals ({t,s} - {p,v}) over one bracket: one equation, not four'

    # T2 the seed's layer 1 with no potential assumed: what its own definitions imply
    L = seed_layer(T, V, S, P)
    chi = br(T, P)*br(S, V)/(br(S, P)*br(T, V))
    if not zero(L['lp']/L['lv'] - chi):
        raise ValueError('lambda_p / lambda_v')
    if not zero(L['zt']/L['ls'] - chi):
        raise ValueError('z_t / lambda_s')
    if not zero(L['lp']*L['ls'] - L['lv']*L['zt']):
        raise ValueError('closed identity with z_t')
    if not zero(L['lp']*L['ls']*L['lt'] + P*V*L['lv']):
        raise ValueError('identity with lambda_t as written')
    if zero(L['lp']*L['ls'] - L['lv']*L['lt']):
        raise ValueError('lambda_t as written should not close the layer')
    res['layer_1'] = {'holds': 'lambda_p/lambda_v = z_t/lambda_s ; lambda_p lambda_s = lambda_v z_t ; lambda_p lambda_s lambda_t = -P V lambda_v',
                      'does_not_hold': 'lambda_p lambda_s = lambda_v lambda_t (with lambda_t = -P (dV/dP)_T as written)'}

    # T3 the general defect of two cut operations: frame part and weight part; it pairs with (value, first layer)
    f, g, pw, qw, h = [sp.Function(n)(a, b) for n in ('f', 'g', 'p_w', 'q_w', 'h')]
    DS = lambda x: sp.diff(x, a) + pw*x
    DV = lambda x: f*sp.diff(x, b) + g*sp.diff(x, a) + qw*x
    beta = sp.diff(f, a)/f
    alpha = sp.diff(g, a) - beta*g
    Fw = sp.diff(qw, a) - f*sp.diff(pw, b) - g*sp.diff(pw, a) - alpha*pw - beta*qw
    if not zero(DS(DV(h)) - DV(DS(h)) - alpha*DS(h) - beta*DV(h) - Fw*h):
        raise ValueError('general defect')
    U = sp.Function('U')(a, b)
    Tu, Pu = DS(U), -DV(U)
    B1, B2 = DV(Tu), -DS(Pu)
    if not zero((B2 - B1) - (alpha*Tu - beta*Pu + Fw*U)):
        raise ValueError('mixed responses')
    res['defect'] = 'B_2 - B_1 = alpha T - beta P + F U : the defect pairs with the value and the first layer only'

    # T4 "the state returns but the responses do not": false with commuting cuts, true with a defect
    lam = sp.Function('lam')(a, b)
    e = sp.Symbol('epsilon')
    closed = lam.subs({a: a + e - e, b: b + e - e}, simultaneous=True) - lam
    if not zero(closed):
        raise ValueError('commuting loop should return every reading')
    Uw = a**2 + b**2 - a*b/2 + 2*b
    a1 = a + e
    b2 = b + e*sp.exp(a1)
    a3 = a1 - e
    b4 = b2 - e*sp.exp(a3)
    back_a = sp.simplify(a3 - a)
    drift = sp.series(Uw.subs({a: a3, b: b4}, simultaneous=True) - Uw, e, 0, 3).removeO()
    if back_a != 0 or zero(drift):
        raise ValueError('defect loop')
    res['return'] = {'commuting': 'every reading returns', 'with_defect': 'the first coordinate returns, the potential does not (UP4-T5)'}
    return res


if __name__ == '__main__':
    out = run()
    with open('UP5_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
