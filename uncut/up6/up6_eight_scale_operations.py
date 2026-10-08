"""UP6: the eight scale operations of the diagram, and where they fail to commute.

The owner's statement: the two forms of the mechanical reading - P (dV/dP)_T and V (dP/dV)_T - are not a slip to
be settled; their mismatch is the reason for the connection of the lambda tower, and their non-commutation may
be what is lost at lambda = 0.

A scale operation is  E[x | y] f = x (df/dx) at fixed y.  Each corner of the diagram has two.
No measurement, no physical constant.  Sources read (unchanged): the seed document; uncut/up1, up3, up4, up5.
Chart (a, b) = (S, V).  Symbolic algebra with sympy.  Python 3.12.
"""
import json

import sympy as sp

a, b = sp.symbols('a b', positive=True)
S, V = a, b


def br(f, g):
    return sp.diff(f, a)*sp.diff(g, b) - sp.diff(f, b)*sp.diff(g, a)


def E(x, y):
    """the scale operation E[x | y] as a function acting on readings"""
    return lambda f: x*br(f, y)/br(x, y)


def comm(D1, D2, f):
    return D1(D2(f)) - D2(D1(f))


def zero(e):
    return sp.simplify(e) == 0


def setup(U):
    T, P = sp.diff(U, S), -sp.diff(U, V)
    ops = {'U_S': E(S, V), 'U_V': E(V, S),          # corner U(S, V)
           'H_S': E(S, P), 'H_P': E(P, S),          # corner H(S, P)
           'F_T': E(T, V), 'F_V': E(V, T),          # corner F(T, V)
           'G_T': E(T, P), 'G_P': E(P, T)}          # corner G(T, P)
    return T, P, ops


def run():
    res = {}
    Ug = sp.Function('U')(a, b)
    T, P, ops = setup(Ug)
    h = sp.Function('h')(a, b)

    # T1 the eight readings: each operation applied to the partner of its own variable; four reciprocal pairs
    lam_v, lam_p = -ops['U_S'](T), -ops['H_S'](T)                 # -S (dT/dS) at fixed V, P
    C_v, C_p = ops['F_T'](S), ops['G_T'](S)
    lam_s, lam_tB = ops['U_V'](P), ops['F_V'](P)                  # V (dP/dV) at fixed S, T
    z_s, lam_tA = -ops['H_P'](V), -ops['G_P'](V)                  # -P (dV/dP) at fixed S, T
    for x, y, prod in ((lam_v, C_v, -S*T), (lam_p, C_p, -S*T), (lam_s, z_s, -P*V), (lam_tB, lam_tA, -P*V)):
        if not zero(x*y - prod):
            raise ValueError('reciprocal pair')
    if not zero(lam_p/lam_v - lam_tB/lam_s):
        raise ValueError('shared ratio')
    res['readings'] = 'eight operations, eight readings, four reciprocal pairs (products -ST, -ST, -PV, -PV); form (A) and form (B) are one pair'

    # T2 inside a corner the two operations commute
    for c in ('U', 'H', 'F', 'G'):
        k1, k2 = [k for k in ops if k.startswith(c + '_')]
        if not zero(comm(ops[k1], ops[k2], h)):
            raise ValueError('corner ' + c)
    res['inside_a_corner'] = 'commute'

    # T3 between corners: the defect is a scale-derivative of a pure ratio
    m = ops['F_V'](S)/S                                           # (d ln S / d ln V) at fixed T
    eT = ops['F_V'](P)/P                                          # (d ln P / d ln V) at fixed T
    if not zero(ops['F_V'](h) - ops['U_V'](h) - m*ops['U_S'](h)):
        raise ValueError('F_V in the U corner')
    if not zero(comm(ops['U_S'], ops['F_V'], h) - ops['U_S'](m)*ops['U_S'](h)):
        raise ValueError('[U_S, F_V]')
    if not zero(comm(ops['U_V'], ops['F_V'], h) - ops['U_V'](m)*ops['U_S'](h)):
        raise ValueError('[U_V, F_V]')
    if not zero(ops['G_P'](h) - ops['F_V'](h)/eT):
        raise ValueError('form A against form B')
    if not zero(comm(ops['F_V'], ops['G_P'], h) + ops['F_V'](sp.log(eT))*ops['G_P'](h)):
        raise ValueError('[F_V, G_P]')
    res['between_corners'] = {'[U_S, F_V]': 'U_S(m) U_S', '[U_V, F_V]': 'U_V(m) U_S',
                              '[F_V, G_P]': '- F_V(ln e_T) G_P',
                              'm': '(d ln S / d ln V) at fixed T', 'e_T': '(d ln P / d ln V) at fixed T'}

    # T4 the cost of UP4, inside ordinary thermodynamics: one cut from the U corner, one from the F corner
    DS, DV = ops['U_S'], ops['F_V']
    first = DS(Ug)
    if not zero(first - S*T):
        raise ValueError('first layer along U_S')
    B1, B2 = DV(DS(Ug)), DS(DV(Ug))
    if not zero((B2 - B1) - ops['U_S'](m)*S*T):
        raise ValueError('mixed responses')
    res['cost'] = 'for the pair (U_S, F_V): B_2 - B_1 = U_S(m) * S T'

    # T5 witnesses
    out = {}
    for name, U in (('power', a**3/b), ('separable_powers', a**3 + 2*sp.sqrt(b)),
                    ('not_a_power', a**2/2 + a*b/3 + b**2/2 + a**2*b/5 + b**3/7)):
        Tn, Pn, on = setup(U)
        mn = on['F_V'](S)/S
        en = on['F_V'](Pn)/Pn
        d1 = sp.simplify(on['U_S'](mn))
        d2 = sp.simplify(on['U_V'](mn))
        d3 = sp.simplify(on['F_V'](sp.log(en)))
        at = {a: sp.Rational(1, 2), b: sp.Rational(1, 2)}
        out[name] = [str(sp.nsimplify(sp.simplify(q.subs(at)))) for q in (d1, d2, d3)]
    if out['power'] != ['0', '0', '0'] or out['separable_powers'] != ['0', '0', '0'] or '0' in out['not_a_power']:
        raise ValueError('witnesses')
    res['witnesses_defects_at_(1/2,1/2)'] = out
    return res


if __name__ == '__main__':
    r = run()
    with open('UP6_RESULT.json', 'w') as fh:
        json.dump(r, fh, indent=1)
    print(json.dumps(r, indent=1))
