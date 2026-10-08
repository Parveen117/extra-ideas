"""UP3: the sequence of diagrams.  Level n+1 is built from level n by the owner's rule

    P_{n+1} = S_n T_n / C_{P,n} ,   V_{n+1} = S_n T_n / C_{V,n} ,   C_{X,n} = T_n (dS_n/dT_n) at fixed X_n ,
    T_{n+1} = P_n V_n / W_{T,n} ,   S_{n+1} = P_n V_n / W_{S,n} ,   W_{X,n} = P_n (dV_n/dP_n) at fixed X_n ,

and the question is what holds at every level with no potential assumed.

No measurement, no physical constant.  Sources read (unchanged):
  Publications research branch  thermo-compass-foundations RESULTS_INDEX (Jacobian brackets; C_P/C_V = K_S/K_T)
                                coherence-first-thermodynamics CF-3 (7)
  Recognition-Kernel-Framework  theorum/42 (lambda-Jacobian tower; cut/transpose grading of the Jacobian)
  response-geometry series      RMG2 T5 (corner layer lambda_t = V Delta/a, lambda_s = -Vc, lambda_p = -S Delta/c,
                                lambda_v = -Sa; rank-one identity; closure only shown for separable U), RMG6 (memory weight r^2)
  extra-ideas uncut/up1, up2
A state is a point of a surface with coordinates (a, b); {f, g} = f_a g_b - f_b g_a.  Symbolic algebra with sympy.  Python 3.12.
"""
import itertools
import json

import sympy as sp

a, b = sp.symbols('a b', real=True)


def br(f, g):
    return sp.diff(f, a)*sp.diff(g, b) - sp.diff(f, b)*sp.diff(g, a)


def step(T, S, P, V):
    """one application of the rule, written with brackets: (dT/dS) at fixed X = {T, X}/{S, X}"""
    Pn = S*br(T, P)/br(S, P)
    Vn = S*br(T, V)/br(S, V)
    Tn = V*br(P, T)/br(V, T)
    Sn = V*br(P, S)/br(V, S)
    return Tn, Sn, Pn, Vn


def cross_ratio(T, S, P, V):
    return br(T, P)*br(S, V)/(br(S, P)*br(T, V))


def zero(e):
    return sp.simplify(e) == 0


def from_potential(U, S, V):
    """level 1 of a potential U(S, V) on the surface (a, b) = (S, V): T = U_S, P = -U_V"""
    return sp.diff(U, S), S, -sp.diff(U, V), V


def run():
    res = {}
    T, S, P, V = [sp.Function(n)(a, b) for n in 'TSPV']

    # T1 any four readings, no potential: both axes of the next level have the same ratio of ends
    T2, S2, P2, V2 = step(T, S, P, V)
    chi = cross_ratio(T, S, P, V)
    if not (zero(P2/V2 - chi) and zero(T2/S2 - chi)):
        raise ValueError('ratio of ends')
    if not zero(P2*S2 - V2*T2):
        raise ValueError('cross product')
    if not zero(cross_ratio(P, V, T, S) - chi):
        raise ValueError('exchange of the two axes')
    res['every_level'] = 'P/V = T/S = chi of the level below; P S = V T'

    # T2 the three-term bracket identity, and 1 - chi
    pl = br(T, S)*br(P, V) - br(T, P)*br(S, V) + br(T, V)*br(S, P)
    if not zero(sp.expand(pl)):
        raise ValueError('three-term identity')
    if not zero(1 - chi + br(T, S)*br(P, V)/(br(S, P)*br(T, V))):
        raise ValueError('1 - chi')
    res['one_minus_chi'] = '- {T,S}{P,V} / ({S,P}{T,V})'

    # T3 with a potential (level 1): {T,S} = {P,V}, so 1 - chi is a square over the two diagonal responses
    Sv, Vv = a, b
    U = sp.Function('U')(a, b)
    T1, S1, P1, V1 = from_potential(U, Sv, Vv)
    A, B, C = sp.diff(U, a, 2), sp.diff(U, a, b), sp.diff(U, b, 2)
    if not zero(br(T1, S1) - br(P1, V1)):
        raise ValueError('area forms of the two axes')
    chi1 = cross_ratio(T1, S1, P1, V1)
    if not zero(chi1 - (A*C - B*B)/(A*C)) or not zero(1 - chi1 - B*B/(A*C)):
        raise ValueError('level-1 ratio')
    L2 = step(T1, S1, P1, V1)
    want = (-Vv*(A*C - B*B)/A, -Vv*C, Sv*(A*C - B*B)/C, Sv*A)          # RMG2-T5's corner layer (signs of P = -U_V)
    if not all(zero(x - w) for x, w in zip(L2, want)):
        raise ValueError('level 2 is not the corner layer')
    res['level_1'] = '1 - chi = B^2/(AC) >= 0 for a positive response: chi <= 1'

    # T4 every level from the second on is (chi s, s, chi v, v); its area forms and its own ratio
    ch, s, v = [sp.Function(n)(a, b) for n in ('chi', 's', 'v')]
    Tq, Sq, Pq, Vq = ch*s, s, ch*v, v
    if not (zero(br(Tq, Sq) - s*br(ch, s)) and zero(br(Pq, Vq) - v*br(ch, v))):
        raise ValueError('area forms')
    eps = sp.Symbol('epsilon')
    if not zero(br(Tq, Sq) - eps*br(Pq, Vq) - br(ch, (s*s - eps*v*v)/2)):
        raise ValueError('potential condition')
    Aq = br(s, v)
    pp, qq = s*br(ch, v)/(ch*Aq), v*br(ch, s)/(ch*Aq)
    chi_next = cross_ratio(Tq, Sq, Pq, Vq)
    if not zero(1 - chi_next + pp*qq/((1 - qq)*(1 + pp))):
        raise ValueError('next ratio')
    res['level_n'] = {'area_forms': '{T,S} = s {chi, s} ; {P,V} = v {chi, v}',
                      'has_a_potential_iff': '{chi, s^2 - eps v^2} = 0 with eps = +1 or -1',
                      'next_ratio': '1 - chi\' = - p q / ((1 - q)(1 + p)),  p = s{chi,v}/(chi{s,v}),  q = v{chi,s}/(chi{s,v})'}

    # T5 a potential of pure power form: chi constant, then chi = 1 at every later level
    Um = a**3/b
    lv = from_potential(Um, a, b)
    c1 = sp.simplify(cross_ratio(*lv))
    lv2 = [sp.simplify(q) for q in step(*lv)]
    c2 = sp.simplify(cross_ratio(*lv2))
    lv3 = [sp.simplify(q) for q in step(*lv2)]
    c3 = sp.simplify(cross_ratio(*lv3))
    if not (c1 == sp.Rational(1, 4) and c2 == 1 and c3 == 1):
        raise ValueError('power form')
    if not (zero(lv3[0] - lv3[1]) and zero(lv3[2] - lv3[3])):
        raise ValueError('ends should coincide at level 3')
    res['power_form'] = {'chi_1': str(c1), 'chi_2': str(c2), 'chi_3': str(c3), 'level_3': 'the two ends of each axis coincide'}

    # T6 a potential that is not a power: chi_2 on both sides of 1
    Ug = a**2/2 + a*b/3 + b**2/2 + a**2*b/5 + b**3/7
    g1 = from_potential(Ug, a, b)
    g2 = [sp.simplify(q) for q in step(*g1)]
    c1g = sp.simplify(cross_ratio(*g1))
    c2g = sp.simplify(cross_ratio(*g2))
    seen = {}
    for sa, sb in itertools.product([sp.Rational(n, 2) for n in range(1, 7)], repeat=2):
        hess = sp.hessian(Ug, (a, b)).subs({a: sa, b: sb})
        if hess[0, 0] > 0 and hess.det() > 0:
            v1, v2 = c1g.subs({a: sa, b: sb}), c2g.subs({a: sa, b: sb})
            if not (0 < v1 <= 1):
                raise ValueError('level-1 ratio left (0, 1] inside the stable region')
            key = 'above' if v2 > 1 else 'below' if v2 < 1 else 'one'
            seen.setdefault(key, (str(sa), str(sb), str(v1), str(sp.nsimplify(v2))))
    if 'above' not in seen or 'below' not in seen:
        raise ValueError('expected both sides of 1')
    res['not_a_power'] = seen
    return res


if __name__ == '__main__':
    out = run()
    with open('UP3_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
