"""UP7: why the diagram exists - a turn and one self-dagger cut generate it; nothing thermodynamic is used.

Sources read (unchanged):
  Recognition-Kernel-Framework theorum/morphic_algebra/main.tex   EMK-0 (cut kappa, curvature iota, recognition chi = iota kappa),
        EMK-3 (4-cycle T -kappa-> V -iota-> S -kappa-> P -iota-> T; status CONDITIONAL), Maxwell corollary (status INCOMPLETE),
        EMK-Rep (the representation fixed for computations), EMK-12 open question 3 (potentials from the 4-cycle)
  Recognition-Kernel-Framework theorum/48, theorum/55             certified relations K^2 = I, R^2 = -I, RK = -KR; dagger R -> -R, K -> K
  response-geometry RMG2 T1(a), RMG9 T1                           X = pK + qS + rR, X^2 = (p^2 + q^2 - r^2)
  extra-ideas uncut/up4 (w), up5 (four relations are one equation), up6 (eight scale operations)
Exact rational matrices on the primitive carrier (two real dimensions).  sympy.  Python 3.12.
"""
import itertools
import json

import sympy as sp

I2 = sp.eye(2)
R = sp.Matrix([[0, -1], [1, 0]])          # the turn: R^2 = -1
K = sp.Matrix([[1, 0], [0, -1]])          # a self-dagger cut
S = R*K                                   # the other one, a quarter turn of the mirror away


def dagger(M):
    """native dagger on the EMK span: coefficients of 1, K, S kept, coefficient of R reversed"""
    a = (M[0, 0] + M[1, 1])/2
    u = (M[0, 0] - M[1, 1])/2
    v = (M[0, 1] + M[1, 0])/2
    w = (M[1, 0] - M[0, 1])/2
    return a*I2 + u*K + v*S - w*R


def cut(u, v, w):
    return u*K + v*S + w*R


def group_generated(gens):
    seen, frontier = [I2], [I2]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                x = g*h
                if not any(x == y for y in seen):
                    seen.append(x)
                    new.append(x)
        if len(seen) > 64:
            raise ValueError('group does not close')
        frontier = new
    return seen


def run():
    res = {}
    u, v, w = sp.symbols('u v w', real=True)

    # T1 every cut on the carrier: X = uK + vS + wR with u^2 + v^2 - w^2 = 1; its relation to the turn
    X = cut(u, v, w)
    if sp.simplify(X*X - (u*u + v*v - w*w)*I2) != sp.zeros(2, 2):
        raise ValueError('square of a traceless element')
    if sp.simplify(X*R + R*X + 2*w*I2) != sp.zeros(2, 2):
        raise ValueError('anticommutator with the turn')
    if sp.simplify(dagger(X) - X + 2*w*R) != sp.zeros(2, 2):
        raise ValueError('dagger')
    res['cuts'] = 'kappa = uK + vS + wR, u^2 + v^2 - w^2 = 1 ; kappa R + R kappa = -2w ; kappa-dagger = kappa - 2wR'

    # T2 the four-step walk: chi = iota kappa ; (chi + w)^2 = 1 + w^2 ; four steps miss by -2w chi
    for turn in (R, -R):
        chi = turn*X
        s = 1 if turn == R else -1
        on_shell = {u: sp.sqrt(1 + w*w - v*v)}
        if sp.simplify(((chi + s*w*I2)**2 - (1 + w*w)*I2).subs(on_shell)) != sp.zeros(2, 2):
            raise ValueError('walk square')
        if sp.simplify((chi*chi - I2 + 2*s*w*chi).subs(on_shell)) != sp.zeros(2, 2):
            raise ValueError('four-step miss')
    res['walk'] = '(iota kappa)^2 = 1 - 2w (iota kappa) (sign of w follows the sense of the turn): closes iff w = 0'

    # T3 w = 0: the group has eight elements, four mirrors; the compass
    turn, kap = -R, S                      # clockwise turn; the cut is the mirror through the F and H corners
    G = group_generated([turn, kap])
    mirrors = [g for g in G if g*g == I2 and g != I2 and g != -I2]
    if len(G) != 8 or len(mirrors) != 4:
        raise ValueError('group order')
    T = sp.Matrix([0, 1])
    V, Sv, P = kap*T, turn*kap*T, kap*turn*kap*T
    if not (V == sp.Matrix([1, 0]) and Sv == sp.Matrix([0, -1]) and P == sp.Matrix([-1, 0]) and turn*P == T):
        raise ValueError('compass walk')
    ends = {'T': T, 'V': V, 'S': Sv, 'P': P}
    corners = {'F': T + V, 'U': Sv + V, 'H': Sv + P, 'G': T + P}
    name = lambda vec, table: [k for k, x in table.items() if x == vec][0]
    chi = turn*kap
    leg_thermal, leg_mech, both = chi, kap*chi*kap, -I2
    acts = {}
    for label, g in (('thermal', leg_thermal), ('mechanical', leg_mech), ('both', both), ('axis exchange', kap)):
        acts[label] = {'ends': ''.join(name(g*x, ends) for x in ends.values()),
                       'corners': ''.join(name(g*x, corners) for x in corners.values())}
    want = {'thermal': {'ends': 'SVTP', 'corners': 'UFGH'}, 'mechanical': {'ends': 'TPSV', 'corners': 'GHUF'},
            'both': {'ends': 'SPTV', 'corners': 'HGFU'}, 'axis exchange': {'ends': 'VTPS', 'corners': 'FGHU'}}
    if acts != want:
        raise ValueError('action on the compass: ' + json.dumps(acts))
    res['group'] = {'order': 8, 'mirrors': 4,
                    'images_of_T_V_S_P_and_F_U_H_G': acts}

    # T4 the eight flags (end, neighbouring end) = the eight scale operations of UP6: one orbit, no stabiliser
    flags = [(a, b) for a, b in itertools.permutations(ends, 2) if ends[a] + ends[b] != sp.zeros(2, 1)]
    img = {(name(g*ends['T'], ends), name(g*ends['V'], ends)) for g in G}
    if len(flags) != 8 or img != set(flags):
        raise ValueError('flags')
    res['flags'] = 'eight, permuted simply transitively by the eight elements'

    # T5 w != 0: the walk never returns
    wv = sp.Rational(3, 4)
    kw = cut(sp.Rational(5, 4), 0, wv)
    if kw*kw != I2:
        raise ValueError('witness is not a cut')
    cw = R*kw
    ev = sorted(cw.eigenvals().keys())
    if ev != [sp.Integer(-2), sp.Rational(1, 2)]:
        raise ValueError('eigenvalues')
    if (cw**2 - I2) != -2*wv*cw or any((cw**n == I2 or cw**n == -I2) for n in range(1, 13)):
        raise ValueError('witness walk')
    res['open_walk'] = {'w': '3/4', 'eigenvalues_of_iota_kappa': ['-2', '1/2'], 'growth_per_four_steps': '4 and 1/4'}

    # T6 the response element has its own cut; it is self-dagger exactly when the mixed responses agree
    A, B1, B2, C = sp.symbols('A B_1 B_2 C', real=True)
    L = sp.Matrix([[A, B1], [B2, C]])
    Xl = L - (A + C)/2*I2
    D = ((A - C)/2)**2 + B1*B2
    if sp.simplify(Xl*Xl - D*I2) != sp.zeros(2, 2):
        raise ValueError('response cut')
    if sp.simplify(Xl*R + R*Xl + (B2 - B1)*I2) != sp.zeros(2, 2):
        raise ValueError('response cut against the turn')
    res['response'] = 'L = mean + sqrt(D) kappa_L ; kappa_L R + R kappa_L = -(B_2 - B_1)/sqrt(D) : self-dagger iff B_1 = B_2'

    # T7 the representation the Morphic manuscript fixes for computations does not give a four-cycle
    km = sp.Matrix([[1, 1], [0, 0]])
    cm = R*km
    stated = sp.Matrix([[0, -1], [1, 1]])                     # the matrix the manuscript prints for iota kappa
    if not (km*km == km and cm*cm == cm and cm != I2):
        raise ValueError('manuscript representation')
    if not (stated != cm and stated**3 == -I2 and stated**2 != I2 and stated**2 != -I2):
        raise ValueError('manuscript printed matrix')
    res['manuscript_representation'] = ('kappa^2 = kappa and (iota kappa)^2 = iota kappa : the walk stops, it does not return ; '
                                        'the matrix printed there for iota kappa is a different one, of order six')
    return res


if __name__ == '__main__':
    out = run()
    with open('UP7_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
