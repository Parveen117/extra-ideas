"""UP4: the cost of assuming that the two cut operations commute.

UP1-UP3 used partial derivatives in a chart: d_S d_V = d_V d_S.  That is a flatness assumption.  Here the two
operations D_S, D_V are only assumed to satisfy  [D_S, D_V] = alpha D_S + beta D_V  and Leibniz; everything is
recomputed and the terms that the assumption had hidden are written out.

No measurement, no physical constant.  Sources read (unchanged):
  response-geometry series      DB1 (derivation algebras [D_a, D_b] = f_abc D_c: identities from bracket + Leibniz only)
                                RMG9 T1, T3 (L = m(1 + uK + vS) + wR; det L = det L_sym + w^2; cycle 2w Area)
  Recognition-Kernel-Framework  theorum/55 T3, theorum/57 (the loop of two flows leaves h^2 [D_1, D_2] to leading order)
  Publications                  cut-first-equivalence (classical = memoryless sector; loop residue is the obstruction)
  extra-ideas                   uncut/up1, up2, up3; physics/tp1 (order defect of a frame)
A concrete model is used only to check the identities: D_S = d_a, D_V = f d_b + g d_a with free f, g.
Symbolic algebra with sympy.  Python 3.12.
"""
import json

import sympy as sp

a, b = sp.symbols('a b', real=True)
f = sp.Function('f')(a, b)
g = sp.Function('g')(a, b)


def DS(h, ff=f, gg=g):
    return sp.diff(h, a)


def DV(h, ff=f, gg=g):
    return ff*sp.diff(h, b) + gg*sp.diff(h, a)


def defect(ff=f, gg=g):
    """[D_S, D_V] = alpha D_S + beta D_V"""
    beta = sp.diff(ff, a)/ff
    alpha = sp.diff(gg, a) - gg*beta
    return alpha, beta


def zero(e):
    return sp.simplify(e) == 0


def response(U, ff=f, gg=g):
    T, P = DS(U, ff, gg), -DV(U, ff, gg)
    A = DS(T, ff, gg)
    B1 = DV(T, ff, gg)          # change of T along the V operation
    B2 = -DS(P, ff, gg)         # change of -P along the S operation
    C = -DV(P, ff, gg)
    return T, P, A, B1, B2, C


def run():
    res = {}
    h = sp.Function('h')(a, b)
    alpha, beta = defect()
    if not zero(DS(DV(h)) - DV(DS(h)) - alpha*DS(h) - beta*DV(h)):
        raise ValueError('bracket relation')

    # T1 the mixed responses differ by the defect times the first layer
    U = sp.Function('U')(a, b)
    T, P, A, B1, B2, C = response(U)
    if not zero((B2 - B1) - (alpha*T - beta*P)):
        raise ValueError('mixed responses')
    res['mixed'] = 'B_2 - B_1 = [D_S, D_V] U = alpha T - beta P = 2w'

    # T2 the ratio of ends is still shared by the two axes, but 1 - chi is no longer a square
    s, v = sp.symbols('s v')                      # a displacement: s along D_S, v along D_V
    dT, dP = A*s + B1*v, -B2*s - C*v
    slope_T_at_P = sp.simplify((dT.subs(v, sp.solve(dP, v)[0]))/s)
    slope_T_at_V = A
    slope_P_at_T = sp.simplify((dP.subs(s, sp.solve(dT, s)[0]))/v)
    slope_P_at_S = -C
    chi_x = sp.simplify(slope_T_at_P/slope_T_at_V)
    chi_y = sp.simplify(slope_P_at_T/slope_P_at_S)
    det = A*C - B1*B2
    if not (zero(chi_x - det/(A*C)) and zero(chi_y - det/(A*C))):
        raise ValueError('ratio of ends')
    Bbar, w = (B1 + B2)/2, (B2 - B1)/2
    if not zero(1 - chi_x - (Bbar**2 - w**2)/(A*C)):
        raise ValueError('1 - chi')
    if not zero(det - (A*C - Bbar**2) - w**2):
        raise ValueError('determinant split')
    res['ratio'] = '1 - chi = (Bbar^2 - w^2)/(AC) ; det = det(symmetric part) + w^2 ; both axes still share chi'

    # T3 first-order statements of UP3 do not use the assumption
    Tq, Sq, Pq, Vq = [sp.Function(n)(a, b) for n in 'TSPV']
    brD = lambda x, y: DS(x)*DV(y) - DV(x)*DS(y)
    three = brD(Tq, Sq)*brD(Pq, Vq) - brD(Tq, Pq)*brD(Sq, Vq) + brD(Tq, Vq)*brD(Sq, Pq)
    if not zero(sp.expand(three)):
        raise ValueError('three-term identity')
    res['survives'] = 'ratio of ends shared by both axes; three-term identity'

    # T4 a witness: positive symmetric response, chi > 1
    ff, gg = sp.exp(a), sp.Integer(0)             # [D_S, D_V] = D_V : alpha = 0, beta = 1
    al, be = defect(ff, gg)
    if not (al == 0 and sp.simplify(be) == 1):
        raise ValueError('witness defect')
    Uw = a**2 + b**2 - a*b/2 + 2*b
    vals = [sp.simplify(q.subs({a: 0, b: 0})) for q in response(Uw, ff, gg)]
    Tw, Pw, Aw, B1w, B2w, Cw = vals
    chi_w = (Aw*Cw - B1w*B2w)/(Aw*Cw)
    sym_det = Aw*Cw - ((B1w + B2w)/2)**2
    if not (chi_w == sp.Rational(19, 16) and Aw > 0 and sym_det > 0 and (B2w - B1w) == -Pw):
        raise ValueError('witness')
    res['witness'] = {'A': str(Aw), 'C': str(Cw), 'B1': str(B1w), 'B2': str(B2w), 'P': str(Pw), 'chi': str(chi_w),
                      'symmetric_part': 'positive'}
    # with the first layer absent (P = 0) the same defect costs nothing at that point
    U0 = a**2 + b**2 - a*b/2
    v0 = [sp.simplify(q.subs({a: 0, b: 0})) for q in response(U0, ff, gg)]
    if not (v0[1] == 0 and v0[3] == v0[4]):
        raise ValueError('cost should vanish with the first layer')

    # T5 the loop of the two operations does not close; the potential changes by eps^2 [D_S, D_V] U
    e = sp.Symbol('epsilon')
    a1 = a + e                                   # along D_S
    b2 = b + e*sp.exp(a1)                        # along D_V = exp(a) d_b at the current a
    a3 = a1 - e                                  # back along D_S
    b4 = b2 - e*sp.exp(a3)                       # back along D_V
    gap_b = sp.series(b4 - b, e, 0, 3).removeO()
    if not zero(gap_b - e**2*sp.exp(a)):
        raise ValueError('gap of the loop')
    Uloop = Uw.subs({a: a3, b: b4}, simultaneous=True)
    change = sp.series(Uloop - Uw, e, 0, 3).removeO()
    comm = DS(DV(Uw, ff, gg), ff, gg) - DV(DS(Uw, ff, gg), ff, gg)
    if not zero(sp.expand(change - e**2*comm)):
        raise ValueError('loop residue')
    res['loop'] = 'four steps of size eps end eps^2 [D_S, D_V] away; the potential differs by eps^2 (alpha T - beta P) = 2w eps^2'
    return res


if __name__ == '__main__':
    out = run()
    with open('UP4_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
