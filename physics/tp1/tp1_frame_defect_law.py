"""TP1: gravity in the order defect of the frame field; which quadratic law does the framework's equivalence select?

The frame field of the thesis carries a flat connection (its own: the one in which the frame is constant).
Every loop closes (information invariance in the sense of the paper), and the field is the ORDER DEFECT of the
frame, [e_a, e_b] = c_ab^k e_k   (EMK-C1's frame term; CV1: c = beta', beta/r for the fall frame).
The coordinate-independent quadratic laws in c are   Q = a1 I1 + a2 I2 + a3 I3 ,
    I1 = c_abk c^abk ,   I2 = c_ab^k c_k^b_a (first and one lower index exchanged) ,   I3 = c_b c^b , c_b = c_ab^a .
T1  The frame's own connection has zero curvature; its torsion is minus the order defect.
T2  Equivalence as a requirement (GR2-G1 made exact): the law must not depend on which frame is used at each
    place.  For flat space read in a frame turned or boosted differently at each place, Q must then be a pure
    boundary term.  This holds  iff  a1 : a2 : a3 = 1 : 2 : -4,  i.e. Q = k (I1/4 + I2/2 - I3).
T3  For the fall frame, with that combination,  Q = -R + (boundary term), R the contracted-twice curvature of CV1.
    So the selected law is the law of CV1.
T4  For other coefficients the profile m = r_s/r of the fall frame is kept only on a plane of coefficients;
    (1, 2, -4) lies on it.
Symbolic algebra with sympy.  Python 3.12.
"""
import json
import sys

import sympy as sp
from sympy.calculus.euler import euler_equations

sys.path.insert(0, '../cv1')
import cv1_frame_curvature as cv

ETA = cv.ETA
t, r, th, ph = cv.X
a1, a2, a3 = sp.symbols('a1 a2 a3')


def structure_in(coords, E):
    Einv = E.inv()
    c = [[[0]*4 for _ in range(4)] for _ in range(4)]
    def ap(vec, f):
        return sum(vec[i]*sp.diff(f, coords[i]) for i in range(4))
    for a in range(4):
        for b in range(4):
            comm = [ap(E.row(a), E[b, i]) - ap(E.row(b), E[a, i]) for i in range(4)]
            for k in range(4):
                c[a][b][k] = sp.simplify(sum(comm[i]*Einv[i, k] for i in range(4)))
    return c


def invariants(c):
    e = lambda i: ETA[i, i]
    I1 = sum(e(k)*e(a)*e(b)*c[a][b][k]**2 for k in range(4) for a in range(4) for b in range(4))
    I2 = sum(e(b)*c[a][b][k]*c[k][b][a] for k in range(4) for a in range(4) for b in range(4))
    tr = [sum(c[a][b][a] for a in range(4)) for b in range(4)]
    I3 = sum(e(b)*tr[b]**2 for b in range(4))
    return sp.simplify(I1), sp.simplify(I2), sp.simplify(I3)


def own_connection_is_flat(E, coords):
    """Gamma^rho_{mu nu} = E_a^rho d_nu theta^a_mu : curvature zero, torsion = -order defect."""
    theta = E.inv()            # theta[mu, a]: coframe components (columns a)
    G = [[[sp.simplify(sum(E[a, rho]*sp.diff(theta[mu, a], coords[nu]) for a in range(4))) for nu in range(4)] for mu in range(4)] for rho in range(4)]
    for rho in range(4):
        for mu in range(4):
            for al in range(4):
                for be in range(al+1, 4):
                    val = sp.diff(G[rho][mu][be], coords[al]) - sp.diff(G[rho][mu][al], coords[be])
                    val += sum(G[rho][s][al]*G[s][mu][be] - G[rho][s][be]*G[s][mu][al] for s in range(4))
                    if sp.simplify(val) != 0:
                        return False
    return True


def flat_space_in_a_turned_frame():
    """flat space; frame boosted along x by eta and turned in the x-y plane by alpha, both varying with place."""
    import random
    T, Xc, Y, Z = sp.symbols('T X Y Z', real=True)
    coords = (T, Xc, Y, Z)
    eta = sp.Function('eta')(T, Xc, Y, Z)
    alpha = sp.Function('alpha')(T, Xc, Y, Z)
    boost = sp.Matrix([[sp.cosh(eta), sp.sinh(eta), 0, 0], [sp.sinh(eta), sp.cosh(eta), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    turn = sp.Matrix([[1, 0, 0, 0], [0, sp.cos(alpha), sp.sin(alpha), 0], [0, -sp.sin(alpha), sp.cos(alpha), 0], [0, 0, 0, 1]])
    rng = random.Random(7)
    rows, cases = [], 0
    for E, fields in ((boost, [eta]), (turn, [alpha]), (boost*turn, [eta, alpha])):
        c = structure_in(coords, sp.simplify(E))
        I1, I2, I3 = invariants(c)
        Q = a1*I1 + a2*I2 + a3*I3          # the frame has unit volume
        for eq in euler_equations(Q, fields, list(coords)):
            expr = sp.simplify(eq.lhs - eq.rhs)
            cases += 1
            atoms = sorted(expr.atoms(sp.Derivative), key=str)
            for _ in range(8):             # the expression must vanish for every field: sample it
                sub = {d: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)) for d in atoms}
                val = expr.subs(sub).subs({eta: sp.Rational(rng.randint(-5, 5), 7), alpha: sp.Rational(rng.randint(-5, 5), 7)})
                rows.append([sp.diff(val, s) for s in (a1, a2, a3)])
    # exact: (1, 2, -4) annihilates every Euler expression identically; exact rank at eta = alpha = 0
    exact_zero, exact_rows = True, []
    for E, fields in ((boost, [eta]), (turn, [alpha]), (boost*turn, [eta, alpha])):
        c = structure_in(coords, sp.simplify(E))
        I1, I2, I3 = invariants(c)
        for eq in euler_equations(a1*I1 + a2*I2 + a3*I3, fields, list(coords)):
            expr = sp.simplify(eq.lhs - eq.rhs)
            if sp.simplify(expr.subs({a1: 1, a2: 2, a3: -4})) != 0:
                exact_zero = False
            atoms = sorted(expr.atoms(sp.Derivative), key=str)
            for _ in range(6):
                sub = {d: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)) for d in atoms}
                val = sp.simplify(expr.subs(sub).subs({eta: 0, alpha: 0}))
                exact_rows.append([sp.diff(val, s_) for s_ in (a1, a2, a3)])
    exact_rank = sp.Matrix(exact_rows).rank()
    import numpy as np
    A = np.array([[float(sp.N(v, 30)) for v in row] for row in rows])
    sv = sorted(np.linalg.svd(A, compute_uv=False).tolist(), reverse=True)
    null_direction = np.linalg.svd(A)[2][-1]
    null_direction = (null_direction/null_direction[0]).tolist()
    # exact direction: check that (1, 2, -4) is annihilated symbolically and that the rank is two
    vec = sp.Matrix([1, 2, -4])
    residual = max(abs(sp.N(sum(row[i]*vec[i] for i in range(3)), 30)) for row in rows)
    return dict(euler_expressions=cases, samples=len(rows), singular_values=[float(x) for x in sv],
                null_direction=null_direction, residual_on_1_2_minus4=float(residual),
                exact_zero_on_1_2_minus4=exact_zero, exact_rank_at_identity=exact_rank)


def fall_frame():
    beta = cv.beta
    E = cv.frame(beta)
    c = cv.structure(E)
    I1, I2, I3 = invariants(c)
    flat = own_connection_is_flat(E, cv.X)
    # T3: Q_GR + R is a boundary term  <=>  its Euler-Lagrange expression in beta vanishes identically
    w = cv.connection(c)
    R = cv.curvature(E, c, w)
    Ric = cv.contracted(R)
    scalar = sp.simplify(sum(ETA[i, i]*Ric[i, i] for i in range(4)))
    # pointwise: Q_GR + R = k x divergence of the trace vector of the order defect
    tr = [sum(c[a][b][a] for a in range(4)) for b in range(4)]
    vol = r**2*sp.sin(th)
    V = [sum(ETA[b, b]*tr[b]*E[b, mu] for b in range(4)) for mu in range(4)]
    div = sp.simplify(sum(sp.diff(vol*V[mu], cv.X[mu]) for mu in range(4))/vol)
    Qgr0 = sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3
    identity = {k: str(sp.simplify(Qgr0 + scalar - k*div)) for k in (2, -2)}
    bfun = sp.Function('b')(r)
    def reduced(expr):
        return sp.simplify((r**2*expr).subs(beta, bfun).doit())
    Qgr = sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3
    el_sum = sp.simplify(euler_equations(reduced(Qgr + scalar), [bfun], [r])[0].lhs)
    el_minus = sp.simplify(euler_equations(reduced(Qgr - scalar), [bfun], [r])[0].lhs)
    # T4: general coefficients, the profile m = r_s / r
    rs = sp.Symbol('r_s', positive=True)
    Qgen = a1*I1 + a2*I2 + a3*I3
    el = euler_equations(reduced(Qgen), [bfun], [r])[0].lhs
    on_profile = sp.simplify(el.subs(bfun, sp.sqrt(rs/r)).doit())
    plane = sp.factor(sp.simplify(on_profile*r**sp.Rational(3, 2)/sp.sqrt(rs)))
    return dict(own_connection_flat=flat, I1=str(I1), I2=str(I2), I3=str(I3), scalar_curvature=str(scalar),
                euler_of_Q_plus_R=str(el_sum), euler_of_Q_minus_R=str(el_minus), Q_plus_R_minus_k_div=identity,
                profile_condition=str(plane),
                profile_condition_at_selected=str(sp.simplify(plane.subs({a1: sp.Rational(1, 4), a2: sp.Rational(1, 2), a3: -1}))))


def run():
    return dict(turned_frame=flat_space_in_a_turned_frame(), fall_frame=fall_frame())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('TP1_RESULT.json', 'w'), indent=1)
    print(res['turned_frame'])
    for k, v in res['fall_frame'].items():
        print(k, ':', v)
