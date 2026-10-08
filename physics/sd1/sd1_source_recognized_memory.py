"""SD1: T24 Theorem 6.1, S = R + D and F = R - D, on the three objects of the physics line.

T24: an event lift Z, a cut pair P + Q = 1; S = Z*Z (source), R = Z*PZ (recognized), D = Z*QZ (memory), F = R - D.
All four are Gram forms: quadratic.  Here Z is (a) a record of turns (GE2-T2), (b) a single reading (IN1),
(c) the stretch block of the fall frame (SE1); (d) the order defect of any frame (GF1).
Exact rational arithmetic and sympy.  Python 3.12."""
import json, os, sys
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1'))
sys.path.insert(0, os.path.join(HERE, '..', 'gf1'))
_cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'tp1'))
import cv1_frame_curvature as cv
import tp1_frame_defect_law as tp
os.chdir(_cwd)
import gf1_no_memory_in_every_frame as gf
z = lambda e: sp.simplify(e) == 0


def run():
    out = {}
    I = sp.I
    # (a) record of turns: carrier = weighted direct sum over records s; Z f = (u_s f)_s ; P = J J*, J f = (f)_s , J*(f_s) = sum p_s f_s
    us = [sp.Rational(3, 5) + sp.Rational(4, 5)*I, sp.Rational(5, 13) + sp.Rational(12, 13)*I, I]       # a + bR, R as the unit
    ps = [sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 6)]
    pair = lambda f, g_: sum(p*sp.conjugate(a)*b for p, a, b in zip(ps, f, g_))                          # the weighted pairing
    Zf = [u for u in us]                                                                                # Z applied to f = 1
    mean = sum(p*u for p, u in zip(ps, us))
    PZf = [mean for _ in us]                                                                            # P Z f = J J* Z f
    QZf = [a - b for a, b in zip(Zf, PZf)]
    S, R, D = pair(Zf, Zf), pair(Zf, PZf), pair(Zf, QZf)
    out['a_cut_pair'] = z(pair(PZf, QZf)) and z(pair(PZf, PZf) - R) and z(pair(QZf, QZf) - D)            # P, Q orthogonal projections
    out['a_S_is_R_plus_D'] = z(S - 1) and z(S - R - D) and sp.simplify(R) >= 0 and sp.simplify(D) >= 0
    out['a_D_is_GE2_memory'] = z(D - (1 - sp.Abs(mean)**2)) and z(R - sp.expand(mean*sp.conjugate(mean)))
    # (b) a single reading psi on the cut-complex carrier; cut i: P_i = (1 + C_i)/2
    C = [sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]])]
    psi = sp.Matrix([sp.Rational(2, 3) - I, 3])
    n = (psi.H*psi)[0]
    ok, Fs = True, []
    for Ci in C:
        P = (sp.eye(2) + Ci)/2; Qp = sp.eye(2) - P
        ok &= (P*P - P == sp.zeros(2)) and (P*Qp == sp.zeros(2))
        Ri, Di = sp.simplify((psi.H*P*psi)[0]), sp.simplify((psi.H*Qp*psi)[0])
        Fi = sp.simplify((psi.H*Ci*psi)[0])
        ok &= z(n - Ri - Di) and z(Fi - (Ri - Di)) and Ri >= 0 and Di >= 0
        Fs.append((Ri, Di, Fi))
    out['b_S_is_R_plus_D_each_cut'] = ok
    # IN1-T2: S^2 = sum F_i^2 for a single reading; and 4 R_i D_i = sum of the other cuts' F^2: lost to one cut = observed by the others
    out['b_single_reading'] = z(n**2 - sum(f[2]**2 for f in Fs)) and all(z(4*Fs[i][0]*Fs[i][1] - sum(Fs[j][2]**2 for j in range(3) if j != i)) for i in range(3))
    # (c) the stretch block K: Z = K, cut i: P_i = e_i e_i^T
    K = sp.Matrix(3, 3, lambda i, j: sp.Symbol('k%d%d' % (min(i, j), max(i, j))))
    ok, Ss, As = True, [], []
    for i in range(3):
        e = sp.zeros(3, 1); e[i] = 1
        P = e*e.T; Qp = sp.eye(3) - P
        Si = sp.expand((e.T*K.T*K*e)[0]); Ri = sp.expand((e.T*K.T*P*K*e)[0]); Di = sp.expand((e.T*K.T*Qp*K*e)[0])
        ok &= z(Si - Ri - Di) and z(Ri - K[i, i]**2) and z(Di - sum(K[j, i]**2 for j in range(3) if j != i))
        Ss.append(Si); As.append(K[i, i])
    out['c_S_is_R_plus_D_each_cut'] = ok                                                                  # = LN1-N1
    e2 = sum(K[i, i]*K[j, j] - K[i, j]**2 for i in range(3) for j in range(i + 1, 3))
    # the energy law of SE1: 2 e2 = (sum of signed readings)^2 - sum of sources ; empty space: sum S_i = (tr Z)^2
    out['c_law_in_S_and_trace'] = z(sp.expand(2*e2 - (sum(As)**2 - sum(Ss))))
    out['c_total_memory_is_cross_terms'] = z(sp.expand(2*e2 - (2*sum(As[i]*As[j] for i in range(3) for j in range(i + 1, 3)) - sum(Ss[i] - As[i]**2 for i in range(3)))))
    # (d) any frame: Q = G - I3 , G = I1/4 + I2/2 blind to turn rates, I3 = (uncut, trace reading)^2
    c, syms = gf.general_defect()
    I1, I2, I3 = tp.invariants(c)
    G = sp.expand(sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2)
    w1, w2, w3 = sp.symbols('w1 w2 w3')
    sub = {s: 0 for s in syms.values()}
    sub.update({syms[(0, 1, 2)]: w3, syms[(0, 2, 1)]: -w3, syms[(0, 2, 3)]: w1, syms[(0, 3, 2)]: -w1, syms[(0, 3, 1)]: w2, syms[(0, 1, 3)]: -w2})
    tr = [sum(c[a][b][a] for a in range(4)) for b in range(4)]
    out['d_Gram_part_blind_to_turn'] = z(G.subs(sub)) and z(sp.expand(I3 - sum(cv.ETA[b, b]*tr[b]**2 for b in range(4))))
    # on one plane the Gram part equals the square of the trace reading: Q = 0
    keep = {syms[(0, 1, 0)], syms[(0, 1, 1)]}
    one = {s: 0 for s in syms.values() if s not in keep}
    out['d_one_plane_G_equals_trace_square'] = z(sp.expand((G - I3).subs(one)))
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out), open(os.path.join(HERE, 'SD1_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
