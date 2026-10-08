"""MM1: one memory.  The record memory of turns (GE2-T2), of readings (IN1-T5) and of the frame's stretch (SE1) is one
quadratic form: e2, the determinant in a plane (EMK-1), taken as  mean of e2 - e2 of the mean.
Gravity's law on flat slices is: no such memory between the cuts, except the energy of a source.

Exact rational arithmetic and sympy, with the functions of CV1 and TP1.  Python 3.12."""
import json, os, sys
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1'))
_cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'tp1'))
import cv1_frame_curvature as cv
import tp1_frame_defect_law as tp
os.chdir(_cwd)
z = lambda e: sp.simplify(e) == 0


def e2(M):
    return sp.expand((M.trace()**2 - (M*M).trace())/2)


def run():
    out = {}
    # U1: for any records X_s with weights of sum one:  mean of e2 - e2 of the mean = sum_{s<t} p_s p_t e2(X_s - X_t)
    ok = True
    for n in (2, 3):
        Xs = [sp.Matrix(n, n, lambda i, j: sp.Symbol('x%d_%d%d' % (s, i, j))) for s in range(3)]
        p1, p2 = sp.symbols('p1 p2')
        ps = [p1, p2, 1 - p1 - p2]
        mean = sum((ps[s]*Xs[s] for s in range(3)), sp.zeros(n))
        lhs = sum(ps[s]*e2(Xs[s]) for s in range(3)) - e2(mean)
        rhs = sum(ps[s]*ps[t]*e2(Xs[s] - Xs[t]) for s in range(3) for t in range(s + 1, 3))
        ok &= z(sp.expand(lhs - rhs))
    out['U1_memory_identity'] = ok
    # e2 of a 2x2 block is its determinant (EMK-1)
    A2 = sp.Matrix(2, 2, sp.symbols('m0:4'))
    out['U1_e2_is_det_in_a_plane'] = z(e2(A2) - A2.det())
    # U2: turns a + bR as blocks: e2 = det = 1, so the memory is 1 - det S : GE2-T2 as used in OR1
    turn = lambda a, b: sp.Matrix([[a, -b], [b, a]])
    us = [turn(sp.Rational(3, 5), sp.Rational(4, 5)), turn(sp.Rational(5, 13), sp.Rational(12, 13)), turn(0, 1)]
    ps = [sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 6)]
    S = sum((ps[s]*us[s] for s in range(3)), sp.zeros(2))
    out['U2_turns'] = all(e2(u) == 1 for u in us) and (sum(ps[s]*e2(us[s]) for s in range(3)) - e2(S) == 1 - S.det())
    # U3: readings rho = psi psi^dagger on the cut-complex carrier (IN1): e2 = det = 0, memory = det of the record = (n^2 - r.r)/4
    I = sp.I
    C = [sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]])]      # K, RK, iota R
    def reading(psi):
        psi = sp.Matrix(psi); rho = psi*psi.H
        n = rho.trace(); r = [sp.simplify((rho*Ci).trace()) for Ci in C]
        return rho, n, r
    r1, n1, v1 = reading([1, sp.Rational(1, 2) + I])
    r2, n2, v2 = reading([sp.Rational(2, 3) - I, 3])
    p = sp.Rational(2, 5)
    rec = p*r1 + (1 - p)*r2
    nn = p*n1 + (1 - p)*n2; rr = [p*a + (1 - p)*b for a, b in zip(v1, v2)]
    out['U3_single_reading_has_none'] = z(e2(r1)) and z(e2(r2)) and z(n1**2 - sum(a**2 for a in v1))
    out['U3_record_of_readings'] = z(e2(rec) - (nn**2 - sum(a**2 for a in rr))/4) and z(e2(rec) - p*(1 - p)*(n1*n2 - sum(a*b for a, b in zip(v1, v2)))/2)
    # U4: stretch: K = sum lam_i P_i, P_i = single-cut readings (any directions): e2(P) = 0, e2(K) = sum_{i<j} lam_i lam_j (1 - overlap^2)
    lam = sp.symbols('l1:4', real=True)
    vs = [sp.Matrix([1, 2, 2]), sp.Matrix([2, -1, 3]), sp.Matrix([0, 1, -4])]
    Ps = [v*v.T/(v.dot(v)) for v in vs]
    Kc = sum((lam[i]*Ps[i] for i in range(3)), sp.zeros(3))
    ov2 = lambda i, j: (vs[i].dot(vs[j]))**2/(vs[i].dot(vs[i])*vs[j].dot(vs[j]))
    out['U4_single_cut_has_none'] = all(z(e2(P)) for P in Ps)
    out['U4_stretch_as_record'] = z(e2(Kc) - sum(lam[i]*lam[j]*(1 - ov2(i, j)) for i in range(3) for j in range(i + 1, 3)))
    # U5: LN1's form: 2 e2(K) = (sum of seen)^2 - sum of the first tower layer ; trace and tower: tr T(K) - T(tr K) = -2 lambda e2
    Ks = sp.Matrix(3, 3, lambda i, j: sp.Symbol('k%d%d' % (min(i, j), max(i, j))))
    seen = [Ks[i, i] for i in range(3)]
    layer = [(Ks*Ks)[i, i] for i in range(3)]
    lost = [sp.expand(layer[i] - seen[i]**2) for i in range(3)]
    out['U5_seen_and_lost'] = z(2*e2(Ks) - (sum(seen)**2 - sum(layer))) and z(2*e2(Ks) - (2*sum(seen[i]*seen[j] for i in range(3) for j in range(i + 1, 3)) - sum(lost)))
    lm = sp.Symbol('lambda')
    T_ = lambda H: H + lm*H*H
    out['U5_trace_commutes_with_the_tower'] = z((T_(Ks)).trace() - (Ks.trace() + lm*Ks.trace()**2) + 2*lm*e2(Ks))
    # U6: the three invariants of TP1 on flat slices: Q = A tr K^2 + B (turn)^2 + G (tr K)^2
    T, x, y, zz = sp.symbols('T x y z', real=True); co = (T, x, y, zz)
    u = [sp.Function('u%d' % i)(x, y, zz) for i in (1, 2, 3)]
    E = sp.Matrix([[1, u[0], u[1], u[2]], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    c = tp.structure_in(co, E)
    I1, I2, I3 = tp.invariants(c)
    a1, a2, a3 = sp.symbols('a1 a2 a3')
    d_ = lambda i, j: sp.diff(u[i], co[j + 1])
    K = sp.Matrix(3, 3, lambda i, j: (d_(i, j) + d_(j, i))/2)
    Om = sp.Matrix(3, 3, lambda i, j: (d_(i, j) - d_(j, i))/2)
    trK2, Om2, tr2 = (K*K).trace(), -(Om*Om).trace(), K.trace()**2
    A, B, G = 2*a1 + a2, 2*a1 - a2, a3
    out['U6_general_law_on_flat_slices'] = z(sp.expand(a1*I1 + a2*I2 + a3*I3 - (A*trK2 + B*Om2 + G*tr2)))
    # no memory for a pure turn (K = 0): B = 0 ; no memory for a stretch along a single cut (tr K^2 = (tr K)^2): A + G = 0
    sol = sp.solve([B, A + G], [a2, a3], dict=True)[0]
    out['U6_two_conditions_select_TP1'] = (sol[a2]/a1 == 2) and (sol[a3]/a1 == -4)
    plane = json.load(open(os.path.join(HERE, '..', 'tp1', 'TP1_RESULT.json')))['fall_frame']['profile_condition']
    out['U6_single_cut_condition_is_TP1_T4_plane'] = z(sp.sympify(plane)/(sp.Rational(9, 2)*sp.Symbol('r')) - (A + G))
    # U7: the volume of a cell of the falling frame:  V''/V = theta^2 + u.grad(theta) = 2 e2(K) - R_00
    w = cv.connection(c)
    ap = lambda vec, f: sum(vec[i]*sp.diff(f, co[i]) for i in range(4))
    def riem(a_, b_, k_, l_):
        val = ap(E.row(k_), w[a_][b_][l_]) - ap(E.row(l_), w[a_][b_][k_])
        val += sum(w[a_][e][k_]*w[e][b_][l_] - w[a_][e][l_]*w[e][b_][k_] for e in range(4))
        val -= sum(c[k_][l_][e]*w[a_][b_][e] for e in range(4))
        return val
    R00 = sp.simplify(sum(riem(a_, 0, a_, 0) for a_ in range(4)))
    th = K.trace()
    acc = th**2 + sum(u[i]*sp.diff(th, co[i + 1]) for i in range(3))
    e2K = sum(K[i, i]*K[j, j] - K[i, j]**2 for i in range(3) for j in range(i + 1, 3))
    out['U7_volume_acceleration'] = z(sp.expand(acc - 2*e2K + R00))
    # radial fall beta^2 = r_s/r: a cell labelled by r0 has volume factor linear in the frame's time
    t, r0, rs = sp.symbols('t r0 r_s', positive=True)
    r = (r0**sp.Rational(3, 2) - sp.Rational(3, 2)*sp.sqrt(rs)*t)**sp.Rational(2, 3)
    out['U7_fall_law'] = z(sp.diff(r, t) + sp.sqrt(rs/r))
    vol = sp.simplify(r**2*sp.diff(r, r0))
    out['U7_volume_linear_in_time'] = z(sp.diff(vol, t, 2)) and not z(sp.diff(vol, t))
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out), open(os.path.join(HERE, 'MM1_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
