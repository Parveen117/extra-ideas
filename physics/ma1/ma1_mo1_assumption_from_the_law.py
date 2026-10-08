"""MA1: MO1's assumption (the locally pure frames share one flat space and one time) from the line's own law.

Family: the general static, spherically symmetric frame in falling form, with a free radial stretch A(r):
    e^0 = dt ,  e^1 = A(r) (dr + b(r) dt) ,  e^2 = r dtheta ,  e^3 = r sin(theta) dphi .
MO1 assumed A = 1. Here A is left free and the law of CV1/TP1 (contracted curvature of the frame zero in empty
space; TP1's selected combination = -R + boundary term) is applied, with the functions of those stages.
Symbolic (sympy)."""
import json, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1')); sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
import tp1_frame_defect_law as tp
cv = tp.cv
t, r, th, ph = cv.X
ETA = cv.ETA


def run():
    out = {}
    z = lambda e: sp.simplify(e) == 0
    A = sp.Function('A')(r); b = sp.Function('b')(r)
    E = sp.Matrix([[1, -b, 0, 0], [0, 1/A, 0, 0], [0, 0, 1/r, 0], [0, 0, 0, 1/(r*sp.sin(th))]])   # rows: frame vectors
    # T1: no loss: a static frame N dT, S dr (N < 1) is this family with A = N S, b^2 = (1 - N^2)/(N S)^2
    N, S_ = sp.symbols('N S', positive=True)
    g_ = S_*sp.sqrt(1 - N**2)/N                       # dt = dT + g dr
    A0 = N*S_; b0 = N**2*g_/A0**2
    out['T1_general_static_frame'] = z(1 - A0**2*b0**2 - N**2) and z(A0**2*b0 - N**2*g_) and z(N**2*g_**2 - S_**2 + A0**2)
    # T2: TP1's identity holds on the whole family: (1/4) I1 + (1/2) I2 - I3 = -R - 2 div(trace of the order defect)
    c = tp.structure_in(cv.X, E)
    I1, I2, I3 = tp.invariants(c)
    Q = sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3
    w = cv.connection(c); Rm = cv.curvature(E, c, w); Ric = cv.contracted(Rm)
    scal = sp.simplify(sum(ETA[i, i]*Ric[i, i] for i in range(4)))
    tr = [sum(c[a][k][a] for a in range(4)) for k in range(4)]
    vol = A*r**2*sp.sin(th)
    V = [sum(ETA[k, k]*tr[k]*E[k, mu] for k in range(4)) for mu in range(4)]
    div = sum(sp.diff(vol*V[mu], cv.X[mu]) for mu in range(4))/vol
    out['T2_TP1_identity_on_the_family'] = z(Q + scal + 2*div)
    # T3: the mixed component of the contracted curvature is -2 b A'/(r A^2): zero in empty space  <=>  A' = 0
    out['T3_mixed_component'] = z(Ric[0, 1] + 2*b*sp.diff(A, r)/(r*A**2))
    # T4: with A constant every component vanishes exactly when  A^2 (r b^2)' = 1 - A^2
    a0 = sp.symbols('a0', positive=True)
    cond = a0**2*sp.diff(r*b**2, r) - (1 - a0**2)
    RicA = Ric.subs(A, a0).doit().applyfunc(sp.simplify)
    out['T4_angular_component'] = z(RicA[2, 2] - cond/(r**2*a0**2))
    rs = sp.symbols('r_s', positive=True)
    bsol = sp.sqrt((1 - a0**2)/a0**2 + rs/r)
    out['T4_all_components_vanish_on_the_solution'] = all(z(RicA[i, j].subs(b, bsol).doit()) for i in range(4) for j in range(4))
    # T5: far condition: the frame is at rest far away (b -> 0)  =>  A = 1  =>  flat slices, one time, and b^2 = r_s / r
    far = sp.limit(bsol**2, r, sp.oo)
    out['T5_far_condition_gives_A_1'] = sp.solve(sp.Eq(far, 0), a0) == [1]
    out['T5_profile'] = z(bsol.subs(a0, 1)**2 - rs/r)
    # T6: otherwise the same field read by frames with a speed far away: b_far^2 = (1 - A^2)/A^2
    out['T6_other_members'] = z(far - (1 - a0**2)/a0**2)
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(out, open(os.path.join(HERE, 'MA1_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
