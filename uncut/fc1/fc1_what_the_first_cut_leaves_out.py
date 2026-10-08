"""FC1: what a first cut that is its own dagger leaves out, and the mass relation with that part restored.

The owner's question: today's physics and mathematics stand on a first cut made at lambda = 0.  Was something
missed at that cut - a sine, a vertical component, a dimension - by the frame or by flatness?  E = m c^2 also
assumes flatness; is it the complete equation?

A reading is rho = a + bK + cS + dR on the primitive carrier; a change of frame acts as rho -> G rho G-dagger
(native dagger: K, S kept, R reversed), det G = 1.
Sources read (unchanged): Publications emk-ugd-algebra EMK-1; theorum/48, 55 (relations, dagger); RMG9 T3 (quadratic
readings blind to w, loops read 2w Area); extra-ideas physics in1 (n^2 - r.r), gb1 (rho -> M rho M-dagger), dm1;
uncut up4, up7, tt1, ci1.
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


def dagger(M):
    a, b, c, d = parts(M)
    return a*I2 + b*K + c*S - d*R


def zero(e):
    return sp.simplify(sp.expand(e)) == 0


def run():
    res = {}
    a, b, c, d = sp.symbols('a b c d', real=True)
    rho = a*I2 + b*K + c*S + d*R
    p, q, r = sp.symbols('p q r', real=True)
    G = sp.Matrix([[p, q], [r, (1 + q*r)/p]])                 # every frame change with p != 0: det = 1
    if not zero(G.det() - 1) or dagger(G) != G.T:
        raise ValueError('frame change')

    # T1 under every change of frame: the turn-part d is kept, and a^2 - b^2 - c^2 is kept, separately
    a2, b2, c2, d2 = parts(sp.simplify(G*rho*dagger(G)))
    if not zero(d2 - d):
        raise ValueError('turn-part is not invariant')
    if not zero(a2*a2 - b2*b2 - c2*c2 - (a*a - b*b - c*c)):
        raise ValueError('flat mass relation is not invariant')
    if not zero(rho.det() - (a*a - b*b - c*c) - d*d):
        raise ValueError('determinant')
    res['invariants'] = 'two, separately: m^2 = a^2 - b^2 - c^2 (the self-dagger part) and W = d (the turn-part); det = m^2 + W^2'

    # T2 a boost mixes a with the momentum and never touches W; at rest det = M^2 gives a = M cos, W = M sin
    eta = sp.Symbol('eta', real=True)
    B = sp.cosh(eta/2)*I2 + sp.sinh(eta/2)*K
    ab, bb, cb, db = parts(sp.simplify((B*rho*dagger(B)).applyfunc(lambda e: sp.expand(e.rewrite(sp.exp)))))
    if not (sp.simplify((ab - (a*sp.cosh(eta) + b*sp.sinh(eta))).rewrite(sp.exp)) == 0
            and sp.simplify((bb - (b*sp.cosh(eta) + a*sp.sinh(eta))).rewrite(sp.exp)) == 0 and zero(cb - c) and zero(db - d)):
        raise ValueError('boost')
    Mm, th = sp.symbols('M theta', positive=True)
    rest = Mm*(sp.cos(th)*I2 + sp.sin(th)*R)
    if not zero(sp.simplify(rest.det() - Mm**2)) or parts(rest)[0] != Mm*sp.cos(th) or parts(rest)[3] != Mm*sp.sin(th):
        raise ValueError('rest element')
    res['rest'] = 'energy read at rest = M cos(theta) ; turn-part = M sin(theta) ; M^2 = E_0^2 + W^2'

    # T3 no quadratic reading sees W; a loop does
    x, y = sp.symbols('x y', real=True)
    X = sp.Matrix([x, y])
    if not zero((X.T*rho*X)[0, 0] - (X.T*(rho - d*R)*X)[0, 0]):
        raise ValueError('quadratic reading')
    J = rho*X
    curl = sp.diff(J[1], x) - sp.diff(J[0], y)
    if not zero(curl - 2*d):
        raise ValueError('loop reading')
    res['who_sees_it'] = 'X.rho X does not contain W ; the return around a loop is 2 W x Area'

    # T4 W adds like a charge: for two readings, det(rho_1 + rho_2) = m_12^2 + (W_1 + W_2)^2
    a1, b1, c1, d1, a2_, b2_, c2_, d2_ = sp.symbols('a1 b1 c1 d1 a2 b2 c2 d2', real=True)
    r1 = a1*I2 + b1*K + c1*S + d1*R
    r2 = a2_*I2 + b2_*K + c2_*S + d2_*R
    m12 = (a1 + a2_)**2 - (b1 + b2_)**2 - (c1 + c2_)**2
    if not zero((r1 + r2).det() - m12 - (d1 + d2_)**2):
        raise ValueError('pair')
    res['pairs'] = 'turn-parts add; opposite turn-parts cancel and the pair is flat'

    # T5 exact witness: M = 5, W = 3: energy at rest 4; boosted with speed 3/5: energy 5, momentum 3; W still 3
    w0 = 4*I2 + 3*R
    Bq = sp.Matrix([[sp.sqrt(2), 0], [0, 1/sp.sqrt(2)]])       # exp(eta) = 2 : cosh = 5/4, sinh = 3/4
    wb = sp.simplify(Bq*w0*dagger(Bq))
    pa = parts(wb)
    if not (pa == (5, 3, 0, 3) and wb.det() == 25 and pa[0]**2 - pa[1]**2 == 16):
        raise ValueError('witness: ' + str(pa))
    res['witness'] = {'at_rest': 'E = 4, W = 3, M = 5', 'moving': 'E = 5, p = 3, W = 3 ; E^2 - p^2 = 16 ; det = 25'}
    return res


if __name__ == '__main__':
    out = run()
    with open(__file__.replace('fc1_what_the_first_cut_leaves_out.py', 'FC1_RESULT.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
