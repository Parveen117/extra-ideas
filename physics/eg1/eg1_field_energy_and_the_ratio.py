"""EG1: the energy of the light field makes gravity the same way heat does; the ratio law with a source.

Sources: CV1 (contracted curvature of the fall frame), NC1 (native ratio rho), EM1/OB1 (field F = E + iota B,
reading-type products F rho F^dagger), IN1 (reading tensor), TL1 (units follow the clock), QC1 (content q is a count).
G1  Reading of a radial field F = E n.C through the cuts:  (1/2) F F^dagger = u ,  (1/2) F C_a F^dagger = u (2 n_a n.C - C_a):
    energy u = E^2/2, pull along the field, push across it: the pattern (u ; -u, +u, +u).  Its invariant part is zero.
G2  Memory m = r_s/r - q2/r^2 has contracted curvature (q2/r^4) x (1 ; -1, 1, 1): exactly the pattern of G1 with
    u ~ 1/r^4.  Field energy with zero invariant part gravitates (the correction to R3 is confirmed).
G3  For the fall-frame family the energy component of the law is  (r m)'/r^2 = k u  for ANY energy density u, and
        rho + 1/2 = (r m)' / (2 m) = k u r^2 / (2 m) :
    the excess of the native ratio over -1/2 is the energy density, of whatever kind, in units of m / r^2.
G4  Potential follows the clock: content q is a count, so the potential V = energy/q obeys V_A N_A = V_B N_B (TL1-L2).
G5  The response element of the light field in a medium, Hessian of u(D,B) = D^2/2eps + B^2/2mu:
    determinant = 1/(eps mu) = (speed)^2, eigenvalue ratio = mu/eps = (impedance)^2.
Exact rational arithmetic for G1, G4, G5; sympy for G2, G3.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

import sympy as sp

sys.path.insert(0, '../in1')
sys.path.insert(0, '../cv1')
import in1_the_invariant as in1
import cv1_frame_curvature as cv

c, cmul, cadd, conj = in1.c, in1.cmul, in1.cadd, in1.conj
ZC = (F(0), F(0))
CUTS = [[[c(0), c(1)], [c(1), c(0)]], [[c(0), c(0, -1)], [c(0, 1), c(0)]], [[c(1), c(0)], [c(0), c(-1)]]]
ONE = [[c(1), c(0)], [c(0), c(1)]]


def mm(a, b):
    return [[cadd(cmul(a[i][0], b[0][j]), cmul(a[i][1], b[1][j])) for j in range(2)] for i in range(2)]


def lin(coeffs, mats):
    out = [[ZC, ZC], [ZC, ZC]]
    for k, m in zip(coeffs, mats):
        out = [[cadd(out[i][j], cmul(c(k), m[i][j])) for j in range(2)] for i in range(2)]
    return out


def dagger(a):
    return [[conj(a[j][i]) for j in range(2)] for i in range(2)]


def stress_pattern():
    rows = []
    for n, E in (((F(3, 13), F(4, 13), F(12, 13)), F(5)), ((F(2, 3), F(-1, 3), F(2, 3)), F(7, 2)), ((F(0), F(0), F(1)), F(1))):
        if sum(x*x for x in n) != 1:
            raise ValueError('direction must be a unit')
        nC = lin(n, CUTS)
        Fm = lin([E], [nC])
        u = E*E/2
        half = lambda M: [[cmul(c(F(1, 2)), x) for x in row] for row in M]
        if half(mm(Fm, dagger(Fm))) != lin([u], [ONE]):
            raise ValueError('energy reading failed')
        for a in range(3):
            got = half(mm(mm(Fm, CUTS[a]), dagger(Fm)))
            want = lin([2*u*n[a]] + [(-u if b == a else F(0)) for b in range(3)], [nC] + CUTS)
            if got != want:
                raise ValueError('stress reading failed')
        # eigen-pattern: along n the cut reading is +u n.C (pull -u in the stress), across it -u C_perp (push +u)
        along = half(mm(mm(Fm, nC), dagger(Fm)))
        if along != lin([u], [nC]):
            raise ValueError('pull along the field failed')
        # invariant (trace) part: energy reading plus the scalar part of sum_a C_a (1/2 F C_a F^dagger) must vanish
        total = half(mm(Fm, dagger(Fm)))
        for a in range(3):
            total = [[cadd(x, y) for x, y in zip(r1, r2)] for r1, r2 in zip(total, mm(CUTS[a], half(mm(mm(Fm, CUTS[a]), dagger(Fm)))))]
        scalar_part = cadd(total[0][0], total[1][1])
        if scalar_part != ZC:
            raise ValueError('invariant (trace) part must vanish')
        rows.append((str(E), str(u)))
    return rows


def curvature_with_field_energy():
    r = cv.r
    rs, q2, k = sp.symbols('r_s q2 k', positive=True)
    m = sp.Function('m')(r)
    E = cv.frame(cv.beta)
    cst = cv.structure(E)
    Ric = cv.contracted(cv.curvature(E, cst, cv.connection(cst)))
    def through_m(expr, mm_):
        e = sp.simplify(expr)
        e = e.subs(sp.diff(cv.beta, r, 2), (sp.diff(mm_, r, 2)/2 - sp.diff(cv.beta, r)**2)/cv.beta)
        e = e.subs(sp.diff(cv.beta, r), sp.diff(mm_, r)/(2*cv.beta))
        return sp.simplify(sp.simplify(e.subs(cv.beta**2, mm_)).subs(cv.beta, sp.sqrt(mm_)).doit())
    charged = rs/r - q2/r**2
    diag = [sp.simplify(through_m(Ric[i, i], charged)) for i in range(4)]
    want = [q2/r**4, -q2/r**4, q2/r**4, q2/r**4]
    if any(sp.simplify(a - b) != 0 for a, b in zip(diag, want)):
        raise ValueError('contracted curvature of the charged memory is not the field-energy pattern')
    # G3: energy component and the ratio
    Rm = [through_m(Ric[i, i], m) for i in range(4)]
    scalar = Rm[0] - Rm[1] - Rm[2] - Rm[3]
    energy_component = sp.simplify(Rm[0] - scalar/2)
    if sp.simplify(energy_component - sp.diff(r*m, r)/r**2) != 0:
        raise ValueError('energy component is not (r m)\' / r^2')
    rho_m = sp.simplify(r*sp.diff(m, r)/(2*m))                    # NC1: rho = d ln tanh(eta) / d ln r, tanh^2(eta) = m
    excess = sp.simplify(rho_m + sp.Rational(1, 2) - sp.diff(r*m, r)/(2*m))
    if excess != 0:
        raise ValueError('rho + 1/2 is not (r m)\'/(2 m)')
    rho_charged = sp.simplify((r*sp.diff(charged, r)/(2*charged)) + sp.Rational(1, 2))
    return dict(contracted=[str(x) for x in diag], energy_component=str(energy_component),
                ratio_excess_charged=str(sp.factor(rho_charged)))


def potential_and_response():
    rows = []
    for na, nb, va in ((F(3, 5), F(4, 5), F(12)), (F(5, 13), F(12, 13), F(7, 3))):
        q = 3                                   # content: a count, the same at both places (QC1)
        energy_a = q*va
        energy_b = energy_a*na/nb               # TL1-L1/L2
        vb = energy_b/q
        if vb*nb != va*na:
            raise ValueError('potential must follow the clock')
        rows.append((str(va), str(vb)))
    resp = []
    for eps, mu in ((F(2), F(8)), (F(9, 4), F(1)), (F(1), F(1))):
        a, cc = 1/eps, 1/mu                     # Hessian of u(D, B)
        det, ratio = a*cc, a/cc
        if det != 1/(eps*mu) or ratio != mu/eps:
            raise ValueError('response dictionary failed')
        resp.append((str(eps), str(mu), str(det), str(ratio)))
    return dict(potential=rows, response=resp)


def illustration():
    G_N, eps0, cl, e = 6.67430e-11, 8.8541878128e-12, 299792458.0, 1.602176634e-19
    import math
    rq = lambda Q: math.sqrt(G_N*Q*Q/(4*math.pi*eps0))/cl**2
    return dict(length_of_one_electron_charge_m=rq(e), length_of_one_coulomb_m=rq(1.0),
                electron_r_s_m=2*G_N*9.1093837e-31/cl**2)


def run():
    return dict(stress=stress_pattern(), curvature=curvature_with_field_energy(), dictionary=potential_and_response(),
                illustration=illustration())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('EG1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
