"""GW1: waves of the frame's order defect -- a pair of readings, content two, and the same count.

Sources: TP1 (law Q = I1/4 + I2/2 - I3 on the order defect of the frame; selected by equivalence),
WQ1 (a wave mode is a pair; recordable <=> area = 2 pi kappa n <=> E = n kappa omega), QC1-T3 (content q silent
iff q Theta in 2 pi Z; the covariance has content 2), DC1-K1 / R43.3 (ledgers that exchange share one unit),
EM1 (the light wave), EG1/TL1 (energy is exchanged between the sides).
Frame: flat, plus a small transverse turn-free change travelling along z,
    e1 = (1 - a/2) d_x - (b/2) d_y ,   e2 = -(b/2) d_x + (1 + a/2) d_y ,   e0 = d_t ,   e3 = d_z ,
with a(t,z), b(t,z) of first order.
V1  To second order the law is   Q2 = (1/2) [ (a_t)^2 - (a_z)^2 + (b_t)^2 - (b_z)^2 ] , pointwise:
    a difference of two squares, like E^2 - B^2 for light; the rate term is positive.
V2  Stationarity gives  a_tt - a_zz = 0  and the same for b: two waves, at the cone speed, uncoupled.
V3  Turning the frame about z by theta turns (a, b) by 2 theta: content two.  The pattern returns after half a turn.
V4  One mode a = A(t) cos(kz): the pair (A, dA/dt) is an oscillator with omega = k and energy
    (1/4) [ (dA/dt)^2 + k^2 A^2 ] averaged over a period in z: a pair in the sense of WQ1.
V5  Hence (WQ1) a recordable mode has E = n kappa omega; and since these waves exchange energy with the other
    sides, R43.3 gives them the same unit kappa: E = n h nu with the h of light.
Symbolic algebra with sympy.  Python 3.12.
"""
import json
import sys

import sympy as sp
from sympy.calculus.euler import euler_equations

sys.path.insert(0, '../tp1')
sys.path.insert(0, '../cv1')
import tp1_frame_defect_law as tp

T, X, Y, Z = sp.symbols('T X Y Z', real=True)
COORDS = (T, X, Y, Z)
eps = sp.Symbol('epsilon')


def wave_frame(a, b):
    return sp.Matrix([[1, 0, 0, 0],
                      [0, 1 - eps*a/2, -eps*b/2, 0],
                      [0, -eps*b/2, 1 + eps*a/2, 0],
                      [0, 0, 0, 1]])


def second_order_law(a, b):
    c = tp.structure_in(COORDS, wave_frame(a, b))
    I1, I2, I3 = tp.invariants(c)
    Q = sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3
    ser = sp.series(sp.simplify(Q), eps, 0, 3).removeO()
    return sp.simplify(ser.coeff(eps, 0)), sp.simplify(ser.coeff(eps, 1)), sp.simplify(ser.coeff(eps, 2))


def run():
    a = sp.Function('a')(T, Z)
    b = sp.Function('b')(T, Z)
    q0, q1, q2 = second_order_law(a, b)
    if q0 != 0 or q1 != 0:
        raise ValueError('the law must start at second order on a flat frame')
    form = (sp.diff(a, T)**2 - sp.diff(a, Z)**2 + sp.diff(b, T)**2 - sp.diff(b, Z)**2)/2
    sign = None
    for s in (1, -1):
        # equal up to a boundary term  <=>  same Euler-Lagrange expressions
        diff_el = [sp.simplify(e.lhs) for e in euler_equations(q2 - s*form, [a, b], [T, Z])]
        if all(x == 0 for x in diff_el):
            sign = s
    if sign != 1:
        raise ValueError('second-order law is not a difference of squares of the two rates with a positive rate term')
    pointwise = sp.simplify(q2 - sign*form)
    el = [sp.simplify(e.lhs) for e in euler_equations(sign*q2, [a, b], [T, Z])]
    want = [sp.simplify(-(sp.diff(a, T, 2) - sp.diff(a, Z, 2))), sp.simplify(-(sp.diff(b, T, 2) - sp.diff(b, Z, 2)))]
    if [sp.simplify(x - y) for x, y in zip(el, want)] != [0, 0]:
        raise ValueError('stationarity is not the wave equation at the cone speed')
    # V3: turn the frame about z by theta
    th = sp.Symbol('theta', real=True)
    a0, b0 = sp.symbols('a0 b0')
    S = sp.Matrix([[-a0/2, -b0/2], [-b0/2, a0/2]])          # the change of (e1, e2)
    R = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
    S2 = sp.simplify(R*S*R.T)
    a1, b1 = sp.simplify(-2*S2[0, 0]), sp.simplify(-2*S2[0, 1])
    turn2 = (sp.simplify(a1 - (a0*sp.cos(2*th) - b0*sp.sin(2*th))) == 0 and
             sp.simplify(b1 - (a0*sp.sin(2*th) + b0*sp.cos(2*th))) == 0)
    half_turn = sp.simplify(S2.subs(th, sp.pi) - S) == sp.zeros(2, 2)
    quarter = sp.simplify(S2.subs(th, sp.pi/2) + S) == sp.zeros(2, 2)
    if not (turn2 and half_turn and quarter):
        raise ValueError('content two failed')
    # V4: one mode
    k = sp.Symbol('k', positive=True)
    A = sp.Function('A')(T)
    mode = (sign*q2).subs(b, 0).doit().subs(a, A*sp.cos(k*Z)).doit()
    Lz = 2*sp.pi/k
    L_mode = sp.simplify(sp.integrate(sp.expand(sign*form.subs(b, 0).doit().subs(a, A*sp.cos(k*Z)).doit()), (Z, 0, Lz))/Lz)
    want_L = (sp.diff(A, T)**2 - k**2*A**2)/4
    if sp.simplify(L_mode - want_L) != 0:
        raise ValueError('mode Lagrangian failed')
    p = sp.diff(L_mode, sp.diff(A, T))
    energy = sp.simplify(p*sp.diff(A, T) - L_mode)
    if sp.simplify(energy - (sp.diff(A, T)**2 + k**2*A**2)/4) != 0:
        raise ValueError('mode energy failed')
    # the pair of WQ1: w = A sqrt(k/2), p_w = (dA/dt)/sqrt(2k) gives E = (k/2)(w^2 + p_w^2), omega = k
    w = A*sp.sqrt(k/2)
    pw = sp.diff(A, T)/sp.sqrt(2*k)
    if sp.simplify(energy - (k/2)*(w**2 + pw**2)) != 0:
        raise ValueError('pair form failed')
    return dict(order0=str(q0), order1=str(q1), sign=sign, pointwise_difference=str(pointwise),
                wave_equations=[str(x) for x in el], content=2, returns_after_half_turn=bool(half_turn),
                mode_energy=str(energy), mode_rate='k')


def illustration():
    import math
    h, G, c = 6.62607015e-34, 6.67430e-11, 299792458.0
    f, strain = 100.0, 1e-21
    flux = c**3/(16*math.pi*G)*0.5*(2*math.pi*f*strain)**2       # one polarisation, amplitude 'strain'
    return dict(frequency_Hz=f, energy_per_count_J=h*f, flux_W_per_m2=flux, counts_per_m2_per_s=flux/(h*f))


if __name__ == '__main__':
    res = run()
    res['illustration'] = illustration()
    json.dump(res, open('GW1_RESULT.json', 'w'), indent=1)
    for kk, v in res.items():
        print(kk, ':', v)
