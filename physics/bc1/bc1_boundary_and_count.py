"""BC1: is the exact counting at the boundary of the tower the counting of waves?

Sources: RMG2-T1(c) ((R H)^2 = -Delta 1: complex inside the disc, dual on the boundary), RMG2 holonomy
(Hol = Exp(Theta J_H), J_H = R H / sqrt(Delta)), LT1 (the tower ends on the boundary; invariant -> 0),
WQ1 (a wave is a pair; area in whole units; E = n kappa omega; count frame-independent), QC1-T2/T3 (return = bracket;
content q silent iff q Theta in 2 pi Z), TL1-L1 (a turn carried between places scales by the clock factors),
GR1 (N^2 = 1 - m), SY1 (three sectors).
B1  (R H)^2 = -Delta.  With the un-normalised generator the return is
        Exp(s R H) = cos(s d) 1 + (sin(s d)/d) R H ,   d = sqrt(Delta):
    a turn through the angle s d.  A whole turn needs s = 2 pi / d.
B2  On the boundary (Delta = 0):  Exp(s R H) = 1 + s R H.  Returns add, (1 + s1 N)(1 + s2 N) = 1 + (s1 + s2) N,
    and no s other than 0 comes back to 1: there is no whole turn, hence no silent content and no step.
B3  For the gravity element at clock factor N the trace-normalised invariant is Delta_n = N^2, so d = N:
    the turn per unit of far parameter is the clock factor (TL1-L1).  The whole-turn unit, seen from far, is
    2 pi / N; the step of a wave's energy, seen from far, is kappa omega N.  Both lose their finiteness at N = 0.
B4  The bracket of a pair is itself additive: in the three-by-three unipotent form the loop
        (1 + a E12)(1 + b E23)(1 - a E12)(1 - b E23) = 1 + a b E13
    returns the enclosed area in a central, additive slot.  Whole numbers appear only when that additive
    quantity is read as a turn (B1).
Conclusion: counting is an additive quantity (dual sector) read through a turn (circular sector).  The interior
supplies the period, the boundary supplies the addition; on the boundary alone the period is infinite.
Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json

R = [[F(0), F(-1)], [F(1), F(0)]]
ID = [[F(1), F(0)], [F(0), F(1)]]


def mul(a, b):
    n = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def lin(x, a, y, b):
    return [[x*a[i][j] + y*b[i][j] for j in range(len(a))] for i in range(len(a))]


def turn(c, s, d, RH):
    """cos 1 + (sin/d) R H"""
    return lin(c, ID, s/d, RH)


def run():
    # B1
    rows = []
    for (a, b, cc), d in (((F(5), F(3), F(5)), F(4)), ((F(13, 4), F(3), F(13, 4)), F(5, 4)), ((F(2), F(0), F(1, 2)), F(1))):
        H = [[a, b], [b, cc]]
        delta = a*cc - b*b
        if delta != d*d:
            raise ValueError('fixture: Delta must be a rational square')
        RH = mul(R, H)
        if mul(RH, RH) != lin(-delta, ID, F(0), ID):
            raise ValueError('(R H)^2 = -Delta failed')
        pts = ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(-7, 25), F(24, 25)))
        for (c1, s1) in pts:
            for (c2, s2) in pts:
                lhs = mul(turn(c1, s1, d, RH), turn(c2, s2, d, RH))
                rhs = turn(c1*c2 - s1*s2, s1*c2 + c1*s2, d, RH)
                if lhs != rhs:
                    raise ValueError('returns must compose as turns')
        if turn(F(1), F(0), d, RH) != ID or turn(F(-1), F(0), d, RH) != lin(F(-1), ID, F(0), ID):
            raise ValueError('whole and half turn failed')
        rows.append((str(delta), str(d)))
    # B2
    Hb = [[F(9, 5), F(12, 5)], [F(12, 5), F(16, 5)]]            # rank one: Delta = 0
    if Hb[0][0]*Hb[1][1] - Hb[0][1]**2 != 0:
        raise ValueError('boundary fixture failed')
    Nb = mul(R, Hb)
    if mul(Nb, Nb) != [[F(0), F(0)], [F(0), F(0)]]:
        raise ValueError('boundary generator must square to zero')
    shear = lambda s: lin(F(1), ID, s, Nb)
    for s1, s2 in ((F(1, 3), F(5, 7)), (F(-2), F(9, 4))):
        if mul(shear(s1), shear(s2)) != shear(s1 + s2):
            raise ValueError('boundary returns must add')
    never = all(shear(F(k, 7)) != ID for k in range(1, 400))
    if not never:
        raise ValueError('a boundary return must not come back')
    # B3
    clock = []
    kappa, omega = F(1, 2), F(6)
    for beta in (F(3, 5), F(4, 5), F(12, 13), F(99, 101)):
        Hn = [[1 + beta, F(0)], [F(0), 1 - beta]]                # trace-normalised unit block in its principal frame
        delta_n = Hn[0][0]*Hn[1][1]
        N2 = 1 - beta*beta
        if delta_n != N2:
            raise ValueError('normalised invariant must be the clock factor squared')
        root = {F(16, 25): F(4, 5), F(9, 25): F(3, 5), F(25, 169): F(5, 13), F(400, 10201): F(20, 101)}[N2]
        step_far = kappa*omega*root                               # WQ1 step carried to far away (TL1-L1)
        clock.append((str(beta), str(root), str(step_far)))
    if not all(F(a[2]) > F(b[2]) for a, b in zip(clock, clock[1:])):
        raise ValueError('the far step must shrink toward the boundary')
    # B4
    def E(i, j):
        m = [[F(0)]*3 for _ in range(3)]
        m[i][j] = F(1)
        return m
    I3 = [[F(int(i == j)) for j in range(3)] for i in range(3)]
    add = lambda x, y: [[x[i][j] + y[i][j] for j in range(3)] for i in range(3)]
    sc = lambda k, m: [[k*v for v in row] for row in m]
    areas = []
    for a, b in ((F(2), F(3)), (F(1, 2), F(-5, 3)), (F(7), F(1, 7))):
        A, Am = add(I3, sc(a, E(0, 1))), add(I3, sc(-a, E(0, 1)))
        B, Bm = add(I3, sc(b, E(1, 2))), add(I3, sc(-b, E(1, 2)))
        loop = mul(mul(mul(A, B), Am), Bm)
        if loop != add(I3, sc(a*b, E(0, 2))):
            raise ValueError('loop of a pair must return its area in the centre')
        areas.append(str(a*b))
    # the central slot adds
    C1, C2 = add(I3, sc(F(3), E(0, 2))), add(I3, sc(F(5, 2), E(0, 2)))
    if mul(C1, C2) != add(I3, sc(F(11, 2), E(0, 2))):
        raise ValueError('areas must add')
    return dict(turn_fixtures=rows, boundary_adds=True, boundary_never_returns=never, clock=clock, loop_areas=areas)


if __name__ == '__main__':
    res = run()
    json.dump(res, open('BC1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
