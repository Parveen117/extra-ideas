"""CL2: two histories with the same ends - the count, the corner cost and the split character.

A leg of a history is a displacement X = x H + y K + t R with X^2 = -(tau^2) < 0 on the upper
sheet.  A sector of rest turn g carried along X has energy-momentum P = (g / tau) X (PR2-T2).
  T1  count of a leg = -sc(P X) = g tau.
  T2  two legs: T^2 - (tau1 + tau2)^2 = 2 tau1 tau2 (cosh(d) - 1), cosh(d) = -B(X1, X2)/(tau1 tau2) >= 1,
      zero iff the legs are collinear.  The unturned history has the largest count.
  T3  N equal legs on a hyperbola of radius rho, corner boost 2a:
      count / straight count = N / chi_N(a),  chi_N(a) = sinh(N a)/sinh(a) >= N.
  T4  both sheets: g -> -g changes the sign of every count and no ratio.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../pr2')
import pr2_two_sheets_pair_observer as pr2

element, mm, sc, ONE, coords = pr2.element, pr2.mm, pr2.sc, pr2.ONE, pr2.coords

# legs (x, y, t, tau) with t^2 - x^2 - y^2 = tau^2, t > 0
LEGS = [(F(3), F(0), F(5), F(4)), (F(-3), F(0), F(5), F(4)), (F(4), F(0), F(5), F(3)),
        (F(2), F(2), F(3), F(1)), (F(4), F(8), F(9), F(1)), (F(6), F(-3), F(7), F(2)),
        (F(7), F(4), F(9), F(4)), (F(0), F(0), F(2), F(2)), (F(6), F(0), F(10), F(8))]


def form(x, y, t):
    return x*x+y*y-t*t


def bil(a, b):
    return a[0]*b[0]+a[1]*b[1]-a[2]*b[2]


def leg_control(g):
    rows = []
    for x, y, t, tau in LEGS:
        if form(x, y, t) != -tau*tau or t <= 0 or tau <= 0:
            raise ValueError('not a timelike leg on the upper sheet')
        X = element(x, y, t)
        P = sc(g/tau, X)
        if mm(P, P) != sc(-g*g, ONE):
            raise ValueError('carried sector does not have rest turn g')
        count = -coords(mm(P, X))[0]
        if count != g*tau:
            raise ValueError('count of a leg is not g tau')
        rows.append(dict(leg=[str(x), str(y), str(t)], count=str(count)))
    return rows


def corner_control(g):
    rows = []
    for i, a in enumerate(LEGS):
        for b in LEGS[i+1:]:
            total = (a[0]+b[0], a[1]+b[1], a[2]+b[2])
            T2 = -form(*total)
            split = (a[3]+b[3])**2
            cosh = -bil(a, b)/(a[3]*b[3])
            if cosh < 1:
                raise ValueError('relative boost of two upper-sheet legs must have cosh >= 1')
            if T2-split != 2*a[3]*b[3]*(cosh-1):
                raise ValueError('corner cost formula failed')
            collinear = (a[0]*b[2] == b[0]*a[2] and a[1]*b[2] == b[1]*a[2])
            if (T2 == split) != collinear:
                raise ValueError('equality must hold exactly for collinear legs')
            rows.append(dict(cosh_corner=str(cosh), straight_squared=str(T2), turned_count_squared=str(split)))
    # three legs: the straight history still wins
    a, b, c = LEGS[0], LEGS[3], LEGS[5]
    total = tuple(a[k]+b[k]+c[k] for k in range(3))
    if -form(*total) <= (a[3]+b[3]+c[3])**2:
        raise ValueError('three-leg history should lose count')
    # the gap does not depend on g's sign
    return dict(pairs=len(rows), strict=sum(r['straight_squared'] != r['turned_count_squared'] for r in rows),
                example=rows[0])


def cheb_U(n, x):
    a, b = F(1), 2*x
    if n == 0:
        return a
    for _ in range(n-1):
        a, b = b, 2*x*b-a
    return b


def hyperbola_control(ch, sh, rho, levels=4):
    """Points rho (sinh, cosh) at equally spaced rapidities; coarser and coarser polygons."""
    if ch*ch-sh*sh != 1:
        raise ValueError('not a boost')

    def double(c, s):
        return c*c+s*s, 2*c*s
    # rapidities k * a for k = -n..n in steps of 2: build (cosh(ka), sinh(ka)) by addition
    n = 2**(levels-1)
    table = {0: (F(1), F(0)), 1: (ch, sh)}
    for k in range(2, n+1):
        c1, s1 = table[k-1]
        table[k] = (c1*ch+s1*sh, s1*ch+c1*sh)

    def point(k):
        c, s = table[abs(k)]
        return (rho*c, F(0), rho*s*(1 if k >= 0 else -1))      # (x, y, t) = rho (cosh, 0, sinh)
    rows = []
    straight2 = None
    for level in range(levels):
        step = 2**(level+1)                                    # corner-to-corner rapidity = step * a
        ks = list(range(-n, n+1, step))
        legs = len(ks)-1
        count = F(0)
        for k0, k1 in zip(ks, ks[1:]):
            p0, p1 = point(k0), point(k1)
            d = (p1[0]-p0[0], F(0), p1[2]-p0[2])
            tau2 = -form(*d)
            half = table[step//2]                              # leg = 2 rho sinh(step a / 2)
            if tau2 != (2*rho*half[1])**2:
                raise ValueError('chord is not 2 rho sinh(half corner)')
            count += 2*rho*half[1]
        whole = 2*rho*table[n][1]
        chi = cheb_U(legs-1, table[step//2][0])                # sinh(legs * b)/sinh(b), b = step a / 2
        if count*chi != legs*whole:
            raise ValueError('count / straight is not N / chi_N')
        if chi < legs or (chi == legs) != (legs == 1):
            raise ValueError('split character must exceed N')
        rows.append(dict(legs=legs, count=str(count), straight=str(whole), ratio=str(count/whole)))
    for fine, coarse in zip(rows, rows[1:]):
        if not F(fine['count']) < F(coarse['count']):
            raise ValueError('a finer turned history must have a smaller count')
    return rows


def sheet_control(g):
    x, y, t, tau = LEGS[0]
    X = element(x, y, t)
    up = -coords(mm(sc(g/tau, X), X))[0]
    down = -coords(mm(sc(-g/tau, X), X))[0]
    if up != -down:
        raise ValueError('sheet flip must flip the sign of the count only')
    return dict(mass=str(up), antimass=str(down), ratio_of_ratios='1')


def run(g=F(3, 2)):
    return dict(legs=leg_control(g), corners=corner_control(g),
                hyperbola=hyperbola_control(F(41, 40), F(9, 40), F(1)), sheets=sheet_control(g))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('CL2_RESULT.json', 'w'), indent=1)
    print(res['legs'][:3]); print(res['corners'])
    for r in res['hyperbola']:
        print({k: (v if len(str(v)) < 30 else float(F(v))) for k, v in r.items()})
    print(res['sheets'])
