"""PH3: the return angle is produced by sharing a probe between two flat readings.

E reading: the probe's extensive offsets (ds, dv) are frozen      -> transport 0.
M reading: the probe's intensive offsets (dT, -dP) are frozen     -> transport H^-1 dH.
Share lambda of M: A = lambda H^-1 dH, curvature -lambda(1-lambda)[X_s, X_v].
Exact part: stdlib rationals.  Fluid part: numpy + CoolProp.  Python 3.12.
"""
from fractions import Fraction as F
import json
import math
import sys


# ------------------------------------------------------------- exact part
def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def add(a, b, s=1):
    return [[a[i][j]+s*b[i][j] for j in range(2)] for i in range(2)]


def scal(c, a):
    return [[c*x for x in row] for row in a]


def inv(a):
    d = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[a[1][1]/d, -a[0][1]/d], [-a[1][0]/d, a[0][0]/d]]


def field(x, y):
    """Hessian of U = x^4 + y^4 + x^2 y^2 + x^3 y/3 + 2x^2 + 3y^2 and its derivatives."""
    H = [[12*x*x+2*y*y+2*x*y+4, 4*x*y+x*x], [4*x*y+x*x, 12*y*y+2*x*x+6]]
    Hx = [[24*x+2*y, 4*y+2*x], [4*y+2*x, 4*x]]
    Hy = [[4*y+2*x, 4*x], [4*x, 24*y]]
    Hxy = [[F(2), F(4)], [F(4), F(0)]]
    return H, Hx, Hy, Hxy


def exact_controls():
    rows = []
    for x, y in ((F(1, 2), F(1, 3)), (F(-2, 5), F(3, 4)), (F(1), F(-1, 7))):
        H, Hx, Hy, Hxy = field(x, y)
        if H[0][0] <= 0 or H[0][0]*H[1][1]-H[0][1]**2 <= 0:
            raise ValueError('response is not positive')
        Hi = inv(H)
        X, Y = mm(Hi, Hx), mm(Hi, Hy)
        dyX = add(mm(Hi, Hxy), mm(mm(Hi, Hy), mm(Hi, Hx)), -1)
        dxY = add(mm(Hi, Hxy), mm(mm(Hi, Hx), mm(Hi, Hy)), -1)
        comm = add(mm(X, Y), mm(Y, X), -1)
        flat = add(add(dxY, dyX, -1), comm)
        if any(v for row in flat for v in row):
            raise ValueError('intensive reading is not flat')
        if not any(v for row in comm for v in row):
            raise ValueError('witness has no curvature')
        for lam in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
            curv = add(scal(lam, add(dxY, dyX, -1)), scal(lam*lam, comm))
            if curv != scal(-lam*(1-lam), comm):
                raise ValueError('share law failed')
        split = add(dxY, dyX, -1)                       # F_1 + F_2 of the split protocols
        if split != scal(F(4), scal(F(-1, 4), comm)):
            raise ValueError('split-leg sum rule failed')
        rows.append(dict(point=[str(x), str(y)], commutator=[[str(v) for v in r] for r in comm]))
    # polygons: both pure readings return exactly
    pts = [(F(0), F(0)), (F(1, 2), F(0)), (F(3, 4), F(2, 3)), (F(-1, 5), F(1)), (F(-1, 2), F(1, 4))]
    Hs = [field(*p)[0] for p in pts]
    W = [[F(1), F(0)], [F(0), F(1)]]
    for a, b in zip(Hs, Hs[1:]+Hs[:1]):
        W = mm(mm(inv(b), a), W)
    if W != [[F(1), F(0)], [F(0), F(1)]]:
        raise ValueError('intensive reading has a return on a polygon')
    return dict(points=rows, share_values=['0', '1/4', '1/2', '3/4', '1'], polygon_vertices=len(pts),
                pure_readings_return_exactly=True)


# ------------------------------------------------------------- fluid part
def fluid(n, shrink=1.0):
    sys.path.insert(0, '../ph1')
    import ph1_real_fluid_holonomy as p
    H, _, _ = p.real_fluid('CarbonDioxide')
    Tm, dT, rm, dr = 330.0, 20.0*shrink, 450.0, 250.0*shrink
    T, rho = p.loop(n, Tm-dT, Tm+dT, rm-dr, rm+dr)
    return [H(float(t), float(r)) for t, r in zip(T, rho)], p


def frame_angle(W, H0):
    import numpy as np
    ev, O = np.linalg.eigh(H0)
    root = O@np.diag(np.sqrt(ev))@O.T
    rot = root@W@np.linalg.inv(root)
    return (math.atan2(rot[1, 0]-rot[0, 1], rot[0, 0]+rot[1, 1]),
            float(np.linalg.norm(rot@rot.T-np.eye(2))), rot)


def alternating(Hs):
    """Even segments: E reading (offsets frozen). Odd segments: M reading."""
    import numpy as np
    W = np.eye(2)
    n = len(Hs)
    for k in range(n):
        if k % 2:
            W = np.linalg.solve(Hs[(k+1) % n], Hs[k])@W
    return W


def shared(Hs, lam):
    import numpy as np
    W = np.eye(2)
    n = len(Hs)
    for k in range(n):
        h0, h1 = Hs[k], Hs[(k+1) % n]
        X = np.linalg.solve((h0+h1)/2, h1-h0)
        W = (np.eye(2)-lam*X+lam*lam*X@X/2)@W
    return W


def run():
    import numpy as np
    out = dict(exact=exact_controls(), fluid='CarbonDioxide', cycle='PH1 Table 1')
    Hs, p = fluid(16384)
    theta = p.holonomy(Hs)['theta_transport']
    out['theta_PH1'] = theta
    out['alternation'] = []
    for n in (8, 32, 128, 512, 2048, 16384):
        Hn, _ = fluid(n)
        ang, defect, _ = frame_angle(alternating(Hn), Hn[0])
        out['alternation'].append(dict(segments=n, angle=ang, orthogonality_defect=defect))
    out['share'] = []
    Hsmall, _ = fluid(16384, shrink=0.05)
    theta_small = frame_angle(shared(Hsmall, 0.5), Hsmall[0])[0]
    for lam in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        big = frame_angle(shared(Hs, lam), Hs[0])
        small = frame_angle(shared(Hsmall, lam), Hsmall[0])
        out['share'].append(dict(share=lam, law=4*lam*(1-lam), angle=big[0], angle_over_theta=big[0]/theta,
                                 orthogonality_defect=big[1], small_cycle_ratio=small[0]/theta_small))
    # a concrete probe: pure entropy offset of 1 J/(kg K) at the start state
    W = shared(Hs, 0.5)
    w0 = np.array([1.0, 0.0])
    w1 = W@w0
    e0, e1 = float(w0@Hs[0]@w0)/2, float(w1@Hs[0]@w1)/2
    out['probe'] = dict(start_offset=[1.0, 0.0], end_offset_ds_J_per_kgK=float(w1[0]),
                        end_offset_dv_m3_per_kg=float(w1[1]), second_order_energy_start=e0,
                        second_order_energy_end=e1)
    return out


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PH3_RESULT.json', 'w'), indent=1)
    print('exact controls pass; theta PH1 =', res['theta_PH1'])
    for r in res['alternation']:
        print('  alternation %6d segments: angle %+.6f  defect %.2e' % (r['segments'], r['angle'], r['orthogonality_defect']))
    for r in res['share']:
        print('  share %.2f  law %.4f  angle/Theta %.4f  defect %.2e  small-cycle ratio %.4f' % (
            r['share'], r['law'], r['angle_over_theta'], r['orthogonality_defect'], r['small_cycle_ratio']))
    print('  probe', res['probe'])
