"""QC1: closure of the return under a reading of content q.

Exact part (stdlib): the two readings as a canonical pair, the return's
curvature as their bracket, the content spectrum of the return, the silent
circles and their weight ladder.  Fluid part (numpy + CoolProp): CO2 cycles
approaching the critical point and the first covariance-silent cycle.
Python 3.12.
"""
from fractions import Fraction as F
import json
import math
import sys

sys.path.insert(0, '../ph3')
import ph3_two_readings_transport as ph3

mm, add, scal, inv = ph3.mm, ph3.add, ph3.scal, ph3.inv


def vec(a, w):
    return [a[0][0]*w[0]+a[0][1]*w[1], a[1][0]*w[0]+a[1][1]*w[1]]


# ---------------------------------------------------------------- T1, T2
def pair_controls():
    rows = []
    for x, y in ((F(1, 2), F(1, 3)), (F(-2, 5), F(3, 4)), (F(1), F(-1, 7))):
        H, Hx, Hy, _, _ = ph3.field(x, y)
        Hi = inv(H)
        w = [F(3, 7), F(-2, 5)]
        p = vec(H, w)
        for A in (Hx, Hy):
            B = scal(F(-1), mm(mm(Hi, A), Hi))            # derivative of H^-1
            # E reading: w fixed, p moves by A w.  M reading: p fixed, w moves by B p.
            dw = scal(F(1, 2), [vec(B, p)])[0]
            dp = scal(F(1, 2), [vec(A, w)])[0]
            native = scal(F(-1, 2), [vec(mm(Hi, A), w)])[0]
            if dw != native:
                raise ValueError('equal share is not the native transport on the graph')
            graph = [a+b for a, b in zip(vec(A, w), vec(H, dw))]   # d(Hw) along the motion
            if graph != dp:
                raise ValueError('equal share leaves the equilibrium graph')
            if sum(a*b for a, b in zip(dp, w))+sum(a*b for a, b in zip(p, dw)) != 0:
                raise ValueError('pairing p.w is not conserved')
        Bx = scal(F(-1), mm(mm(Hi, Hx), Hi))
        By = scal(F(-1), mm(mm(Hi, Hy), Hi))
        bracket = add(mm(Hx, By), mm(Hy, Bx), -1)          # A_s B_v - A_v B_s
        X, Y = mm(Hi, Hx), mm(Hi, Hy)
        native_curv = scal(F(-1, 4), add(mm(X, Y), mm(Y, X), -1))
        if scal(F(1, 4), mm(mm(Hi, bracket), H)) != native_curv:
            raise ValueError('curvature is not the bracket of the two readings')
        if not any(v for row in bracket for v in row):
            raise ValueError('witness readings commute')
        rows.append(dict(point=[str(x), str(y)], bracket=[[str(v) for v in r] for r in bracket]))
    return rows


# -------------------------------------------------- T3: content spectrum
def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def poly_mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            key = (i+k, j+l)
            val = cmul(a, b)
            old = out.get(key, (F(0), F(0)))
            out[key] = (old[0]+val[0], old[1]+val[1])
    return {k: v for k, v in out.items() if v != (F(0), F(0))}


def poly_pow(p, n):
    out = {(0, 0): (F(1), F(0))}
    for _ in range(n):
        out = poly_mul(out, p)
    return out


def content_control(c, s, degree):
    """z^a zbar^b under the rotation (c, s) is multiplied by (c + i s)^(a-b)."""
    z = {(1, 0): (F(1), F(0)), (0, 1): (F(0), F(1))}
    zb = {(1, 0): (F(1), F(0)), (0, 1): (F(0), F(-1))}
    y1 = {(1, 0): (c, F(0)), (0, 1): (-s, F(0))}
    y2 = {(1, 0): (s, F(0)), (0, 1): (c, F(0))}
    zr = {k: v for k, v in poly_mul({(0, 0): (F(1), F(0))}, y1).items()}
    for k, v in y2.items():
        old = zr.get(k, (F(0), F(0)))
        zr[k] = (old[0]-v[1], old[1]+v[0])                 # + i y2
    zbr = {k: (v[0], -v[1]) for k, v in zr.items()}
    contents = []
    for a in range(degree+1):
        b = degree-a
        left = poly_mul(poly_pow(zr, a), poly_pow(zbr, b))
        ph = (F(1), F(0))
        for _ in range(abs(a-b)):
            ph = cmul(ph, (c, s if a >= b else -s))
        right = poly_mul({(0, 0): ph}, poly_mul(poly_pow(z, a), poly_pow(zb, b)))
        if left != right:
            raise ValueError('content law failed')
        contents.append(a-b)
    return sorted(contents)


# ---------------------------------------------------- T4: silent circles
def ladder(q, count):
    """Circles about isotropy: Theta = -pi (cosh d - 1). Content q is silent iff q Theta in 2 pi Z."""
    k = F(q, 2)
    rows = []
    for n in range(1, count+1):
        cosh = 1+F(2*n, q)
        turns = -F(q)*(cosh-1)/2                           # q Theta / (2 pi)
        if turns != -n or k*cosh != k+n:
            raise ValueError('weight ladder failed')
        rows.append(dict(n=n, cosh_d=str(cosh), height=str(k*cosh), klein_radius_squared=str(1-1/cosh**2)))
    return dict(content=q, index_k=str(k), circles=rows)


def exact():
    return dict(pair=pair_controls(),
                contents={str(m): content_control(F(3, 5), F(4, 5), m) for m in (1, 2, 3, 4)},
                ladders=[ladder(q, 4) for q in (1, 2, 3, 4)])


# ---------------------------------------------------------- fluid part
def fluid():
    import numpy as np
    sys.path.insert(0, '../ph1')
    import ph1_real_fluid_holonomy as p
    H, crit, _ = p.real_fluid('CarbonDioxide')
    Tc = crit['Tc']

    def cycle(delta, n=6000):
        T, rho = p.loop(n, Tc+delta, 350.0, 200.0, 700.0)
        return [H(float(t), float(r)) for t, r in zip(T, rho)]

    table = []
    for delta in (5.87, 3.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.05, 0.03, 0.02, 0.01, 0.001):
        table.append(dict(kelvin_above_critical=delta, theta=p.holonomy(cycle(delta))['theta_loop']))
    lo, hi = 0.015, 0.04                                    # theta(lo) < -pi < theta(hi)
    f = lambda d: p.holonomy(cycle(d))['theta_loop']+math.pi
    flo, fhi = f(lo), f(hi)
    silent = None
    if flo < 0 < fhi:
        for _ in range(40):
            mid = math.sqrt(lo*hi)
            if f(mid) < 0:
                lo = mid
            else:
                hi = mid
        silent = math.sqrt(lo*hi)
        Hs = cycle(silent, 24000)
        W = ph3.shared(Hs, 0.5)
        C = np.array([[2.0, 0.3], [0.3, 1.0]])
        C = np.linalg.inv(Hs[0])@C@np.linalg.inv(Hs[0])*1e-3
        moved = W@C@W.T
        w = np.array([1.0, 0.0])
        silent = dict(kelvin_above_critical=silent, theta=p.holonomy(Hs)['theta_loop'],
                      covariance_relative_change=float(np.linalg.norm(moved-C)/np.linalg.norm(C)),
                      first_moment_relative_change=float(np.linalg.norm(W@w-w)/np.linalg.norm(w)))
    return dict(fluid='CarbonDioxide', family='PH1 ellipse with lowest temperature T_c + delta',
                table=table, covariance_silent_cycle=silent)


if __name__ == '__main__':
    out = dict(exact=exact(), fluid=fluid())
    json.dump(out, open('QC1_RESULT.json', 'w'), indent=1)
    print('contents', out['exact']['contents'])
    for lad in out['exact']['ladders']:
        print('q=%d k=%s' % (lad['content'], lad['index_k']), [(c['cosh_d'], c['height']) for c in lad['circles']])
    for r in out['fluid']['table']:
        print('  delta %.3f K  theta %+.4f' % (r['kelvin_above_critical'], r['theta']))
    print('  silent', out['fluid']['covariance_silent_cycle'])
