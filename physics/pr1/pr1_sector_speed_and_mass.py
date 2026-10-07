"""PR1: sector speed, mass and curvature from one coin turn; the relativity of a sector.

A  Native walk (R28 source roles H, K, R = KH; address shift S; D = S P + S^-1 Q):
   U = D Exp(theta R):  U + U^-1 = cos(theta) (S + S^-1);  cone speed cos(theta);
   rest turn theta;  [H, Exp(theta R)] = -2 sin(theta) K;  speed^2 + (curvature/2)^2 = 1.
B  Response space is velocity space: collinear responses add by the relativistic law;
   non-collinear ones leave the rotation of PH2.
C  Continuum sector law d_t = c A d_x + g R: square law, covariance under its own boost.
D  The flip between the two light-like readings in its two parities (QC3):
   turn g R -> mass;  record g (S - 1) -> telegraph law, diffusion c^2 / (2 g).
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../ph2')
import ph2_strain_cycle_rotation as ph2

ONE = [[F(1), F(0)], [F(0), F(1)]]
H = [[F(1), F(0)], [F(0), F(-1)]]          # cut role (EMK K-type generator)
K = [[F(0), F(1)], [F(1), F(0)]]           # exchange role (EMK S-type generator)
R = [[F(0), F(-1)], [F(1), F(0)]]          # R = K H


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def add(a, b, s=1):
    return [[a[i][j]+s*b[i][j] for j in range(2)] for i in range(2)]


def sc(c, a):
    return [[c*x for x in row] for row in a]


def roles_control():
    if mm(K, H) != R or mm(R, R) != sc(F(-1), ONE) or add(mm(H, K), mm(K, H)) != sc(F(0), ONE):
        raise ValueError('source role relations failed')
    if add(mm(H, K), mm(K, H), -1) != sc(F(-2), R):
        raise ValueError('R28 curvature [H, K] = -2R failed')
    return True


# ------------------------------------------------ Laurent matrices in the shift
def lmul(a, b):
    out = [[{}, {}], [{}, {}]]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for p, x in a[i][k].items():
                    for q, y in b[k][j].items():
                        out[i][j][p+q] = out[i][j].get(p+q, F(0))+x*y
            out[i][j] = {p: v for p, v in out[i][j].items() if v}
    return out


def const(m):
    return [[({0: m[i][j]} if m[i][j] else {}) for j in range(2)] for i in range(2)]


def ladd(a, b, s=1):
    out = [[{}, {}], [{}, {}]]
    for i in range(2):
        for j in range(2):
            d = dict(a[i][j])
            for p, v in b[i][j].items():
                d[p] = d.get(p, F(0))+s*v
            out[i][j] = {p: v for p, v in d.items() if v}
    return out


def scalar(poly):
    clean = {p: v for p, v in poly.items() if v}
    return [[dict(clean), {}], [{}, dict(clean)]]


SHIFT = [[{1: F(1)}, {}], [{}, {-1: F(1)}]]            # D = S P + S^-1 Q
UNSHIFT = [[{-1: F(1)}, {}], [{}, {1: F(1)}]]


def walk_control(c, s):
    if c*c+s*s != 1:
        raise ValueError('not a turn')
    coin = add(sc(c, ONE), sc(s, R))                     # Exp(theta R)
    coin_inv = add(sc(c, ONE), sc(-s, R))
    U = lmul(SHIFT, const(coin))
    Uinv = lmul(const(coin_inv), UNSHIFT)
    if lmul(U, Uinv) != const(ONE):
        raise ValueError('walk is not invertible as declared')
    if ladd(U, Uinv) != scalar({1: c, -1: c}):
        raise ValueError('U + U^-1 = cos(theta)(S + S^-1) failed')
    curvature = add(mm(H, coin), mm(coin, H), -1)
    if curvature != sc(-2*s, K):
        raise ValueError('[H, Exp(theta R)] = -2 sin(theta) K failed')
    half = s                                             # |curvature| / 2
    if c*c+half*half != 1:
        raise ValueError('speed^2 + (curvature/2)^2 = 1 failed')
    # dispersion cos w = c cos k on rational address turns; group speed below the cone speed
    rows = []
    for t in (F(1, 7), F(1, 3), F(1, 2), F(1), F(3, 2)):
        ck, sk = (1-t*t)/(1+t*t), 2*t/(1+t*t)
        Ew = c*ck                                        # even part of the event turn
        sin2w = 1-Ew*Ew
        if sin2w != s*s+c*c*sk*sk:
            raise ValueError('event odd part identity failed')
        if (1-Ew) != (1-c)+c*(1-ck):
            raise ValueError('defect composition failed')
        v2 = c*c*sk*sk/sin2w                             # (d w / d k)^2
        if v2 > c*c or (v2 == c*c) != (sk == 1 or s == 0 or c == 0):
            raise ValueError('group speed exceeds the cone speed')
        rows.append(dict(cos_k=str(ck), cos_w=str(Ew), group_speed_squared=str(v2)))
    rest = add(sc(c, ONE), sc(s, R))                     # at k = 0 the walk is the coin turn itself
    return dict(cos_theta=str(c), sin_theta=str(s), cone_speed=str(c), curvature_over_two=str(s),
                rest_turn_even_part=str(rest[0][0]), modes=rows)


def hadamard_control():
    """R28's coin (H + K)/sqrt2 without the root: Ut = D (H + K), Ut^2 - (S - S^-1) Ut - 2 = 0."""
    Ut = lmul(SHIFT, const(add(H, K)))
    lhs = ladd(ladd(lmul(Ut, Ut), lmul(scalar({1: F(1), -1: F(-1)}), Ut), -1), const(sc(F(2), ONE)), -1)
    if lhs != const(sc(F(0), ONE)):
        raise ValueError('R28 walk identity failed')
    # V = Ut^2 / 2:  V + V^-1 = 2 - (1/2)(2 - S^2 - S^-2), so the cone speed squared is 1/2 (R38)
    V2 = lmul(Ut, Ut)                                    # = 2 V
    V4 = lmul(V2, V2)                                    # = 4 V^2
    want = ladd(lmul(scalar({0: F(2), 2: F(1), -2: F(1)}), V2), const(sc(F(4), ONE)), -1)
    # 4V^2 + 4 = 2 (2 - (1/2)(2 - S^2 - S^-2)) * 2V  <=>  4V^2 = (2 + S^2 + S^-2)(2V) - 4
    if V4 != want:
        raise ValueError('two-event wave identity failed')
    return dict(cone_speed_squared='1/2', curvature_over_two_squared='1/2')


# --------------------------------------------- B: response space = velocity space
def response(u, v):
    return add(ONE, add(sc(u, H), sc(v, K)))


def velocity_control():
    u1, u2 = F(1, 3), F(2, 5)
    P = mm(response(u1, 0), response(u2, 0))
    ratio = (P[0][0]-P[1][1])/(P[0][0]+P[1][1])
    if ratio != (u1+u2)/(1+u1*u2) or ph2.polar_tangent(P)[0] != 0:
        raise ValueError('collinear responses do not add by the relativistic law')
    w1, w2 = (F(1, 3), F(1, 4)), (F(-1, 5), F(1, 2))
    P = mm(response(*w1), response(*w2))
    num, den = ph2.polar_tangent(P)
    cross = w1[0]*w2[1]-w1[1]*w2[0]
    dotp = w1[0]*w2[0]+w1[1]*w2[1]
    if F(num, 1)/den not in (cross/(1+dotp), -cross/(1+dotp)) or num == 0:
        raise ValueError('rotation of two responses is not cross/(1 + dot)')
    PtP = mm([list(r) for r in zip(*P)], P)
    detP = P[0][0]*P[1][1]-P[0][1]*P[1][0]
    cosh3 = (PtP[0][0]+PtP[1][1])/(2*detP)
    r1, r2 = w1[0]**2+w1[1]**2, w2[0]**2+w2[1]**2
    if cosh3 != ((1+r1)*(1+r2)+4*dotp)/((1-r1)*(1-r2)):
        raise ValueError('composition of rapidities failed')
    return dict(collinear_sum=str(ratio), rotation_tangent=str(abs(F(num, 1)/den)), composed_cosh=str(cosh3))


# ------------------------------------------------- C: continuum sector law
def sector_control(c, g, ch, sh):
    """Operator d_t - c A d_x - g R with A = H; boost B = Exp(eta A / 2), (ch, sh) = cosh, sinh of eta/2."""
    if ch*ch-sh*sh != 1:
        raise ValueError('not a boost')
    A = H
    if mm(A, A) != ONE or add(mm(A, R), mm(R, A)) != sc(F(0), ONE):
        raise ValueError('propagation and mass generators must anticommute')
    square = mm(add(sc(c, A), sc(g, R)), add(sc(c, A), sc(g, R)))    # symbols xi = 1 for d_x
    if square != sc(c*c-g*g, ONE):
        raise ValueError('square law failed')
    B = add(sc(ch, ONE), sc(sh, A))
    C, Sh = ch*ch+sh*sh, 2*ch*sh                                     # cosh eta, sinh eta
    if mm(B, B) != add(sc(C, ONE), sc(Sh, A)) or mm(mm(B, A), B) != add(sc(C, A), sc(Sh, ONE)):
        raise ValueError('boost does not act on the time and space coefficients as a boost')
    if mm(mm(B, R), B) != R:
        raise ValueError('mass term is not invariant')
    # B (d_t - c A d_x - g R) B = (C d_t - c Sh d_x) - c A (C d_x - (Sh/c) d_t) - g R:
    # the boost built on the sector's own speed c keeps d_t^2 - c^2 d_x^2.  Substituting
    # d_t = C d_t' - c Sh d_x', d_x = C d_x' - (Sh/c) d_t' into d_t^2 - w^2 d_x^2:
    def boosted(w):
        tt = C*C-w*w*Sh*Sh/(c*c)
        xx = c*c*Sh*Sh-w*w*C*C
        tx = -2*c*C*Sh+2*w*w*C*Sh/c
        return tt, xx, tx
    if boosted(c) != (1, -c*c, 0):
        raise ValueError('sector is not invariant under its own boost')
    if boosted(c/2)[2] == 0 or boosted(2*c)[2] == 0:
        raise ValueError('a sector of another speed must not be invariant under this boost')
    return dict(speed=str(c), mass_rate=str(g), cosh_eta=str(C), sinh_eta=str(Sh),
                second_order='d_t^2 = c^2 d_x^2 - g^2')


# -------------------------------------------------- D: two parities of the flip
def parity_control(c, g, k):
    turn = add(sc(c*k, H), sc(g, R))                     # symbol of c H d_x + g R on a mode, d_x^2 = -k^2
    # (c H d_x + g R)^2 = c^2 d_x^2 - g^2  ->  frequency^2 = c^2 k^2 + g^2
    if mm(add(sc(c, H), sc(g, R)), add(sc(c, H), sc(g, R))) != sc(c*c-g*g, ONE):
        raise ValueError('turn parity square law failed')
    freq2 = c*c*k*k+g*g
    record = add(sc(c, H), sc(g, K))                     # c H d_x + g K, the record generator plus g
    if mm(record, record) != sc(c*c+g*g, ONE):
        raise ValueError('record parity square law failed')
    # (d_t + g)^2 = c^2 d_x^2 + g^2  ->  rates l: l^2 + 2 g l + c^2 k^2 = 0
    disc = g*g-c*c*k*k
    lo, hi = c*c*k*k/(2*g), c*c*k*k/g                    # diffusion rate D k^2 and twice it
    if disc >= 0 and g >= hi:
        # the slow rate g - sqrt(disc) lies in [D k^2, 2 D k^2] with D = c^2/(2g)
        if not (g-hi)**2 <= disc <= (g-lo)**2:
            raise ValueError('slow rate is not the diffusion rate to leading order')
    return dict(speed=str(c), flip_rate=str(g), wave_number=str(k), frequency_squared=str(freq2),
                diffusion_constant=str(c*c/(2*g)), diffusion_times_two_flip_rate=str(c*c))


def run():
    roles_control()
    return dict(walk=[walk_control(c, s) for c, s in ((F(1), F(0)), (F(4, 5), F(3, 5)), (F(3, 5), F(4, 5)),
                                                      (F(5, 13), F(12, 13)), (F(0), F(1)))],
                hadamard=hadamard_control(), velocity=velocity_control(),
                sector=[sector_control(F(1), F(3, 4), F(5, 4), F(3, 4)), sector_control(F(1, 3), F(2), F(13, 12), F(5, 12))],
                parity=[parity_control(F(1), F(4), F(3)), parity_control(F(1), F(5), F(3)), parity_control(F(1, 2), F(10), F(1))])


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PR1_RESULT.json', 'w'), indent=1)
    for w in res['walk']:
        print({k: v for k, v in w.items() if k != 'modes'}, w['modes'][2])
    print(res['hadamard']); print(res['velocity'])
    for r in res['sector']+res['parity']:
        print(r)
