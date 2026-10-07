"""SY1: one structure - the unit EMK block, its scalar/traceless split, and every stage as a case.

Carrier: the certified EMK block (Publications EMK-1 T1, T2)
    M = a I + b K + c R + d RK,   K^2 = I, R^2 = -I, RK = -KR,
    det M = (a^2 - b^2) + (c^2 - d^2) = D_par + D_perp.
S1  E = a, O = M - a:  O^2 = (b^2 + d^2 - c^2) I;  on det M = 1:  E^2 - O^2 = I.
S2  three sectors by the sign of O^2 (circular, dual, split).
S3  defect det(I - M) = 2 (1 - a) on the unit quadric.
S4  powers: M^N = T_N(a) + U_{N-1}(a) O;  det(I - W^{2N}) = det(I - W^2) U_{N-1}(a_W)^2.
S5  composition: sc(M1 M2) = a1 a2 + B(O1, O2);  traceless part a1 O2 + a2 O1 + [O1, O2]/2.
S6  record {M, M^-1} with weights (p, 1-p): mean a + (2p-1) O, variance 4 p (1-p) O^2.
S7  state tensor psi psi^T = (n + j K + sigma RK)/2 is null: D_par = -D_perp = sigma^2/4.
S8  channels: coin  D_par = speed^2, D_perp = (curvature/2)^2;  response  -D_perp/D_par = r^2.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
from itertools import product
import json

I2 = [[F(1), F(0)], [F(0), F(1)]]
K = [[F(1), F(0)], [F(0), F(-1)]]
R = [[F(0), F(-1)], [F(1), F(0)]]


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def lin(*terms):
    out = [[F(0), F(0)], [F(0), F(0)]]
    for c, m in terms:
        for i in range(2):
            for j in range(2):
                out[i][j] += c*m[i][j]
    return out


RK = mm(R, K)


def block(a, b, c, d):
    return lin((a, I2), (b, K), (c, R), (d, RK))


def parts(M):
    a = (M[0][0]+M[1][1])/2
    b = (M[0][0]-M[1][1])/2
    d = (M[0][1]+M[1][0])/2
    c = (M[1][0]-M[0][1])/2
    if block(a, b, c, d) != M:
        raise ValueError('decomposition failed')
    return a, b, c, d


def det(M):
    return M[0][0]*M[1][1]-M[0][1]*M[1][0]


def cheb_T(n, x):
    a, b = F(1), x
    for _ in range(n):
        a, b = b, 2*x*b-a
    return a


def cheb_U(n, x):
    if n < 0:
        return F(0)
    a, b = F(1), 2*x
    for _ in range(n):
        a, b = b, 2*x*b-a
    return a


def relations_control():
    if mm(K, K) != I2 or mm(R, R) != lin((F(-1), I2)) or mm(R, K) != lin((F(-1), mm(K, R))):
        raise ValueError('EMK-1 T1 relations failed')
    if mm(RK, RK) != I2 or lin((1, mm(R, K)), (-1, mm(K, R))) != lin((F(2), RK)):
        raise ValueError('(RK)^2 = I or [R, K] = 2RK failed')
    checked = 0
    grid = [F(-2), F(-1, 2), F(0), F(1, 3), F(1), F(3, 2)]
    for a, b, c, d in product(grid, repeat=4):
        M = block(a, b, c, d)
        if det(M) != (a*a-b*b)+(c*c-d*d):
            raise ValueError('EMK-1 T2 determinant identity failed')
        O = lin((1, M), (-a, I2))
        if mm(O, O) != lin((b*b+d*d-c*c, I2)):
            raise ValueError('S1: O^2 is not the scalar b^2 + d^2 - c^2')
        checked += 1
    return checked


# rational unit blocks in the three sectors
UNITS = {
    'circular': [block(F(3, 5), F(0), F(4, 5), F(0)), block(F(5, 13), F(0), F(-12, 13), F(0)),
                 block(F(3, 5), F(4, 5), F(6, 5), F(2, 5)), block(F(4, 5), F(0), F(3, 5), F(0))],
    'split': [block(F(5, 4), F(3, 4), F(0), F(0)), block(F(5, 4), F(0), F(0), F(3, 4)),
              block(F(13, 12), F(3, 12), F(0), F(4, 12)), block(F(17, 8), F(9, 8), F(0), F(12, 8))],
    'dual': [block(F(1), F(1), F(1), F(0)), block(F(1), F(3, 5), F(1), F(4, 5)), block(F(1), F(0), F(0), F(0))],
}


def unit_control():
    rows = []
    for sector, blocks in UNITS.items():
        for M in blocks:
            a, b, c, d = parts(M)
            if det(M) != 1:
                raise ValueError('not on the unit quadric')
            O = lin((1, M), (-a, I2))
            q = b*b+d*d-c*c
            if lin((a*a, I2), (-1, mm(O, O))) != I2 or q != a*a-1:
                raise ValueError('S1: E^2 - O^2 = I failed')
            kind = 'circular' if q < 0 else 'split' if q > 0 else 'dual'
            if kind != sector:
                raise ValueError('S2: sector misread')
            if det(lin((1, I2), (-1, M))) != 2*(1-a):
                raise ValueError('S3: defect is not 2 (1 - a)')
            P = I2
            for N in range(1, 7):
                P = mm(P, M)
                chi = cheb_U(N-1, a)
                if P != lin((cheb_T(N, a), I2), (chi, O)):
                    raise ValueError('S4: power law failed')
                if kind == 'circular' and abs(chi) > N or kind == 'split' and a > 1 and chi < N or kind == 'dual' and a == 1 and chi != N:
                    raise ValueError('S4: character on the wrong side of N')
                W2N = mm(P, P)                                   # (M)^{2N} with W = M
                if det(lin((1, I2), (-1, W2N))) != det(lin((1, I2), (-1, mm(M, M))))*chi*chi:
                    raise ValueError('S4: defect character law failed')
            Minv = lin((a, I2), (-1, O))
            if mm(M, Minv) != I2:
                raise ValueError('inverse is not a - O')
            for p in (F(1), F(3, 4), F(1, 2), F(0)):
                mean = lin((p, M), (1-p, Minv))
                if mean != lin((a, I2), (2*p-1, O)):
                    raise ValueError('S6: record mean failed')
                dU, dV = lin((1, M), (-1, mean)), lin((1, Minv), (-1, mean))
                var = lin((p, mm(dU, dU)), (1-p, mm(dV, dV)))
                if var != lin((4*p*(1-p)*q, I2)):
                    raise ValueError('S6: record variance is not 4 p (1-p) O^2')
            rows.append(dict(sector=sector, scalar=str(a), O_squared=str(q), defect=str(2*(1-a))))
    return rows


def composition_control():
    flat = [M for blocks in UNITS.values() for M in blocks]
    checked = 0
    for M1 in flat:
        for M2 in flat:
            a1, b1, c1, d1 = parts(M1)
            a2, b2, c2, d2 = parts(M2)
            O1, O2 = lin((1, M1), (-a1, I2)), lin((1, M2), (-a2, I2))
            B = b1*b2+d1*d2-c1*c2
            P = mm(M1, M2)
            comm = lin((F(1, 2), mm(O1, O2)), (F(-1, 2), mm(O2, O1)))
            if P != lin((a1*a2+B, I2), (a1, O2), (a2, O1), (1, comm)):
                raise ValueError('S5: composition law failed')
            if parts(comm)[0] != 0:
                raise ValueError('S5: commutator part must be traceless')
            checked += 1
    return checked


def cases_control():
    """Each earlier stage's law read off the same block."""
    out = {}
    # PR1 / CL1: coin turn c + s R
    c, s = F(3, 5), F(4, 5)
    coin = block(c, F(0), s, F(0))
    a, b, cc, d = parts(coin)
    out['PR1 speed^2 + (curvature/2)^2 = 1'] = [str(a*a-b*b), str(cc*cc-d*d)]            # D_par, D_perp
    if (a*a-b*b, cc*cc-d*d) != (c*c, s*s):
        raise ValueError('coin channels failed')
    comm = lin((1, mm(K, coin)), (-1, mm(coin, K)))
    if comm != lin((-2*s, RK)):
        raise ValueError('curvature of the coin against the cut failed')
    out['CL1 mass^2 = clock curvature = det(I - coin)'] = str(det(lin((1, I2), (-1, coin))))
    # MS1 combination: sc(coin1 coin2) = c1 c2 - s1 s2
    c2, s2 = F(5, 13), F(12, 13)
    if parts(mm(coin, block(c2, F(0), s2, F(0))))[0] != c*c2-s*s2:
        raise ValueError('combination law failed')
    out['MS1 combined speed'] = str(c*c2-s*s2)
    # CL2 corner: two boosts, defect of the relative boost
    b1 = block(F(5, 4), F(3, 4), F(0), F(0))
    b2 = block(F(5, 4), F(-3, 4), F(0), F(0))
    rel = mm(b1, lin((F(5, 4), I2), (F(3, 4), K)))                                      # b1 b2^-1
    cosh = parts(rel)[0]
    if cosh != F(17, 8) or det(lin((1, I2), (-1, rel))) != -2*(cosh-1):
        raise ValueError('corner defect failed')
    out['CL2 corner cosh and split defect'] = [str(cosh), str(det(lin((1, I2), (-1, rel))))]
    # PH2 / PR1-T3: two non-collinear responses leave the R-part of the product
    r1, r2 = block(F(1), F(1, 3), F(0), F(1, 4)), block(F(1), F(-1, 5), F(0), F(1, 2))
    prod = parts(mm(r1, r2))
    cross = F(1, 3)*F(1, 2)-F(1, 4)*F(-1, 5)
    if abs(prod[2]) != abs(cross) or prod[0] != 1+F(1, 3)*F(-1, 5)+F(1, 4)*F(1, 2):
        raise ValueError('rotation left by two responses failed')
    out['PH2 rotation tangent'] = str(abs(prod[2]/prod[0]))
    # PR4 / MS1 state: psi psi^T is null, channels +-sigma^2/4
    psi = (F(3), F(1, 2))
    T = [[psi[i]*psi[j] for j in range(2)] for i in range(2)]
    a, b, cc, d = parts(T)
    n, j, sigma = 2*a, 2*b, 2*d
    if det(T) != 0 or cc != 0 or (a*a-b*b, cc*cc-d*d) != (sigma*sigma/4, -sigma*sigma/4):
        raise ValueError('state tensor is not null with channels +-sigma^2/4')
    out['MS1 state: speed, memory'] = [str(j/n), str(4*(a*a-b*b)/(n*n))]
    if 1-(j/n)**2 != 4*(a*a-b*b)/(n*n):
        raise ValueError('speed^2 + memory = 1 failed')
    # RMG9-T5: response m(1 + uK + v RK): coupling r^2 = -D_perp/D_par, ratio D_par/det
    u, v = F(1, 3), F(1, 2)
    Hm = block(F(2), 2*u, F(0), 2*v)
    a, b, cc, d = parts(Hm)
    Dpar, Dperp = a*a-b*b, cc*cc-d*d
    r2c = v*v/(1-u*u)
    if -Dperp/Dpar != r2c or Dpar/det(Hm) != 1/(1-r2c):
        raise ValueError('response channels failed')
    out['RMG9 coupling r^2 and ratio D_par/det'] = [str(r2c), str(Dpar/det(Hm))]
    # PR2: traceless X has X^2 = -(det X) I
    X = block(F(0), F(1, 2), F(1), F(1, 3))
    if mm(X, X) != lin((-det(X), I2)):
        raise ValueError('traceless square law failed')
    out['PR2 form of a traceless block'] = str(-det(X))
    # EMK-1 T7: commutes with K iff rotation-free
    for M, free in ((block(F(2), F(1), F(0), F(0)), True), (coin, False), (block(F(1), F(0), F(0), F(1, 2)), False)):
        if (mm(M, K) == mm(K, M)) != free:
            raise ValueError('EMK-1 T7 failed')
    return out


def run():
    return dict(grid_blocks=relations_control(), units=unit_control(), compositions=composition_control(),
                cases=cases_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('SY1_RESULT.json', 'w'), indent=1)
    print(res['grid_blocks'], res['compositions'])
    for r in res['units']:
        print(r)
    for k, v in res['cases'].items():
        print(k, v)
