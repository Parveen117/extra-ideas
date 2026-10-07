"""LN1: what a linear reading in a cut loses, the first non-linear layer of the tower returns.

The owner's statement: an observation keeps part of the information (spectral blindness); the lambda-tower moves
toward the whole without a new observation; the relation between the two is the bridge from linear to non-linear.
Sources: EMK-1 T2 (det = Delta_par + Delta_perp: observed and lost channels), NT-1 (H = ((a+c)/2) + ((a-c)/2) K + b RK),
RMG1 Cor 2.1 (spectral blindness), RMG2 (lambda-tower T_lambda(H) = H + lambda H^2; spectral operations never turn the
axes), RKF Theorem 42 sec. 7 (a higher layer repairs what a lower layer is blind to), QC4 (curvature of shared
readings = -variance), GB1 (scalar-only gravity bends light by half), NC1 (gravity element Exp(psi n.C)), PH3.
A reading in a declared cut sees the diagonal entries a_i of H; the couplings b_ij are lost to it.
N1  (H^2)_ii = a_i^2 + sum_j b_ij^2 :  non-linear layer in the same cut = (linear reading)^2 + exactly what was lost.
    lost_i = (H^2)_ii - (H_ii)^2  -- the variance of the reading.
N2  For every lambda != 0:  lost_i = ( T_lambda(H)_ii - a_i - lambda a_i^2 ) / lambda : every level carries it.
N3  Two cuts: the linear reading gives (a, c); one generation adds b^2; then det = a c - b^2 (EMK-1's two channels),
    the trace, the spectrum and the anisotropy l are all known.  The sign of b is never returned by any generation.
N4  Three cuts: the second layer returns all three b_ij^2; the third layer returns the sign of b12 b13 b23.
N5  Gravity element G = Exp(psi n.C) read in a cut across n: the linear reading is cosh(psi) in both slots
    (blind to the direction); lost = sinh^2(psi); and
        (G^2)_ii - 1 = [cosh^2(psi) - 1] + sinh^2(psi) = 2 sinh^2(psi) : the seen part and the lost part are EQUAL.
    GB1: scalar-only gravity bends light by half; the missing half is the lost part (0.876" + 0.876" = 1.751").
Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import itertools
import json


def mul(a, b):
    n = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def sym(diag, off):
    n = len(diag)
    H = [[F(0)]*n for _ in range(n)]
    for i in range(n):
        H[i][i] = F(diag[i])
    for (i, j), v in off.items():
        H[i][j] = H[j][i] = F(v)
    return H


def lost(H):
    H2 = mul(H, H)
    return [H2[i][i] - H[i][i]**2 for i in range(len(H))]


def run():
    # N1, N2
    H = sym([F(5), F(3)], {(0, 1): F(2)})
    if lost(H) != [F(4), F(4)]:
        raise ValueError('lost part must be b^2')
    H2 = mul(H, H)
    for lam in (F(1, 3), F(-2), F(7)):
        T = [[H[i][j] + lam*H2[i][j] for j in range(2)] for i in range(2)]
        for i in range(2):
            if (T[i][i] - H[i][i] - lam*H[i][i]**2)/lam != lost(H)[i]:
                raise ValueError('every level of the tower must carry the lost part')
    # N3
    two = []
    for a, c, b in ((F(5), F(3), F(2)), (F(7, 2), F(4), F(-3, 2)), (F(2), F(2), F(1))):
        Hp, Hm = sym([a, c], {(0, 1): b}), sym([a, c], {(0, 1): -b})
        b2 = lost(Hp)[0]
        det = a*c - b2
        if b2 != b*b or det != Hp[0][0]*Hp[1][1] - Hp[0][1]**2:
            raise ValueError('two channels failed')
        # no power of H read in the cut tells +b from -b
        Pp, Pm = Hp, Hm
        for _ in range(6):
            if [Pp[0][0], Pp[1][1]] != [Pm[0][0], Pm[1][1]]:
                raise ValueError('the sign of the coupling must stay hidden')
            Pp, Pm = mul(Pp, Hp), mul(Pm, Hm)
        two.append((str(a), str(c), str(b2), str(det)))
    # N4
    off = {(0, 1): F(2), (0, 2): F(-1), (1, 2): F(3)}
    H3 = sym([F(4), F(6), F(5)], off)
    l = lost(H3)                                   # l0 = b01^2 + b02^2, l1 = b01^2 + b12^2, l2 = b02^2 + b12^2
    b01 = (l[0] + l[1] - l[2])/2
    b02 = (l[0] + l[2] - l[1])/2
    b12 = (l[1] + l[2] - l[0])/2
    if (b01, b02, b12) != (F(4), F(1), F(9)):
        raise ValueError('three cuts: second layer must return all squares')
    cube = mul(mul(H3, H3), H3)
    flipped = sym([F(4), F(6), F(5)], {(0, 1): F(2), (0, 2): F(1), (1, 2): F(3)})     # one sign changed
    cube_f = mul(mul(flipped, flipped), flipped)
    same_second = lost(flipped) == l
    differ_third = [cube[i][i] for i in range(3)] != [cube_f[i][i] for i in range(3)]
    if not (same_second and differ_third):
        raise ValueError('third layer must carry the sign of the triple product')
    both = sym([F(4), F(6), F(5)], {(0, 1): F(-2), (0, 2): F(1), (1, 2): F(3)})       # two signs changed: same product
    cube_b = mul(mul(both, both), both)
    if [cube_b[i][i] for i in range(3)] != [cube[i][i] for i in range(3)]:
        raise ValueError('changing two signs must stay hidden at the third layer')
    # N5: gravity element, q = e^psi rational; cut across n
    grav = []
    for q in (F(2), F(3, 2), F(11, 10)):
        ch, sh = (q + 1/q)/2, (q - 1/q)/2
        G = [[ch, sh], [sh, ch]]                    # Exp(psi L) read in the K cut
        if G[0][0] != G[1][1] or lost(G)[0] != sh*sh:
            raise ValueError('gravity element across the cut failed')
        G2 = mul(G, G)
        seen, gone = ch*ch - 1, sh*sh
        if seen != gone or G2[0][0] - 1 != seen + gone:
            raise ValueError('seen and lost parts must be equal')
        grav.append((str(q), str(seen), str(gone), str(G2[0][0] - 1)))
    half = 0.8755
    return dict(two_cuts=two, three_cuts=dict(squares=[str(b01), str(b02), str(b12)], sign_at_layer=3),
                gravity=grav, light_at_sun_arcsec=dict(seen=half, lost=half, total=2*half))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('LN1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
