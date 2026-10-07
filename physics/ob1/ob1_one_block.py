"""OB1: one block, one law - matter, light, potential, source and force as parts of the cut-complex block.

Block algebra over the cut-complex unit: scalar 1, three cuts C1 = K, C2 = RK, C3 = iota R,
three turns iota C_i, and iota = C1 C2 C3.  One operator D = d_t + sum C_i d_i and its reverse
Dbar = d_t - sum C_i d_i, with D Dbar = d_t^2 - Laplacian (units c = 1).
U1  C_i C_j = delta_ij + iota eps_ijk C_k;  C1 C2 C3 = iota;  D Dbar = wave operator.
U2  Potential (a reading-type block) Apot = phi - A.C.  Then  Dbar Apot = L + F.C  with
    L = d_t phi + div A (scalar) and F = E + iota B,  E = -grad phi - d_t A,  B = curl A.
U3  D (F.C) = (div E + iota div B) + (d_t E - curl B + iota (d_t B + curl E)).C :
    the four field equations are the four grades of one equation; for a field from a potential the
    two iota-parts vanish identically.
U4  The source J = D(F.C) is then a reading-type block (real scalar rho, real vector -j) and
    d_t rho + div j = 0 identically.
U5  Matter: for a column psi, d_t n + div r = psi^dagger (D psi) + (D psi)^dagger psi with
    (n; r) = (psi^dagger psi; psi^dagger C psi): a reading of D psi = 0 is a conserved reading-type block.
U6  Force: F rho + rho F^dagger has scalar part E.r and vector part n E + r x B, and leaves n^2 - r.r unchanged.
Exact arithmetic: polynomials in (t, x, y, z) over the Gaussian rationals.  stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json

# ---------------------------------------------------------------- polynomials
ZERO = {}


def pc(re, im=0):
    v = (F(re), F(im))
    return {(0, 0, 0, 0): v} if v != (0, 0) else {}


def var(i, power=1):
    e = [0, 0, 0, 0]
    e[i] = power
    return {tuple(e): (F(1), F(0))}


def padd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        old = out.get(k, (F(0), F(0)))
        new = (old[0]+s*v[0], old[1]+s*v[1])
        if new == (0, 0):
            out.pop(k, None)
        else:
            out[k] = new
    return out


def pmul(a, b):
    out = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(x+y for x, y in zip(k1, k2))
            old = out.get(k, (F(0), F(0)))
            out[k] = (old[0]+v1[0]*v2[0]-v1[1]*v2[1], old[1]+v1[0]*v2[1]+v1[1]*v2[0])
    return {k: v for k, v in out.items() if v != (0, 0)}


def pscale(z, a):
    return pmul(pc(*z) if isinstance(z, tuple) else pc(z), a)


def pconj(a):
    return {k: (v[0], -v[1]) for k, v in a.items()}


def pd(a, i):
    out = {}
    for k, v in a.items():
        if k[i]:
            e = list(k)
            n = e[i]
            e[i] -= 1
            out[tuple(e)] = (v[0]*n, v[1]*n)
    return out


IOTA = (0, 1)

# -------------------------------------------------------------------- blocks
def bconst(m):
    return [[pc(*m[i][j]) for j in range(2)] for i in range(2)]


C = [bconst([[(1, 0), (0, 0)], [(0, 0), (-1, 0)]]),      # K
     bconst([[(0, 0), (1, 0)], [(1, 0), (0, 0)]]),       # RK
     bconst([[(0, 0), (0, -1)], [(0, 1), (0, 0)]])]      # iota R
ONE = bconst([[(1, 0), (0, 0)], [(0, 0), (1, 0)]])


def bmul(a, b):
    return [[padd(pmul(a[i][0], b[0][j]), pmul(a[i][1], b[1][j])) for j in range(2)] for i in range(2)]


def badd(a, b, s=1):
    return [[padd(a[i][j], b[i][j], s) for j in range(2)] for i in range(2)]


def bscale(p, a):
    return [[pmul(p, a[i][j]) for j in range(2)] for i in range(2)]


def bdag(a):
    return [[pconj(a[j][i]) for j in range(2)] for i in range(2)]


def bd(a, i):
    return [[pd(a[r][s], i) for s in range(2)] for r in range(2)]


def parts(a):
    """scalar part and the three vector parts (complex polynomials)."""
    half = pc(F(1, 2))
    s = pmul(half, padd(a[0][0], a[1][1]))
    v = []
    for k in range(3):
        m = bmul(a, C[k])
        v.append(pmul(half, padd(m[0][0], m[1][1])))
    rebuilt = bscale(s, ONE)
    for k in range(3):
        rebuilt = badd(rebuilt, bscale(v[k], C[k]))
    if rebuilt != a:
        raise ValueError('block is not scalar + vector.C')
    return s, v


def vec_block(v):
    out = [[ZERO, ZERO], [ZERO, ZERO]]
    for k in range(3):
        out = badd(out, bscale(v[k], C[k]))
    return out


def D(a, sign=1):
    out = bd(a, 0)
    for k in range(3):
        out = badd(out, bmul(C[k], bd(a, k+1)), sign)
    return out


def grad(p):
    return [pd(p, i) for i in (1, 2, 3)]


def div(v):
    return padd(padd(pd(v[0], 1), pd(v[1], 2)), pd(v[2], 3))


def curl(v):
    return [padd(pd(v[2], 2), pd(v[1], 3), -1), padd(pd(v[0], 3), pd(v[2], 1), -1), padd(pd(v[1], 1), pd(v[0], 2), -1)]


def is_real(p):
    return all(v[1] == 0 for v in p.values())


# ------------------------------------------------------------------ controls
def algebra_control():
    eps = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (1, 0, 2): -1, (2, 1, 0): -1, (0, 2, 1): -1}
    for i in range(3):
        for j in range(3):
            want = bscale(pc(int(i == j)), ONE)
            for (a, b, k), s in eps.items():
                if (a, b) == (i, j):
                    want = badd(want, bscale(pc(0, s), C[k]))
            if bmul(C[i], C[j]) != want:
                raise ValueError('cut products failed')
    if bmul(bmul(C[0], C[1]), C[2]) != bscale(pc(0, 1), ONE):
        raise ValueError('iota is not the product of the three cuts')
    return dict(cuts=3, turns=3, volume='iota = C1 C2 C3')


def potentials():
    t, x, y, z = (var(i) for i in range(4))
    phi = padd(padd(pmul(pmul(t, x), y), pscale(F(2), pmul(z, z))), pscale(F(1, 3), pmul(pmul(x, x), x)))
    A = [padd(pmul(pmul(t, t), y), pscale(F(-1), pmul(x, z))),
         padd(pmul(pmul(x, y), z), pscale(F(1, 2), pmul(t, var(3, 2)))),
         padd(pscale(F(3), pmul(t, x)), pmul(var(2, 2), y))]
    return phi, A


def field_control():
    phi, A = potentials()
    Apot = badd(bscale(phi, ONE), vec_block(A), -1)
    if bdag(Apot) != Apot:
        raise ValueError('potential is not a reading-type block')
    s, v = parts(D(Apot, -1))
    L = padd(pd(phi, 0), div(A))
    E = [padd(pscale(F(-1), g), pd(a, 0), -1) for g, a in zip(grad(phi), A)]
    B = curl(A)
    Fv = [padd(e, pscale(IOTA, b)) for e, b in zip(E, B)]
    if s != L or v != Fv:
        raise ValueError('Dbar A is not L + (E + iota B).C')
    Fb = vec_block(Fv)
    s2, v2 = parts(D(Fb))
    divE, divB = div(E), div(B)
    if s2 != padd(divE, pscale(IOTA, divB)):
        raise ValueError('scalar grade of D F failed')
    want = [padd(padd(pd(e, 0), cb, -1), pscale(IOTA, padd(pd(b, 0), ce)))
            for e, b, cb, ce in zip(E, B, curl(B), curl(E))]
    if v2 != want:
        raise ValueError('vector grade of D F failed')
    if divB != ZERO or any(padd(pd(b, 0), ce) != ZERO for b, ce in zip(B, curl(E))):
        raise ValueError('the iota-parts must vanish for a field from a potential')
    # the source is reading-type and conserved
    rho, minus_j = s2, v2
    if not is_real(rho) or not all(is_real(p) for p in minus_j):
        raise ValueError('source is not a reading-type block')
    j = [pscale(F(-1), p) for p in minus_j]
    if padd(pd(rho, 0), div(j)) != ZERO:
        raise ValueError('source is not conserved')
    # D Dbar = wave operator, on the potential
    box = lambda p: padd(pd(pd(p, 0), 0), padd(padd(pd(pd(p, 1), 1), pd(pd(p, 2), 2)), pd(pd(p, 3), 3)), -1)
    lhs = D(D(Apot, -1))
    rhs = [[box(Apot[i][k]) for k in range(2)] for i in range(2)]
    if lhs != rhs:
        raise ValueError('D Dbar is not the wave operator')
    return dict(lorenz_scalar_terms=len(L), field_components=3, source_real=True, source_conserved=True,
                charge_density_terms=len(rho))


def wave_control():
    """A wave along x: E = (0, p, 0), B = (0, 0, p), p a polynomial of (x - t): D F = 0 and F.F = 0."""
    t, x = var(0), var(1)
    u = padd(x, t, -1)
    p = padd(pmul(pmul(u, u), u), pscale(F(2), u))
    Fv = [ZERO, p, pscale(IOTA, p)]
    if D(vec_block(Fv)) != [[ZERO, ZERO], [ZERO, ZERO]]:
        raise ValueError('wave does not solve D F = 0')
    ff = padd(padd(pmul(Fv[0], Fv[0]), pmul(Fv[1], Fv[1])), pmul(Fv[2], Fv[2]))
    if ff != ZERO:
        raise ValueError('wave is not null')
    sq = bmul(vec_block(Fv), vec_block(Fv))
    if sq != [[ZERO, ZERO], [ZERO, ZERO]]:
        raise ValueError('block square is not F.F')
    back = [ZERO, p, pscale((0, -1), p)]                          # wrong handedness for this direction
    if D(vec_block(back)) == [[ZERO, ZERO], [ZERO, ZERO]]:
        raise ValueError('the reversed sheet should not solve the same law')
    return dict(solves=True, null=True)


def matter_control():
    t, x, y, z = (var(i) for i in range(4))
    psi = [padd(pmul(t, x), pscale((1, 2), z)), padd(pscale((0, 1), pmul(y, y)), pscale(F(3), t))]
    col = [[psi[0], ZERO], [psi[1], ZERO]]                        # column as first column of a block
    rho = bmul(col, bdag(col))
    n2, r = parts(rho)                                            # rho = (n + r.C)/2: parts return n/2, r/2
    n, rvec = pscale(F(2), n2), [pscale(F(2), q) for q in r]
    if not is_real(n) or not all(is_real(q) for q in rvec):
        raise ValueError('matter reading is not reading-type')
    Dpsi = D(col)
    lhs = padd(pd(n, 0), div(rvec))
    m = badd(bmul(bdag(col), Dpsi), bmul(bdag(Dpsi), col))
    if lhs != m[0][0]:
        raise ValueError('reading conservation identity failed')
    null = padd(pmul(n, n), padd(padd(pmul(rvec[0], rvec[0]), pmul(rvec[1], rvec[1])), pmul(rvec[2], rvec[2])), -1)
    if null != ZERO:
        raise ValueError('a single reading must be null')
    return dict(identity='d_t n + div r = psi^dag (D psi) + (D psi)^dag psi', null=True)


def force_control():
    rows = []
    for (E, B, n, r) in (((1, 2, 0), (0, 1, 3), 5, (1, -2, 2)), ((F(1, 2), 0, -1), (2, 1, 1), 3, (0, 1, 0)),
                         ((0, 0, 0), (0, 0, 4), 2, (1, 1, 0)), ((3, 0, 0), (0, 0, 0), 4, (0, 0, 0))):
        Fb = vec_block([pc(e, b) for e, b in zip(E, B)])
        rho = badd(bscale(pc(F(n, 2)), ONE), bscale(pc(F(1, 2)), vec_block([pc(q) for q in r])))
        change = badd(bmul(Fb, rho), bmul(rho, bdag(Fb)))
        s, v = parts(change)
        Er = sum(F(a)*F(b) for a, b in zip(E, r))
        cross = [F(r[1])*F(B[2])-F(r[2])*F(B[1]), F(r[2])*F(B[0])-F(r[0])*F(B[2]), F(r[0])*F(B[1])-F(r[1])*F(B[0])]
        want_v = [n*F(E[k])+cross[k] for k in range(3)]
        if s != pc(Er) or v != [pc(q) for q in want_v]:
            raise ValueError('frame change by the field is not (E.r ; n E + r x B)')
        # with rho = (n + r.C)/2 the change of (n; r) is twice the parts: rate of n^2 - r.r
        if 2*n*(2*Er)-2*sum(F(r[k])*2*want_v[k] for k in range(3)) != 0:
            raise ValueError('the unrecoverable part changed')
        rows.append(dict(power=str(Er), force=[str(q) for q in want_v]))
    return rows


def run():
    return dict(algebra=algebra_control(), field=field_control(), wave=wave_control(), matter=matter_control(),
                force=force_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('OB1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
