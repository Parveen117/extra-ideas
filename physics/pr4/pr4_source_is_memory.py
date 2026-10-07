"""PR4: what a sector can offer the clock - flip rate times proper density, and proper density is memory.

Sector law (c = 1): d_t psi = (A d_x + g R) psi, A = H = diag(1, -1), R = K H, K the exchange.
  n = psi.psi          density            j = psi.A psi        current
  sigma = psi.K psi    scalar             tau = g sigma        the source candidate
T1  d_t n = d_x j   (conserved)           d_t j - d_x n = -2 tau   (dual residue)
T2  (n, j) is a vector under boosts, sigma a scalar, sigma^2 = n^2 - j^2
T3  with p the share of the first light-like reading: j = (2p-1) n, sigma^2 = n^2 * 4p(1-p)
T4  tau is even under the sheet map (psi -> A psi, g -> -g); every boost-invariant bilinear is a multiple of sigma
T5  in the accelerated frame of PR3 the same two identities hold in local units
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../pr3')
import pr3_accelerated_frame as pr3

fadd, scale, du, dv, d_s, rho, clean = pr3.fadd, pr3.scale, pr3.du, pr3.dv, pr3.d_s, pr3.rho, pr3.clean


def fmul(f, g):
    out = {}
    for (a, b), c in f.items():
        for (d, e), h in g.items():
            out[(a+d, b+e)] = out.get((a+d, b+e), F(0))+c*h
    return clean(out)


def d_t(f):
    return fadd(du(f), dv(f), -1)


def d_x(f):
    return fadd(du(f), dv(f))


def bilinears(psi):
    p1, p2 = psi
    n = fadd(fmul(p1, p1), fmul(p2, p2))
    j = fadd(fmul(p1, p1), fmul(p2, p2), -1)
    sigma = scale(F(2), fmul(p1, p2))
    return n, j, sigma


def rhs(psi, g):
    """(A d_x + g R) psi,  R psi = (-psi2, psi1)."""
    p1, p2 = psi
    return (fadd(d_x(p1), scale(g, p2), -1), fadd(scale(F(-1), d_x(p2)), scale(g, p1)))


def dot2(a, b):
    return fadd(fmul(a[0], b[0]), fmul(a[1], b[1]))


def current_control(g):
    for psi in pr3.TESTS:
        n, j, sigma = bilinears(psi)
        dpsi = rhs(psi, g)                                    # d_t psi on solutions
        dn = scale(F(2), dot2(psi, dpsi))
        dj = scale(F(2), dot2((psi[0], scale(F(-1), psi[1])), dpsi))
        if dn != d_x(j):
            raise ValueError('density is not conserved')
        if fadd(dj, d_x(n), -1) != scale(-2*g, sigma):
            raise ValueError('dual residue is not -2 tau')
        if g and not sigma:
            raise ValueError('witness has no source')
    # at g = 0 both currents are conserved: the flat (light-like) sector
    for psi in pr3.TESTS:
        n, j, _ = bilinears(psi)
        dpsi = rhs(psi, F(0))
        if scale(F(2), dot2((psi[0], scale(F(-1), psi[1])), dpsi)) != d_x(n):
            raise ValueError('massless sector should have no residue')
    return dict(conserved='d_t n = d_x j', residue='d_t j - d_x n = -2 g sigma')


def frame_control(g):
    """PR3 frame: d_eta chi = (A d_s + g rho R) chi; local units d_rho = (1/rho) d_s."""
    for chi in pr3.TESTS:
        n, j, sigma = bilinears(chi)
        dchi = pr3.G(chi, g)
        dn = scale(F(2), dot2(chi, dchi))
        dj = scale(F(2), dot2((chi[0], scale(F(-1), chi[1])), dchi))
        if dn != d_s(j):
            raise ValueError('frame conservation failed')
        # (1/rho) d_eta j - d_rho n = -2 g sigma   <=>   d_eta j - d_s n = -2 g rho sigma
        if fadd(dj, d_s(n), -1) != rho(sigma, 1, -2*g):
            raise ValueError('frame residue is not -2 tau in local units')
    return dict(residue_in_frame='(1/rho) d_eta j - d_rho n = -2 g sigma')


def invariance_control():
    rows = []
    for psi, b in (((F(3), F(1, 2)), F(2)), ((F(-1, 3), F(5)), F(3)), ((F(2), F(0)), F(1, 2))):
        n, j, sigma = psi[0]**2+psi[1]**2, psi[0]**2-psi[1]**2, 2*psi[0]*psi[1]
        if sigma*sigma != n*n-j*j:
            raise ValueError('sigma^2 = n^2 - j^2 failed')
        q = (psi[0]*b, psi[1]/b)                              # boost Exp(eta A / 2), e^{eta/2} = b
        ch, sh = (b*b+1/(b*b))/2, (b*b-1/(b*b))/2
        n2, j2, s2 = q[0]**2+q[1]**2, q[0]**2-q[1]**2, 2*q[0]*q[1]
        if (n2, j2) != (n*ch+j*sh, j*ch+n*sh) or s2 != sigma:
            raise ValueError('(n, j) is not a vector or sigma not a scalar')
        p = psi[0]**2/n
        if j != (2*p-1)*n or sigma*sigma != n*n*4*p*(1-p):
            raise ValueError('proper density is not density times root of the share memory')
        rows.append(dict(share=str(p), memory=str(4*p*(1-p)), sigma_over_n_squared=str(sigma*sigma/(n*n))))
    if rows[-1]['memory'] != '0':
        raise ValueError('a pure reading must have no memory and no source')
    return rows


def uniqueness_control():
    """B^T M B = M for B = diag(b, 1/b), b = 2 and 3, forces M11 = M22 = 0: bilinear = (M12 + M21) psi1 psi2."""
    free = []
    for i in range(2):
        for k in range(2):
            ok = True
            for b in (F(2), F(3)):
                w = (b, 1/b)
                if w[i]*w[k] != 1:
                    ok = False
            if ok:
                free.append((i, k))
    if free != [(0, 1), (1, 0)]:
        raise ValueError('invariant bilinears are not spanned by the exchange')
    return dict(invariant_entries=['12', '21'], bilinear='multiple of sigma')


def sheet_control(g):
    for psi in pr3.TESTS:
        n, j, sigma = bilinears(psi)
        flipped = (psi[0], scale(F(-1), psi[1]))              # A psi
        n2, j2, s2 = bilinears(flipped)
        if (n2, j2) != (n, j) or s2 != scale(F(-1), sigma):
            raise ValueError('sheet map parities failed')
        if scale(-g, s2) != scale(g, sigma):
            raise ValueError('tau is not sheet-even')
        # and the flipped field solves the law with -g
        a = rhs(flipped, -g)
        b = rhs(psi, g)
        if a != (b[0], scale(F(-1), b[1])):
            raise ValueError('sheet map does not carry the law to its mirror')
    return dict(n='even', j='even', sigma='odd', tau='even')


def run(g=F(3, 2)):
    return dict(currents=current_control(g), frame=frame_control(g), invariance=invariance_control(),
                uniqueness=uniqueness_control(), sheet=sheet_control(g))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PR4_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
