"""FR1: the frame is what cuts - one invariant, many readings; how many cuts a block carries.

Block M = a + b K + c R + d RK (EMK-1).  A cut is a traceless element C with C^2 = I.
T1  Turning the cut inside the split plane, K_phi = cos(phi) K + sin(phi) RK, changes the
    channels: D_par(phi) = a^2 - b_phi^2, D_perp(phi) = c^2 - d_phi^2, with the sum fixed.
    Cut-independent: a, c, b^2 + d^2.  Under all frame changes of a unit block only E = a survives.
T2  Cuts are the elements with form +1 (one sheet); turns have form -1 (two sheets); states are null.
T3  Number of mutually anticommuting cuts: 2 over the rationals (K, RK), 3 once the cut-complex
    unit iota is adjoined (K, RK, iota R), and no fourth.
T4  With three cuts used for propagation no linear element anticommutes with all of them: a single
    block has no mass term.  The conjugation Jc = (iota R) o conj anticommutes with all three and
    Jc^2 = -I:  (sum_i xi_i C_i + g Jc)^2 = (xi.xi - g^2) I.  Mass couples the block to its conjugate sheet.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../sy1')
import sy1_one_structure as sy1

I2, K, R, RK, mm, lin, block, parts, det = sy1.I2, sy1.K, sy1.R, sy1.RK, sy1.mm, sy1.lin, sy1.block, sy1.parts, sy1.det


# ------------------------------------------------------------ T1: turning the cut
def channels(M, cphi, sphi):
    """Channels of M read against the cut K_phi = cos K + sin RK (and its partner R K_phi)."""
    a, b, c, d = parts(M)
    b_phi = b*cphi+d*sphi
    d_phi = -b*sphi+d*cphi
    cut = lin((cphi, K), (sphi, RK))
    if mm(cut, cut) != I2:
        raise ValueError('not a cut')
    # independent route: components by scalar parts of products
    if parts(mm(M, cut))[0] != b_phi or parts(mm(M, mm(R, cut)))[0] != d_phi:
        raise ValueError('component along the turned cut failed')
    return a*a-b_phi*b_phi, c*c-d_phi*d_phi


def cut_control():
    angles = [(F(1), F(0)), (F(4, 5), F(3, 5)), (F(3, 5), F(4, 5)), (F(0), F(1)), (F(-5, 13), F(12, 13))]
    rows = []
    for name, M in (('boost along K', block(F(5, 4), F(3, 4), F(0), F(0))),
                    ('coin turn', block(F(3, 5), F(0), F(4, 5), F(0))),
                    ('general unit block', block(F(3, 5), F(4, 5), F(6, 5), F(2, 5))),
                    ('response', block(F(2), F(2, 3), F(0), F(1)))):
        a, b, c, d = parts(M)
        seen = []
        for cphi, sphi in angles:
            Dpar, Dperp = channels(M, cphi, sphi)
            if Dpar+Dperp != det(M):
                raise ValueError('the sum of the channels must not depend on the cut')
            if not a*a-(b*b+d*d) <= Dpar <= a*a:
                raise ValueError('channel left its range')
            seen.append((Dpar, Dperp))
        rows.append(dict(block=name, invariant_sum=str(det(M)), readings=[[str(x), str(y)] for x, y in seen],
                         cut_dependent=len(set(seen)) > 1))
    if rows[0]['readings'][0] != ['1', '0'] or rows[0]['readings'][3] != ['25/16', '-9/16']:
        raise ValueError('boost read along and across the cut failed')
    if rows[1]['cut_dependent']:
        raise ValueError('a turn must read the same against every turned cut')
    return rows


def frame_control():
    """Under similarity by unit blocks only the scalar part survives; c and b^2 + d^2 do not."""
    M = block(F(3, 5), F(4, 5), F(6, 5), F(2, 5))
    g = block(F(5, 4), F(0), F(0), F(3, 4))                    # a boost
    gi = block(F(5, 4), F(0), F(0), F(-3, 4))
    if mm(g, gi) != I2:
        raise ValueError('frame change is not invertible as declared')
    N = mm(mm(g, M), gi)
    a, b, c, d = parts(M)
    a2, b2, c2, d2 = parts(N)
    if a2 != a or det(N) != det(M) or b2*b2+d2*d2-c2*c2 != b*b+d*d-c*c:
        raise ValueError('scalar part or form changed')
    if c2 == c and b2*b2+d2*d2 == b*b+d*d:
        raise ValueError('a boost should change the turn component')
    return dict(scalar=str(a), turn_component_before=str(c), turn_component_after=str(c2))


def sheets_control():
    rows = []
    for name, X, want in (('cut K', K, 1), ('cut RK', RK, 1), ('turned cut', lin((F(3, 5), K), (F(4, 5), RK)), 1),
                          ('boosted cut', block(F(0), F(5, 4), F(3, 4), F(0)), 1),
                          ('turn R', R, -1), ('boosted turn', block(F(0), F(3, 4), F(5, 4), F(0)), -1),
                          ('state tensor (traceless part)', block(F(0), F(1), F(1), F(0)), 0)):
        q = -det(X)
        if mm(X, X) != lin((q, I2)) or q != want:
            raise ValueError('form of %s failed' % name)
        rows.append(dict(element=name, form=str(q)))
    return rows


# ------------------------------------------ T3 / T4: counting cuts, exactly
def kron(a, b):
    n, m = len(a), len(b)
    return [[a[i//m][j//m]*b[i % m][j % m] for j in range(n*m)] for i in range(n*m)]


def mmn(a, b):
    n = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def addn(a, b, s=1):
    n = len(a)
    return [[a[i][j]+s*b[i][j] for j in range(n)] for i in range(n)]


def scn(c, a):
    return [[c*x for x in row] for row in a]


def nullspace_dim(constraints, basis):
    """Dimension of {X in span(basis): every linear map in constraints sends X to 0}, exactly."""
    rows = []
    for con in constraints:
        images = [con(B) for B in basis]
        n = len(images[0])
        for i in range(n):
            for j in range(n):
                rows.append([img[i][j] for img in images])
    rank, width = 0, len(basis)
    rows = [r for r in rows if any(r)]
    for col in range(width):
        piv = next((k for k in range(rank, len(rows)) if rows[k][col]), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        p = rows[rank][col]
        rows[rank] = [x/p for x in rows[rank]]
        for k in range(len(rows)):
            if k != rank and rows[k][col]:
                f = rows[k][col]
                rows[k] = [x-f*y for x, y in zip(rows[k], rows[rank])]
        rank += 1
    return width-rank


def unit(n, i, j):
    return [[F(int((r, c) == (i, j))) for c in range(n)] for r in range(n)]


def counting_control():
    anti = lambda C, mul: (lambda X: [[x+y for x, y in zip(r1, r2)] for r1, r2 in zip(mul(C, X), mul(X, C))])
    # real 2x2: elements anticommuting with K and RK
    basis2 = [unit(2, i, j) for i in range(2) for j in range(2)]
    dim_real = nullspace_dim([anti(K, mm), anti(RK, mm)], basis2)
    if dim_real != 1 or add_zero(lin((1, mm(R, K)), (1, mm(K, R)))) or add_zero(lin((1, mm(R, RK)), (1, mm(RK, R)))):
        raise ValueError('over the rationals only R anticommutes with both cuts')
    if mm(R, R) == I2:
        raise ValueError('R must not be a cut')
    # cut-complex carrier as real 4x4: A (x) 1 for block elements, iota = 1 (x) j, conj = 1 (x) diag(1, -1)
    one, j = [[F(1), F(0)], [F(0), F(1)]], [[F(0), F(-1)], [F(1), F(0)]]
    C1, C2 = kron(K, one), kron(RK, one)
    iota = kron(I2, j)
    C3 = mmn(kron(R, one), iota)                               # iota R
    conj = kron(I2, [[F(1), F(0)], [F(0), F(-1)]])
    I4 = kron(I2, one)
    cuts = [C1, C2, C3]
    for a, A in enumerate(cuts):
        if mmn(A, A) != I4:
            raise ValueError('not a cut on the cut-complex carrier')
        for B in cuts[a+1:]:
            if addn(mmn(A, B), mmn(B, A)) != scn(F(0), I4):
                raise ValueError('cuts do not anticommute')
    basis4 = [unit(4, i, k) for i in range(4) for k in range(4)]
    commute_iota = lambda X: addn(mmn(iota, X), mmn(X, iota), -1)
    anticommute_iota = lambda X: addn(mmn(iota, X), mmn(X, iota))
    cons = [anti(C, mmn) for C in cuts]
    dim_linear = nullspace_dim(cons+[commute_iota], basis4)     # complex-linear solutions
    dim_all = nullspace_dim(cons, basis4)                       # all real-linear solutions
    dim_anti = nullspace_dim(cons+[anticommute_iota], basis4)   # antilinear solutions
    if (dim_linear, dim_all, dim_anti) != (0, 2, 2):
        raise ValueError('expected no linear fourth element and a two-dimensional conjugate family')
    Jc = mmn(C3, conj)
    for C in cuts:
        if addn(mmn(Jc, C), mmn(C, Jc)) != scn(F(0), I4):
            raise ValueError('conjugation does not anticommute with the cuts')
    if mmn(Jc, Jc) != scn(F(-1), I4) or addn(mmn(Jc, iota), mmn(iota, Jc)) != scn(F(0), I4):
        raise ValueError('conjugation must square to -1 and be antilinear')
    for xi, g in (((F(1), F(2), F(2)), F(3)), ((F(1, 2), F(-1, 3), F(5)), F(2, 7))):
        G = scn(g, Jc)
        for c, C in zip(xi, cuts):
            G = addn(G, scn(c, C))
        if mmn(G, G) != scn(sum(c*c for c in xi)-g*g, I4):
            raise ValueError('square law with the conjugate mass failed')
    return dict(real_block_cuts=2, real_block_third_element='R, a turn (square -1)',
                cut_complex_cuts=3, linear_fourth=dim_linear, conjugate_family=dim_anti,
                square_law='(sum xi_i C_i + g Jc)^2 = xi.xi - g^2')


def add_zero(M):
    return any(v for row in M for v in row)


def run():
    return dict(cuts=cut_control(), frame=frame_control(), sheets=sheets_control(), counting=counting_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('FR1_RESULT.json', 'w'), indent=1)
    for r in res['cuts']:
        print(r['block'], r['invariant_sum'], r['readings'][:4], r['cut_dependent'])
    print(res['frame']); print(res['sheets']); print(res['counting'])
