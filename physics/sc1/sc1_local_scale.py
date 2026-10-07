"""SC1: the scale factor made local - what it gives, and why it is not a field.

GB1 left the scale factor of the block unused.  Make it local, as the phase is made local for light.
S1  Law with a self-dagger potential W = w0 + w.C:  (D + W) psi = 0.  Then
        d_t n + div r = psi^dag((D+W)psi) + ((D+W)psi)^dag psi - 2 (w0 n + w.r):
    the count is not kept; W is a local rate of gain or loss.
S2  With the phase potential, iota q A (A self-dagger), the extra term vanishes identically: the count is kept.
S3  D(f chi) = (D f) chi + f D chi for a scalar f: an exact potential W = -(D f)/f only relabels the count.
    For W = D sigma the scale 'field' (grad w0 - d_t w ; curl w) vanishes identically.
S4  With a scale field that is not zero, two histories with the same ends give different scale factors:
    on a record the count differs by the square of the ratio and the invariant by its fourth power.
    Unit-block histories never change the invariant.  A single (null) reading stays null under any scale.
Exact arithmetic (polynomials over the Gaussian rationals; Gaussian rational blocks).  stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../ob1')
sys.path.insert(0, '../in1')
import ob1_one_block as ob
import in1_the_invariant as in1

pc, var, padd, pmul, pscale, pd, pconj = ob.pc, ob.var, ob.padd, ob.pmul, ob.pscale, ob.pd, ob.pconj
ZERO, ONE, C, IOTA = ob.ZERO, ob.ONE, ob.C, ob.IOTA


def column():
    t, x, y, z = (var(i) for i in range(4))
    psi = [padd(pmul(t, x), pscale((1, 2), z)), padd(pscale((0, 1), pmul(y, y)), pscale(F(3), t))]
    return [[psi[0], ZERO], [psi[1], ZERO]]


def reading(col):
    rho = ob.bmul(col, ob.bdag(col))
    n2, r = ob.parts(rho)
    return pscale(F(2), n2), [pscale(F(2), q) for q in r]


def potential_block(w0, w):
    return ob.badd(ob.bscale(w0, ONE), ob.vec_block(w))


def balance_control():
    t, x, y, z = (var(i) for i in range(4))
    col = column()
    n, r = reading(col)
    w0 = padd(pmul(x, y), pscale(F(1, 2), t))
    w = [pmul(t, z), padd(x, pscale(F(-2), y)), pmul(y, y)]
    W = potential_block(w0, w)
    if ob.bdag(W) != W:
        raise ValueError('scale potential must be self-dagger')
    law = ob.badd(ob.D(col), ob.bmul(W, col))
    lhs = padd(pd(n, 0), ob.div(r))
    m = ob.badd(ob.bmul(ob.bdag(col), law), ob.bmul(ob.bdag(law), col))
    loss = padd(pmul(w0, n), padd(padd(pmul(w[0], r[0]), pmul(w[1], r[1])), pmul(w[2], r[2])))
    if lhs != padd(m[0][0], pscale(F(2), loss), -1):
        raise ValueError('count balance with a scale potential failed')
    if loss == ZERO:
        raise ValueError('witness has no loss')
    # phase potential: the count is kept
    A = potential_block(padd(pmul(t, t), z), [x, pmul(y, z), pscale(F(3), t)])
    phase_law = ob.badd(ob.D(col), ob.bscale(pc(0, 2), ob.bmul(A, col)))           # q = 2
    m2 = ob.badd(ob.bmul(ob.bdag(col), phase_law), ob.bmul(ob.bdag(phase_law), col))
    if lhs != m2[0][0]:
        raise ValueError('a phase potential must keep the count')
    return dict(scale_potential='count changes at the rate -2 (w0 n + w.r)', phase_potential='count kept')


def exact_control():
    t, x, y, z = (var(i) for i in range(4))
    col = column()
    f = padd(padd(pmul(x, t), pscale(F(2), pmul(y, z))), pc(5))
    fcol = [[pmul(f, col[i][j]) for j in range(2)] for i in range(2)]
    Df = ob.D(ob.bscale(f, ONE))
    if ob.D(fcol) != ob.badd(ob.bmul(Df, col), ob.bscale(f, ob.D(col))):
        raise ValueError('product rule for a scalar factor failed')
    sigma = padd(pmul(pmul(x, x), t), pmul(y, z))
    w0, w = pd(sigma, 0), ob.grad(sigma)
    Es = [padd(pd(w0, i+1), pd(w[i], 0), -1) for i in range(3)]
    Bs = ob.curl(w)
    if any(p != ZERO for p in Es+Bs):
        raise ValueError('an exact scale potential must have no scale field')
    w_curved = [pscale(F(-1, 2), y), pscale(F(1, 2), x), ZERO]
    if ob.curl(w_curved) == [ZERO, ZERO, ZERO]:
        raise ValueError('witness scale potential should have a field')
    return dict(product_rule=True, exact_potential_field='0', witness_field='1 (constant, along the third cut)')


def history_control():
    """Link scales around a square: bottom 2, right 3, left 1, top 2 (all read forwards)."""
    rho = in1.madd(in1.scale(F(1, 3), in1.tensor((in1.c(3), in1.c(1, 2)))),
                   in1.scale(F(2, 3), in1.tensor((in1.c(1, -1), in1.c(2)))))
    n0, r0 = in1.readings(rho)
    q0 = in1.form(n0, r0)

    def scaled(factor):
        M = in1.scale(factor, in1.ONE)
        return in1.readings(in1.mm(in1.mm(M, rho), in1.dagger(M)))
    path_a, path_b = F(2)*F(3), F(1)*F(2)
    na, ra = scaled(path_a)
    nb, rb = scaled(path_b)
    ratio = path_a/path_b
    if na/nb != ratio**2 or in1.form(na, ra)/in1.form(nb, rb) != ratio**4:
        raise ValueError('history dependence of count and invariant failed')
    # unit-block histories with the same ends... or any: the invariant never changes
    g1 = [[in1.c(2), in1.Z], [in1.Z, in1.c(F(1, 2))]]
    g2 = [[in1.c(F(5, 4)), in1.c(F(3, 4))], [in1.c(F(3, 4)), in1.c(F(5, 4))]]
    for seq in ((g1, g2), (g2, g1)):
        cur = rho
        for g in seq:
            cur = in1.mm(in1.mm(g, cur), in1.dagger(g))
        n, r = in1.readings(cur)
        if in1.form(n, r) != q0:
            raise ValueError('a unit-block history changed the invariant')
    # a single reading stays null under any scale
    single = in1.tensor((in1.c(2, 1), in1.c(1, -3)))
    M = in1.scale(F(7, 3), in1.ONE)
    n, r = in1.readings(in1.mm(in1.mm(M, single), in1.dagger(M)))
    if in1.form(n, r) != 0:
        raise ValueError('a single reading must stay null')
    return dict(two_histories=[str(path_a), str(path_b)], count_ratio=str(ratio**2), invariant_ratio=str(ratio**4),
                unit_block_histories='invariant unchanged', single_reading='null under any scale')


def run():
    return dict(balance=balance_control(), exact=exact_control(), histories=history_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('SC1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
