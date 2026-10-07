"""IN1: the invariant - observed and lost together keep the uncut, in every frame.

Cut-complex two-component carrier; cuts C1 = K, C2 = RK, C3 = iota R (FR1-T3).
A reading tensor is rho = (n + r1 C1 + r2 C2 + r3 C3)/2 with real n, r_i (self-dagger).
T1  det rho = (n^2 - r1^2 - r2^2 - r3^2)/4.
T2  single reading rho = psi psi^dagger:  det = 0,  n^2 = r1^2 + r2^2 + r3^2.
    Lost to cut i  =  observed by the other two:  1 - (r_i/n)^2 = sum_{j != i} (r_j/n)^2.
T3  two cuts are complete for a real psi (r3 = 0) and incomplete for a cut-complex psi; three always are.
T4  every frame change rho -> g rho g^dagger with det g = 1 keeps n^2 - r.r.
T5  a record p rho_1 + (1-p) rho_2 of two single readings has
        n^2 - r.r = 2 p (1-p) (n1 n2 - r1.r2) >= 0,
    zero iff p in {0, 1} or the two readings are parallel:  observed^2 + unrecoverable = uncut^2.
Exact arithmetic over the Gaussian rationals, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json

Z, O, IOTA = (F(0), F(0)), (F(1), F(0)), (F(0), F(1))


def c(re, im=0):
    return (F(re), F(im))


def cadd(a, b):
    return (a[0]+b[0], a[1]+b[1])


def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def conj(a):
    return (a[0], -a[1])


def mm(a, b):
    out = [[Z, Z], [Z, Z]]
    for i in range(2):
        for j in range(2):
            s = Z
            for k in range(2):
                s = cadd(s, cmul(a[i][k], b[k][j]))
            out[i][j] = s
    return out


def dagger(a):
    return [[conj(a[j][i]) for j in range(2)] for i in range(2)]


def scale(x, a):
    return [[cmul(c(x), v) for v in row] for row in a]


def madd(a, b):
    return [[cadd(a[i][j], b[i][j]) for j in range(2)] for i in range(2)]


def det(a):
    d = cadd(cmul(a[0][0], a[1][1]), cmul(c(-1), cmul(a[0][1], a[1][0])))
    return d


def trace(a):
    return cadd(a[0][0], a[1][1])


ONE = [[O, Z], [Z, O]]
C1 = [[O, Z], [Z, c(-1)]]                      # K
C2 = [[Z, O], [O, Z]]                          # RK
C3 = [[Z, c(0, -1)], [c(0, 1), Z]]             # iota R
CUTS = [C1, C2, C3]


def tensor(psi):
    return [[cmul(psi[i], conj(psi[j])) for j in range(2)] for i in range(2)]


def readings(rho):
    n = trace(rho)
    r = [trace(mm(rho, C)) for C in CUTS]
    if n[1] != 0 or any(x[1] != 0 for x in r):
        raise ValueError('readings of a self-dagger tensor must be real')
    rebuilt = scale(F(1, 2), madd(scale(n[0], ONE), madd(scale(r[0][0], C1), madd(scale(r[1][0], C2), scale(r[2][0], C3)))))
    if rebuilt != rho:
        raise ValueError('tensor is not (n + r.C)/2')
    return n[0], [x[0] for x in r]


def form(n, r):
    return n*n-sum(x*x for x in r)


STATES = [(c(3), c(1, 0)), (c(1), c(0, 1)), (c(2, 1), c(1, -3)), (c(1, 2), c(0)), (c(1, 1), c(1, -1)),
          (c(F(1, 2)), c(F(2, 3), F(1, 5)))]


def single_control():
    rows = []
    for psi in STATES:
        rho = tensor(psi)
        n, r = readings(rho)
        d = det(rho)
        if d != Z or 4*d[0] != form(n, r) or form(n, r) != 0:
            raise ValueError('single reading is not null')
        for i in range(3):
            lost = 1-(r[i]/n)**2
            if lost != sum((r[j]/n)**2 for j in range(3) if j != i):
                raise ValueError('lost to one cut is not what the other two observe')
        real = all(x[1] == 0 for x in psi)
        two = r[0]**2+r[1]**2 == n*n
        if real and (r[2] != 0 or not two):
            raise ValueError('a real state must be complete with two cuts')
        rows.append(dict(n=str(n), readings=[str(x) for x in r], real=real, two_cuts_complete=two))
    if all(r['two_cuts_complete'] for r in rows):
        raise ValueError('some cut-complex state must need the third cut')
    return rows


FRAMES = [[[c(2), Z], [Z, c(F(1, 2))]],                                     # boost along C1
          [[c(F(3, 5)), c(F(-4, 5))], [c(F(4, 5)), c(F(3, 5))]],            # turn
          [[O, IOTA], [Z, O]],                                              # shear with iota
          [[c(F(3, 5), F(4, 5)), Z], [Z, c(F(3, 5), F(-4, 5))]],            # phase turn about C1
          [[c(F(5, 4)), c(F(3, 4))], [c(F(3, 4)), c(F(5, 4))]]]             # boost along C2


def frame_control():
    rows = []
    tensors = [tensor(p) for p in STATES[:3]]
    tensors.append(madd(scale(F(1, 3), tensor(STATES[0])), scale(F(2, 3), tensor(STATES[2]))))
    for rho in tensors:
        n0, r0 = readings(rho)
        q0 = form(n0, r0)
        cur = rho
        changed = False
        for g in FRAMES:
            if det(g) != O:
                raise ValueError('frame change must have determinant one')
            cur = mm(mm(g, cur), dagger(g))
            n, r = readings(cur)
            if form(n, r) != q0:
                raise ValueError('n^2 - r.r is not frame independent')
            changed = changed or (n, r) != (n0, r0)
        if not changed:
            raise ValueError('readings should change with the frame')
        rows.append(dict(invariant=str(q0), first=[str(n0)]+[str(x) for x in r0],
                         last=[str(n)]+[str(x) for x in r]))
    return rows


def record_control():
    rows = []
    pairs = [(STATES[0], STATES[1]), (STATES[2], STATES[4]), (STATES[3], STATES[5]), (STATES[0], (c(6), c(2)))]
    for a, b in pairs:
        ra, rb = tensor(a), tensor(b)
        na, va = readings(ra)
        nb, vb = readings(rb)
        gap = na*nb-sum(x*y for x, y in zip(va, vb))
        if gap < 0:
            raise ValueError('two single readings must have n1 n2 >= r1.r2')
        for p in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
            rho = madd(scale(p, ra), scale(1-p, rb))
            n, r = readings(rho)
            q = form(n, r)
            if q != 2*p*(1-p)*gap or q != 4*det(rho)[0]:
                raise ValueError('record invariant is not 2 p (1-p) (n1 n2 - r1.r2)')
            if q < 0:
                raise ValueError('unrecoverable part must be nonnegative')
            observed = sum(x*x for x in r)
            if observed+q != n*n:
                raise ValueError('observed + unrecoverable != uncut')
        rows.append(dict(gap=str(gap), invariant_at_half=str(gap/2), parallel=(gap == 0)))
    if not rows[-1]['parallel'] or rows[0]['parallel']:
        raise ValueError('parallel readings must lose nothing; non-parallel ones must')
    return rows


def run():
    return dict(single=single_control(), frames=frame_control(), records=record_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('IN1_RESULT.json', 'w'), indent=1)
    for r in res['single']:
        print(r)
    for r in res['frames']:
        print(r)
    print(res['records'])
