"""QC4: return is memory.

T1  curvature of a mean of flat readings = minus the wedge-variance of the readings
    (PH3's share law is the two-reading case; four compass readings checked).
T2  GE2's record memory M = I - S^dagger S is the variance of the records.
T3  for a loop arrow read without its orientation, the Wilson density is the
    record defect 1 - E and the squared field strength is the record memory -O^2,
    with M = (1 - E)(1 + E).
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../ph3')
sys.path.insert(0, '../qc3')
import ph3_two_readings_transport as ph3
import qc3_even_odd_arrow as qc3

mm, add, scal, inv = ph3.mm, ph3.add, ph3.scal, ph3.inv
ZERO = [[F(0), F(0)], [F(0), F(0)]]


def comm(a, b):
    return add(mm(a, b), mm(b, a), -1)


# ------------------------------------------------------------ T1: readings
def readings(x, y):
    """Four flat readings of a probe: frames P with frozen coordinates P w.

    U: (ds, dv)   G: (dT, -dP) = H w   F: (dT, dv)   Hc: (ds, -dP)
    Each returns (A_x, A_y, d_x A_y - d_y A_x) from P and its derivatives.
    """
    H, Hx, Hy, Hxy, Hyx = ph3.field(x, y)
    one, zero = F(1), F(0)

    def frame(mask):
        def pick(M, const):
            return [[M[i][j] if mask[i] else const[i][j] for j in range(2)] for i in range(2)]
        ident, nil = [[one, zero], [zero, one]], ZERO
        return (pick(H, ident), pick(Hx, nil), pick(Hy, nil), pick(Hxy, nil), pick(Hyx, nil))

    out = {}
    for name, mask in (('U', (0, 0)), ('G', (1, 1)), ('F', (1, 0)), ('Hc', (0, 1))):
        P, Px, Py, Pxy, Pyx = frame(mask)
        Pi = inv(P)
        Ax, Ay = mm(Pi, Px), mm(Pi, Py)
        dxAy = add(mm(Pi, Pyx), mm(mm(Pi, Px), mm(Pi, Py)), -1)
        dyAx = add(mm(Pi, Pxy), mm(mm(Pi, Py), mm(Pi, Px)), -1)
        curl = add(dxAy, dyAx, -1)
        if add(curl, comm(Ax, Ay)) != ZERO:
            raise ValueError('reading %s is not flat' % name)
        out[name] = (Ax, Ay, curl)
    return out, H


def mixture_curvature(read, weights):
    names = list(weights)
    Ax = Ay = curl = ZERO
    for n in names:
        Ax = add(Ax, scal(weights[n], read[n][0]))
        Ay = add(Ay, scal(weights[n], read[n][1]))
        curl = add(curl, scal(weights[n], read[n][2]))
    direct = add(curl, comm(Ax, Ay))
    variance = ZERO
    for i, s in enumerate(names):
        for t in names[i+1:]:
            dx = add(read[s][0], read[t][0], -1)
            dy = add(read[s][1], read[t][1], -1)
            variance = add(variance, scal(weights[s]*weights[t], comm(dx, dy)))
    if direct != scal(F(-1), variance):
        raise ValueError('curvature of the mean is not minus the variance of the readings')
    return direct


def rotation_part(H, curv):
    """sc(R H F): proportional to the native marker f (common factor 1/sqrt(det H))."""
    R = [[F(0), F(-1)], [F(1), F(0)]]
    m = mm(mm(R, H), curv)
    return (m[0][0]+m[1][1])/2


def readings_control():
    rows = []
    for x, y in ((F(1, 2), F(1, 3)), (F(-2, 5), F(3, 4)), (F(1), F(-1, 7))):
        read, H = readings(x, y)
        native = mixture_curvature(read, {'U': F(1, 2), 'G': F(1, 2)})
        base = rotation_part(H, native)
        if base == 0:
            raise ValueError('witness has no return')
        row = dict(point=[str(x), str(y)], shares={})
        for label, w in (('U:G = 3:1', {'U': F(3, 4), 'G': F(1, 4)}),
                         ('F:Hc = 1:1', {'F': F(1, 2), 'Hc': F(1, 2)}),
                         ('four corners equal', {'U': F(1, 4), 'G': F(1, 4), 'F': F(1, 4), 'Hc': F(1, 4)}),
                         ('U:F:G = 1:1:2', {'U': F(1, 4), 'F': F(1, 4), 'G': F(1, 2)})):
            row['shares'][label] = str(rotation_part(H, mixture_curvature(read, w))/base)
        if row['shares']['U:G = 3:1'] != '3/4':
            raise ValueError('two-reading case is not the share law 4 lambda (1 - lambda)')
        if row['shares']['F:Hc = 1:1'] != '0' or row['shares']['four corners equal'] != '1/2':
            raise ValueError('compass law failed: mixed diagonal 0, four corners 1/2')
        if rotation_part(H, mixture_curvature(read, {'F': F(1, 3), 'Hc': F(2, 3)})) != 0:
            raise ValueError('mixed diagonal has a return at unequal share')
        if mixture_curvature(read, {'F': F(1, 2), 'Hc': F(1, 2)}) == ZERO:
            raise ValueError('mixed diagonal should still carry non-rotational curvature')
        edges = sum(rotation_part(H, mixture_curvature(read, {s: F(1, 2), t: F(1, 2)}))
                    for s, t in (('U', 'F'), ('U', 'Hc'), ('G', 'F'), ('G', 'Hc')))
        if edges != base:
            raise ValueError('the four edges do not sum to the diagonal')
        for pure in ('U', 'G', 'F', 'Hc'):
            if mixture_curvature(read, {pure: F(1)}) != ZERO:
                raise ValueError('a pure reading has a return')
        rows.append(row)
    return rows


# ------------------------------------------------- T2: memory is a variance
def memory_control():
    rows = []
    for weights, turns in (((F(1, 2), F(1, 2)), (F(1, 3), F(-1, 3))),
                           ((F(1, 5), F(3, 10), F(1, 2)), (F(1, 8), F(2, 5), F(-1, 3)))):
        recs = []
        for t in turns:
            c2, s2 = (1-t*t)/(1+t*t), 2*t/(1+t*t)
            recs.append((c2*c2-s2*s2, 2*c2*s2))
        S = (sum(p*c for p, (c, s) in zip(weights, recs)), sum(p*s for p, (c, s) in zip(weights, recs)))
        memory = 1-(S[0]**2+S[1]**2)
        variance = sum(p*((c-S[0])**2+(s-S[1])**2) for p, (c, s) in zip(weights, recs))
        if memory != variance or memory <= 0:
            raise ValueError('record memory is not the variance of the records')
        rows.append(dict(records=len(recs), memory=str(memory)))
    return rows


# ----------------------------------------------- T3: the loop arrow, SU(2)
def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])


def loop_control(r):
    links = [qc3.turn(0, r), qc3.turn(1, r/2), qc3.turn(2, -r), qc3.turn(0, r/3)]
    U = qc3.qmul(qc3.qmul(links[0], links[1]), qc3.qmul(qconj(links[2]), qconj(links[3])))
    back = qconj(U)
    E = tuple((a+b)/2 for a, b in zip(U, back))
    O = tuple((a-b)/2 for a, b in zip(U, back))
    if E[1:] != (0, 0, 0) or O[0] != 0:
        raise ValueError('even part is not the scalar, odd part not the vector')
    W = E[0]
    field2 = sum(v*v for v in O[1:])                       # -O^2
    if 1-W*W != field2:
        raise ValueError('record memory is not the squared odd part')
    if field2 != (1-W)*(1+W):
        raise ValueError('memory is not defect times (1 + E)')
    g = qc3.turn(1, F(3, 7))
    moved = qc3.qmul(qc3.qmul(g, U), qconj(g))
    if moved[0] != W or sum(v*v for v in moved[1:]) != field2:
        raise ValueError('defect or memory is not gauge invariant')
    reverse = qconj(U)
    if reverse[0] != W or sum(v*v for v in reverse[1:]) != field2:
        raise ValueError('defect or memory sees the orientation')
    if U[1:] == reverse[1:]:
        raise ValueError('witness loop has no orientation content')
    return dict(turn_size=str(r), wilson_defect=str(1-W), memory=str(field2),
                memory_over_twice_defect=str(field2/(2*(1-W))))


# ------------------------------------- T3 again where E is not a scalar: SO(3)
def cayley(a, b, c):
    A = [[F(0), -c, b], [c, F(0), -a], [-b, a, F(0)]]
    n = a*a+b*b+c*c
    A2 = qc3.mat_mul(A, A)
    I = qc3.ident(3)
    # (I + A)(I - A)^-1 = I + 2(A + A^2)/(1 + n) for a 3x3 antisymmetric A
    return qc3.mat_add(I, qc3.mat_scale(F(2)/(1+n), qc3.mat_add(A, A2)))


def matrix_loop_control():
    U = qc3.mat_mul(cayley(F(1, 2), F(1, 3), F(0)), cayley(F(0), F(1, 5), F(2, 3)))
    Ut = [list(r) for r in zip(*U)]
    if qc3.mat_mul(U, Ut) != qc3.ident(3):
        raise ValueError('loop arrow is not an isometry')
    E = qc3.mat_scale(F(1, 2), qc3.mat_add(U, Ut))
    O = qc3.mat_scale(F(1, 2), qc3.mat_add(U, Ut, -1))
    if qc3.mat_add(qc3.mat_mul(E, E), qc3.mat_mul(O, O), -1) != qc3.ident(3):
        raise ValueError('E^2 - O^2 = I failed')
    M = qc3.mat_add(qc3.ident(3), qc3.mat_mul(E, E), -1)
    if M != qc3.mat_scale(F(-1), qc3.mat_mul(O, O)):
        raise ValueError('memory is not minus the squared odd part')
    if M != qc3.mat_mul(qc3.mat_add(qc3.ident(3), E, -1), qc3.mat_add(qc3.ident(3), E)):
        raise ValueError('memory is not defect times (I + E)')
    if E[0][1] == 0 and E[0][2] == 0 and E[1][2] == 0:
        raise ValueError('witness even part is scalar')
    tr = lambda m: sum(m[i][i] for i in range(3))
    return dict(mean_defect=str(1-tr(E)/3), mean_memory=str(tr(M)/3))


def run():
    return dict(readings=readings_control(), memory=memory_control(),
                loops=[loop_control(r) for r in (F(1), F(1, 4), F(1, 32))],
                matrix_loop=matrix_loop_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('QC4_RESULT.json', 'w'), indent=1)
    for row in res['readings']:
        print(row)
    print(res['memory'])
    for row in res['loops']:
        print({k: (v if len(v) < 40 else str(float(F(v)))) for k, v in row.items()})
    print(res['matrix_loop'])
