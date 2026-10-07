"""QC3: one arrow, two units - the even part is heat, the odd part is phase.

For a symmetric record (turns +v and -v of one native arrow, GE2-T2/T3):
  E = (U_v + U_-v)/2   record mean (heat side),   O = (U_v - U_-v)/2   (phase side),
  E^2 - O^2 = I,   record memory M = I - E^2 = -O^2.
Finite-turn laws for circular contents q and for su(2) contents j = 1/2, 1.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
import json


# ---------------------------------------------------------- circular contents
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


def circular_control(t, contents):
    """Half-turn with tan(theta/4) = t: (c2, s2) = cos, sin of theta/2, rational."""
    c2, s2 = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    c, s = c2*c2-s2*s2, 2*c2*s2                            # the turn theta
    rows = []
    for q in contents:
        re, im = F(1), F(0)                                # (c + s R)^q by multiplication
        for _ in range(q):
            re, im = re*c-im*s, re*s+im*c
        E, O = re, im                                      # even part, odd part (times R)
        if E != cheb_T(q, c):
            raise ValueError('even part is not the q-fold cosine')
        if E*E+O*O != 1:
            raise ValueError('E^2 - O^2 = I failed (O carries R, R^2 = -1)')
        memory = 1-E*E
        if memory != O*O:
            raise ValueError('record memory is not minus the squared odd part')
        character = cheb_U(q-1, c2)                        # sin(q theta/2)/sin(theta/2)
        if 1-E != (1-c)*character**2:
            raise ValueError('heat defect is not unit defect times character squared')
        if O != s*cheb_U(q-1, c):
            raise ValueError('odd part is not unit odd part times character')
        if not 0 <= q*q-character**2:
            raise ValueError('finite-turn rate exceeds the square of the content')
        rows.append(dict(content=q, even=str(E), odd=str(O), memory=str(memory),
                         defect_ratio=str(character**2), limit=q*q))
    return dict(tan_quarter_turn=str(t), rows=rows)


def circular_limit(q):
    """(1 - cos q theta)/(1 - cos theta) increases to q^2 along tan(theta/4) = 1/n."""
    values = []
    for n in (2, 4, 8, 16, 64, 256):
        t = F(1, n)
        c2 = (1-t*t)/(1+t*t)
        values.append(cheb_U(q-1, c2)**2)
    if any(b <= a for a, b in zip(values, values[1:])) or values[-1] >= q*q:
        raise ValueError('rates do not increase to the square of the content')
    return dict(content=q, last_gap=str(q*q-values[-1]))


# -------------------------------------------------------------- su(2) contents
def qmul(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


def turn(axis, r):
    """GE2 (4): q(v) = (1 - r^2/16 - v/2)/(1 + r^2/16), v = r e_axis."""
    den = 1+r*r/16
    out = [(1-r*r/16)/den, F(0), F(0), F(0)]
    out[axis+1] = -(r/2)/den
    return tuple(out)


def left_matrix(q):
    cols = []
    for k in range(4):
        e = tuple(F(int(i == k)) for i in range(4))
        cols.append(qmul(q, e))
    return [[cols[j][i] for j in range(4)] for i in range(4)]


def mat_mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def mat_add(a, b, s=1):
    return [[x+s*y for x, y in zip(r1, r2)] for r1, r2 in zip(a, b)]


def mat_scale(c, a):
    return [[c*x for x in row] for row in a]


def ident(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


PAIRS = list(combinations_with_replacement(range(4), 2))


def quadratic_lift(L):
    """Action f(x) -> f(Lx) on quadratic monomials x_i x_j (10 x 10)."""
    out = [[F(0)]*10 for _ in range(10)]
    for col, (i, j) in enumerate(PAIRS):
        # x_i x_j evaluated at Lx
        poly = {}
        for a in range(4):
            for b in range(4):
                key = (min(a, b), max(a, b))
                poly[key] = poly.get(key, F(0))+L[i][a]*L[j][b]
        for row, key in enumerate(PAIRS):
            out[row][col] = poly.get(key, F(0))
    return out


def su2_control(r):
    den = 1+r*r/16
    c = (1-r*r/16)/den                                     # cos(phi/2)
    h = r*r/2                                              # GE2 (10): half the trace of Q
    # content 1/2: degree-one coefficients, single axis and isotropic record
    E_axis = mat_scale(F(1, 2), mat_add(left_matrix(turn(0, r)), left_matrix(turn(0, -r))))
    O_axis = mat_scale(F(1, 2), mat_add(left_matrix(turn(0, r)), left_matrix(turn(0, -r)), -1))
    if mat_add(mat_mul(E_axis, E_axis), mat_mul(O_axis, O_axis), -1) != ident(4):
        raise ValueError('E^2 - O^2 = I failed on content 1/2')
    memory = mat_add(ident(4), mat_mul(E_axis, E_axis), -1)
    if memory != mat_scale(F(-1), mat_mul(O_axis, O_axis)):
        raise ValueError('memory is not minus the squared odd part on content 1/2')
    S_half = [[F(0)]*4 for _ in range(4)]
    S_one = [[F(0)]*10 for _ in range(10)]
    for axis in range(3):
        for sign in (1, -1):
            L = left_matrix(turn(axis, sign*r))
            S_half = mat_add(S_half, mat_scale(F(1, 6), L))
            S_one = mat_add(S_one, mat_scale(F(1, 6), quadratic_lift(L)))
    if S_half != mat_scale(c, ident(4)):
        raise ValueError('record mean on content 1/2 is not the normalized character')
    lam = (4*c*c-1)/3                                      # chi_1/3
    if mat_mul(mat_add(S_one, ident(10), -1), mat_add(S_one, mat_scale(lam, ident(10)), -1)) != [[F(0)]*10]*10:
        raise ValueError('record mean on degree two is not {1, chi_1/3}')
    if sum(S_one[i][i] for i in range(10)) != 1+9*lam:
        raise ValueError('content 1 does not have multiplicity nine')
    rate_half, rate_one = (1-c)/h, (1-lam)/h
    if rate_half != F(1, 4)/den or rate_one != F(2, 3)/den**2:
        raise ValueError('finite-turn Casimir law failed')
    if rate_one/rate_half != F(8, 3)/den:
        raise ValueError('Casimir ratio failed')
    return dict(turn_size=str(r), rate_content_half=str(rate_half), rate_content_one=str(rate_one),
                limit_half='1/4 = (3/4)/3', limit_one='2/3 = (2)/3', casimir_half='3/4', casimir_one='2',
                casimir_gap='5/4', ratio=str(rate_one/rate_half))


def run():
    return dict(circular=[circular_control(t, range(1, 7)) for t in (F(1, 8), F(1, 3), F(2, 5))],
                circular_limits=[circular_limit(q) for q in (2, 3, 5)],
                su2=[su2_control(r) for r in (F(1), F(1, 4), F(1, 32))])


if __name__ == '__main__':
    res = run()
    json.dump(res, open('QC3_RESULT.json', 'w'), indent=1)
    for row in res['circular'][0]['rows'][:4]:
        print(row)
    print(res['circular_limits'])
    for row in res['su2']:
        print(row)
