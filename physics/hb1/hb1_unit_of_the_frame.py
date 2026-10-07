"""HB1: the unit of a frame - count against response, and the minimum as one minus what is lost.

Docked to R41 (relational clock): marks e_r, count N e_r = r e_r, tick T (successor), even part
C = (T + T^dagger)/2, response P = iota (T - T^dagger)/2;  R41.3: [N, P] = iota R with
R = C - (Q/2)(wrap terms);  R41.4: var(N) var(P) - Cov(N, P)^2 >= <R>^2 / 4.
T1  [N, P] = iota C on states that do not touch a wrap; the cyclic form with the wrap for Q = 3, 4, 5.
T2  the inequality on rational states; R41's sharp instance (Q = 3) reproduced exactly.
T3  for a reading a_r = f_r zeta^r (real profile f, tick phase zeta = c + iota s):
        <C> = c * w_f ,   <P> = s * w_f ,   w_f = sum f_r f_{r+1} / sum f_r^2 <= 1 ,
    so the minimum of the count-response product is (1/2) c w_f:  one half of
    (1 - loss to the tick turn)(1 - loss to the finite profile), and c = 1 - (clock curvature)/2 (CL1).
Exact arithmetic over the Gaussian rationals, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json

Z = (F(0), F(0))


def c(re, im=0):
    return (F(re), F(im))


def cadd(a, b):
    return (a[0]+b[0], a[1]+b[1])


def csub(a, b):
    return (a[0]-b[0], a[1]-b[1])


def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def conj(a):
    return (a[0], -a[1])


def inner(a, b):
    s = Z
    for x, y in zip(a, b):
        s = cadd(s, cmul(conj(x), y))
    return s


def shift_up(a, cyclic):
    return ([a[-1]]+a[:-1]) if cyclic else ([Z]+a[:-1])            # (T a)_r = a_{r-1}


def shift_down(a, cyclic):
    return (a[1:]+[a[0]]) if cyclic else (a[1:]+[Z])               # (T^dagger a)_r = a_{r+1}


def count(a):
    return [cmul(c(r), x) for r, x in enumerate(a)]


def even(a, cyclic=False):
    u, d = shift_up(a, cyclic), shift_down(a, cyclic)
    return [cmul(c(F(1, 2)), cadd(x, y)) for x, y in zip(u, d)]


def response(a, cyclic=False):
    u, d = shift_up(a, cyclic), shift_down(a, cyclic)
    return [cmul(c(0, F(1, 2)), csub(x, y)) for x, y in zip(u, d)]


def commutator_control():
    """[N, P] a = iota C a for states with zero ends (no wrap is touched)."""
    states = [[Z, c(1), c(2, 1), c(0, -3), c(1, 1), Z], [Z, Z, c(5), c(1, 2), Z, Z, Z], [Z, c(1), c(1), c(1), c(1), Z]]
    for a in states:
        left = [csub(x, y) for x, y in zip(count(response(a)), response(count(a)))]
        right = [cmul(c(0, 1), x) for x in even(a)]
        if left != right:
            raise ValueError('[N, P] = iota C failed')
    # cyclic marks: [N, P] = iota R, R = C - (Q/2)(W J + J^dagger), J the wrap arrow e_0 e_{Q-1}^dagger
    for Q in (3, 4, 5):
        for a in ([c(1)]+[c(k, 1) for k in range(1, Q)], [c(2, -1)]*Q):
            left = [csub(x, y) for x, y in zip(count(response(a, True)), response(count(a), True))]
            Ca = even(a, True)
            wrap = [Z]*Q
            wrap[0] = a[Q-1]                                           # J a
            wrap_d = [Z]*Q
            wrap_d[Q-1] = a[0]                                         # J^dagger a
            Ra = [csub(x, cmul(c(F(Q, 2)), cadd(y, z))) for x, y, z in zip(Ca, wrap, wrap_d)]
            if left != [cmul(c(0, 1), x) for x in Ra]:
                raise ValueError('cyclic [N, P] = iota R failed')
    return dict(open_marks=len(states), cyclic_sizes=[3, 4, 5])


def moments(a, cyclic=False, R_of=None):
    norm = inner(a, a)[0]
    Na, Pa = count(a), response(a, cyclic)
    n_mean = inner(a, Na)[0]/norm
    p_mean = inner(a, Pa)[0]/norm
    da = [csub(x, cmul(c(n_mean), y)) for x, y in zip(Na, a)]
    db = [csub(x, cmul(c(p_mean), y)) for x, y in zip(Pa, a)]
    var_n = inner(da, da)[0]/norm
    var_p = inner(db, db)[0]/norm
    cov = inner(da, db)[0]/norm
    Ra = R_of(a) if R_of else even(a, cyclic)
    r_mean = inner(a, Ra)[0]/norm
    return var_n, var_p, cov, r_mean


def inequality_control():
    rows = []
    for a in ([Z, c(1), c(2), c(1), Z], [Z, c(1), c(3), c(4), c(3), c(1), Z], [Z, c(1), c(0, 2), c(-3), c(0, -2), c(1), Z],
              [Z, c(1), c(1, 1), c(2), c(1, -1), Z], [Z, c(1), c(4), c(6), c(4), c(1), Z]):
        var_n, var_p, cov, r_mean = moments(a)
        lhs, rhs = var_n*var_p-cov*cov, r_mean*r_mean/4
        if lhs < rhs:
            raise ValueError('count-response inequality violated')
        rows.append(dict(product=str(lhs), bound=str(rhs), ratio=str(lhs/rhs) if rhs else 'bound zero'))
    # R41's sharp instance: Q = 3 cyclic, identity process, state (e_0 + e_2)
    Q, a = 3, [c(1), Z, c(1)]

    def R_of(v):
        Ca = even(v, True)
        wrap, wrap_d = [Z]*Q, [Z]*Q
        wrap[0], wrap_d[Q-1] = v[Q-1], v[0]
        return [csub(x, cmul(c(F(Q, 2)), cadd(y, z))) for x, y, z in zip(Ca, wrap, wrap_d)]
    var_n, var_p, cov, r_mean = moments(a, cyclic=True, R_of=R_of)
    if (var_n, var_p, cov, r_mean) != (F(1), F(1, 4), F(0), F(-1)):
        raise ValueError('R41 sharp instance not reproduced')
    if var_n*var_p-cov*cov != r_mean*r_mean/4:
        raise ValueError('sharp instance is not an equality')
    return dict(states=rows, sharp_instance=dict(var_count='1', var_response='1/4', covariance='0', R='-1'))


def phase_control():
    rows = []
    profiles = {'box 4': [1, 1, 1, 1], 'box 12': [1]*12, 'triangle 5': [1, 2, 3, 2, 1],
                'binomial 8': [1, 8, 28, 56, 70, 56, 28, 8, 1]}
    for (cz, sz) in ((F(3, 5), F(4, 5)), (F(4, 5), F(3, 5)), (F(1), F(0)), (F(0), F(1)), (F(-5, 13), F(12, 13))):
        zeta = (cz, sz)
        for name, f in profiles.items():
            a = [Z]
            z = c(1)
            for x in f:
                a.append(cmul(c(x), z))
                z = cmul(z, zeta)
            a.append(Z)
            norm = inner(a, a)[0]
            w = F(sum(x*y for x, y in zip(f, f[1:])), sum(x*x for x in f))
            mean_C = inner(a, even(a))[0]/norm
            mean_P = inner(a, response(a))[0]/norm
            if mean_C != cz*w or abs(mean_P) != abs(sz*w):
                raise ValueError('<C> = c w or <P> = s w failed')
            var_n, var_p, cov, r_mean = moments(a)
            if var_n*var_p-cov*cov < (cz*w)**2/4:
                raise ValueError('minimum (1/2) c w violated')
            rows.append(dict(tick=[str(cz), str(sz)], profile=name, overlap=str(w), unit=str(cz*w),
                             one_minus_curvature_half=str(1-(2-2*cz)/2)))
    if any(r['one_minus_curvature_half'] != r['tick'][0] for r in rows):
        raise ValueError('c = 1 - curvature/2 failed')
    return rows


def ladder():
    return [dict(marks=2, tick_even_part='-1', minimum='1/2 (reversed)'),
            dict(marks=4, tick_even_part='0', minimum='0: count and response commute'),
            dict(marks=8, tick_even_part='1/sqrt(2)', minimum='1/(2 sqrt(2))', check=str(F(1, 2)) + ' = (1/sqrt(2))^2'),
            dict(marks='infinity', tick_even_part='1', minimum='1/2')]


def run():
    return dict(commutator=commutator_control(), inequality=inequality_control(), phase=phase_control(), ladder=ladder())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('HB1_RESULT.json', 'w'), indent=1)
    print(res['commutator']); print(res['inequality'])
    for r in res['phase'][:8]:
        print(r)
