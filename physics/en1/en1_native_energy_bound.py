"""EN1: the native mass bound on recognition energy, in all three sectors; for a unit boost it is the Doppler square.

Stdlib only, exact rationals.  Sources read (unchanged):
  RH-Framework theorems/T01_NATIVE_SUMMABILITY_AND_RECOGNITION_COMPLETION.md
      (2.1) gauge D, (2.4) T01-A, (4.3) energy, (5.1) weighted cut-square identity, T01-C
  this line: WQ1 (a mode is a pair, E = (w/2)(w^2+p^2)), EM1, GW1-V1, RB1, LT1.
A scalar is r + t u with u^2 = s, s in {-1, 0, +1} (turn, shear, boost sector).
"""
from fractions import Fraction as F
import itertools
import json
import random


def mul(z, w, s):
    return (z[0]*w[0] + s*z[1]*w[1], z[0]*w[1] + z[1]*w[0])


def D(z):
    return abs(z[0]) + abs(z[1])


def N_own(z, s):
    """z-dagger z in its own sector: r^2 - s t^2"""
    return z[0]**2 - s*z[1]**2


def N_pos(z):
    """the positive energy of the two readings"""
    return z[0]**2 + z[1]**2


def conv(a, x, s, n):
    """path product on the pair groupoid of n objects (arrows i<-j)"""
    out = [[(F(0), F(0)) for _ in range(n)] for _ in range(n)]
    for i, k, j in itertools.product(range(n), repeat=3):
        p = mul(a[i][k], x[k][j], s)
        out[i][j] = (out[i][j][0] + p[0], out[i][j][1] + p[1])
    return out


def mass(a):
    return sum(D(z) for row in a for z in row)


def energy(a):
    return sum(N_pos(z) for row in a for z in row)


def boost(b):
    """unit element (1 + b u)/N of the boost sector, as (numerator, N^2)"""
    return (F(1), F(b)), 1 - F(b)**2


def gain(b, w, s):
    """energy of (1 + b u) w divided by energy of w"""
    return N_pos(mul((F(1), F(b)), w, s)) / N_pos(w)


def run(seed=7, trials=300):
    rng = random.Random(seed)
    rnd = lambda: F(rng.randint(-9, 9), rng.randint(1, 9))
    res = {}
    worst = {}
    for s in (-1, 0, 1):
        best = F(0)
        for _ in range(trials):
            z, w = (rnd(), rnd()), (rnd(), rnd())
            zw = mul(z, w, s)
            if D(zw) > D(z)*D(w):
                raise ValueError('mass gauge not submultiplicative')
            if N_own(zw, s) != N_own(z, s)*N_own(w, s):
                raise ValueError('own radial square not multiplicative')
            if N_pos(zw) > D(z)**2 * N_pos(w):
                raise ValueError('scalar energy bound fails')
        for _ in range(trials // 6):
            n = 3
            a = [[(rnd(), rnd()) for _ in range(n)] for _ in range(n)]
            x = [[(rnd(), rnd()) for _ in range(n)] for _ in range(n)]
            ax = conv(a, x, s, n)
            if mass(ax) > mass(a)*mass(x):
                raise ValueError('T01-A fails')
            if energy(ax) > mass(a)**2 * energy(x):
                raise ValueError('T01-C fails')
            if energy(x):
                best = max(best, energy(ax) / (mass(a)**2 * energy(x)))
        worst[s] = best
    res['largest_ratio_seen'] = {str(s): float(v) for s, v in worst.items()}

    # the own radial square gives no bound outside the turn sector
    x, y = (F(1), F(1)), (F(1), F(-1))
    if not (N_own(x, 1) == 0 == N_own(y, 1) and N_own((x[0] + y[0], x[1] + y[1]), 1) == 4):
        raise ValueError('null pair example fails')

    # unit boost: bound = (1+b)/(1-b), attained exactly on the co-moving null reading
    rows = []
    b = F(1, 3)
    for _ in range(4):
        num, n2 = boost(b)
        bound = D(num)**2 / n2
        if bound != (1 + b)/(1 - b):
            raise ValueError('bound is not the Doppler square')
        hi = gain(b, (F(1), F(1)), 1) / n2
        lo = gain(b, (F(1), F(-1)), 1) / n2
        mid = gain(b, (F(1), F(0)), 1) / n2
        if not (hi == bound and lo == 1/bound and lo < mid < hi):
            raise ValueError('sharpness on null readings fails')
        for _ in range(40):
            w = (rnd(), rnd())
            if N_pos(w) and abs(w[0]) != abs(w[1]):
                g = gain(b, w, 1) / n2
                if not (1/bound < g < bound):
                    raise ValueError('non-null reading reaches the bound')
        rows.append((str(b), str(bound)))
        b2 = 2*b / (1 + b*b)
        if (1 + b2)/(1 - b2) != bound**2:
            raise ValueError('bound does not square along the tower')
        b = b2
    res['boost_bound_along_tower'] = rows

    # unit turn: gain is exactly 1 although the mass bound is larger; shear: bound never reached
    c, sn = F(3, 5), F(4, 5)
    for _ in range(40):
        w = (rnd(), rnd())
        if N_pos(mul((c, sn), w, -1)) != N_pos(w):
            raise ValueError('turn does not keep the energy')
        if N_pos(w) and gain(F(1, 2), w, 0) >= (1 + F(1, 2))**2:
            raise ValueError('shear reaches the bound')
    res['turn'] = {'gain': '1', 'mass_bound': str(D((c, sn))**2)}
    return res


if __name__ == '__main__':
    out = run()
    with open('EN1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
