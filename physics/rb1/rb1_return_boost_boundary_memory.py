"""RB1: the recognised return is a boost; the horizon end keeps boundary memory, the flat end does not.

Stdlib only, exact rationals.  Sources read (unchanged):
  Recognition-Kernel-Framework research/recognition_return  R1 (RR1, RA1-RA3, RC1) and r2 (RD1, RD2, RI1, RI2)
  Publications papers/native-critical-response               CR-1, CR-2
  this folder's line                                           GB1 (gravity element), LT1 (tower generation), MO1
Native algebra: R^2 = -1, K^2 = 1, KR = -RK, L = KR.  Elements are (1, R, K, L) coefficient 4-tuples.
"""
from fractions import Fraction as F
import json

ONE = (F(1), F(0), F(0), F(0))


def mul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    # R^2=-1, K^2=1, L^2=1, KR=L, RK=-L, RL=K, LR=-K, KL=R, LK=-R
    return (a0*b0 - a1*b1 + a2*b2 + a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 + a2*b0 + a1*b3 - a3*b1,
            a0*b3 + a3*b0 + a2*b1 - a1*b2)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scal(s, a):
    return tuple(s*x for x in a)


def el(u, v):
    """u + v L"""
    return (F(u), F(0), F(0), F(v))


R = (F(0), F(1), F(0), F(0))
K = (F(0), F(0), F(1), F(0))
L = mul(K, R)


def inv_line(a):
    """inverse of u + vL when u^2 != v^2"""
    u, v = a[0], a[3]
    if a[1] or a[2] or u*u == v*v:
        raise ValueError('not invertible on the response line')
    d = u*u - v*v
    return el(u/d, -v/d)


def append_cell(F_tail, p, q, b0=F(1), b1=F(1)):
    """R2 RD1: eliminate an invertible tail; edges B_h = b_h R, C_h = c_h K with b0c0 = p, b1c1 = q."""
    B0, C0 = scal(b0, R), scal(p/b0, K)
    B1, C1 = scal(b1, R), scal(q/b1, K)
    D = add(add(ONE, scal(-1, mul(mul(C1, F_tail), B1))), scal(-1, mul(B0, C0)))
    return add(ONE, mul(mul(C0, inv_line(D)), B0))


def chain_return(cells, tail=F(0)):
    """boundary return 1 + x_0 L of a paired chain a_0..a_{n-1} closed by 1 + tail L (R2 RD2)."""
    Fc = el(1, tail)
    for a in reversed(cells):
        Fc = append_cell(Fc, a, a)
        if Fc[0] != 1 or Fc[1] or Fc[2]:
            raise ValueError('left the unit response line')
    return Fc[3]


def cells_of_profile(xs):
    """R2 SY2 inverse design: a_j = x_j / (1 - x_j x_{j+1})."""
    return [xs[j] / (1 - xs[j]*xs[j+1]) for j in range(len(xs) - 1)]


def tower(x):
    """LT1 generation (squaring of the unit element): x -> 2x/(1+x^2)."""
    return 2*x / (1 + x*x)


def interval(cells):
    """R2 RI1: all responses for tails in [0, infinity]; returns (end0, end_inf, width, C, D)."""
    A, B, C, D = F(1), F(0), F(0), F(1)
    for a in cells:
        r = 1/a
        A, B, C, D = B, A + B*r, D, C + D*r
    return B/D, A/C, abs(A/C - B/D), C, D


def run():
    res = {}
    # T1  uniform chain (R1): return = boost with q = x/(1-x^2)
    x = F(2, 3)
    q = x / (1 - x*x)
    Fx = el(1, x)
    C, B = scal(q, K), R
    if add(ONE, mul(mul(mul(C, Fx), B), Fx)) != Fx:
        raise ValueError('first-return law fails')
    x2 = tower(x)
    if 4*q*q*(1 - x2*x2) != x2*x2:      # (2q)^2 = x2^2/(1-x2^2): 2q = proper speed of the next tower level
        raise ValueError('coupling law fails')
    y = F(1, 5)
    prod = mul(el(1, x), el(1, y))
    if prod != scal(1 + x*y, el(1, (x + y)/(1 + x*y))):
        raise ValueError('composition law fails')
    res['uniform'] = {'x': str(x), 'q': str(q), 'unit_norm': str(1 - x*x), 'next_level': str(x2)}

    # T2  gravity element as an R2 chain along the tower, closed by its own next value
    xs = [F(1, 3)]
    for _ in range(7):
        xs.append(tower(xs[-1]))
    cells = cells_of_profile(xs)
    if chain_return(cells, tail=xs[-1]) != xs[0]:
        raise ValueError('profile is not its own chain return')
    recip = [1/a for a in cells]
    for j, r in enumerate(recip):
        if r != (1 - xs[j]**2) / (xs[j]*(1 + xs[j]**2)):
            raise ValueError('reciprocal cell law fails')
    for j in range(len(xs) - 1):
        if 1 - xs[j+1]**2 != ((1 - xs[j]**2)/(1 + xs[j]**2))**2:
            raise ValueError('clock factor does not square')

    # T3  inward: reciprocal sum bounded, interval stays open (retained memory)
    widths = [interval(cells[:n])[2] for n in range(1, len(cells) + 1)]
    if not all(widths[i+1] <= widths[i] for i in range(len(widths) - 1)):
        raise ValueError('intervals not nested')
    e0, einf, w, _, _ = interval(cells)
    S = sum(recip)
    if not (e0 <= xs[0] <= einf or einf <= xs[0] <= e0):
        raise ValueError('profile outside its interval')
    res['inward'] = {'x0': str(xs[0]), 'levels': len(cells), 'reciprocal_sum': float(S),
                     'ends': [float(e0), float(einf)], 'width': float(w),
                     'width_last_two': [float(widths[-2]), float(widths[-1])]}
    if S >= 4 or w < F(1, 10) or widths[-2] - widths[-1] > F(1, 10**15):
        raise ValueError('inward memory not retained')
    # closures: wall (tail 0), the profile's own (tail -> 1), fully open (tail -> infinity)
    res['inward']['closures'] = {'tail_0': float(chain_return(cells, F(0))),
                                 'tail_1': float(chain_return(cells, F(1))),
                                 'tail_own': float(chain_return(cells, xs[-1]))}

    # T4  outward: x_j = 1/(j+3) (radius (j+3)^2 r_s); reciprocal sum diverges, interval closes
    out = [F(1, j + 3) for j in range(402)]
    oc = cells_of_profile(out)
    w50, w100, w400 = (interval(oc[:n])[2] for n in (50, 100, 400))
    if not (w400 < w100 < w50 and w400 < F(1, 10**6)):
        raise ValueError('outward interval does not close')
    if not all(1/a >= 2 for a in oc):
        raise ValueError('outward reciprocal cells not bounded below')
    res['outward'] = {'width_50': float(w50), 'width_100': float(w100), 'width_400': float(w400)}

    # T5  the inverse at closure (CR-1, CR-2): pole only in P-
    Pp, Pm = el(F(1, 2), F(1, 2)), el(F(1, 2), F(-1, 2))
    if mul(Pp, Pm) != (0, 0, 0, 0) or mul(Pp, Pp) != Pp or mul(R, Pp) != mul(Pm, R):
        raise ValueError('channel algebra fails')
    for xx in (F(9, 10), F(99, 100), F(999, 1000)):
        G = inv_line(el(1, xx))
        if G != add(scal(1/(1 + xx), Pp), scal(1/(1 - xx), Pm)):
            raise ValueError('channel split fails')
        if mul(mul(Pp, G), Pp) != scal(1/(1 + xx), Pp):
            raise ValueError('P+ reading is not regular')
    res['closure'] = {'F_at_closure': '2 P+', 'P+_reading_of_inverse_at_closure': '1/2',
                      'P-_coefficient': '1/(1-x)'}
    return res


if __name__ == '__main__':
    out = run()
    with open('RB1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
