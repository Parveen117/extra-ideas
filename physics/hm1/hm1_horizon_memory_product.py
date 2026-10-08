"""HM1: how much an outer reading remembers of a far closure = the product of the returns along the way.

Stdlib only, exact rationals.  Sources read (unchanged):
  Recognition-Kernel-Framework research/recognition_return r2: RD2 (backward recursion), RI1 (R2.6), RI2
  this line: RB1, TS1, CO1 (expanding space as a fall frame, x = hX), LT1.
A profile is a sequence of returns x_0, x_1, ... along the cells; x = fall speed (RB1's dictionary).
"""
from fractions import Fraction as F
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / 'rb1'))
import rb1_return_boost_boundary_memory as rb  # noqa: E402


def returns_along(recips, tail):
    """x_j = 1/(r_j + x_{j+1}) for all j, closed by the given tail; returns [x_0 .. x_n]"""
    xs = [F(tail)]
    for r in reversed(recips):
        d = r + xs[0]
        if d == 0:
            raise ZeroDivisionError('pivot vanishes')
        xs.insert(0, 1/d)
    return xs


def recips_of_profile(xs):
    return [1/xs[j] - xs[j+1] for j in range(len(xs) - 1)]


def product(vals):
    p = F(1)
    for v in vals:
        p *= v
    return p


def tower_profile(x0, levels):
    xs = [F(x0)]
    for _ in range(levels):
        xs.append(rb.tower(xs[-1]))
    return xs


def sensitivity(xs):
    """|d x_0 / d tail| on the profile itself: product of x_j^2 over the cells"""
    return product(x*x for x in xs[:-1])


def run():
    res = {}
    tower = tower_profile(F(1, 3), 6)
    rec = recips_of_profile(tower)
    n = len(rec)

    # T1 exact two-closure identity (R2.6 written with the returns themselves)
    for u, v in ((F(0), F(1)), (F(1, 2), F(7)), (tower[-1], F(3))):
        xu, xv = returns_along(rec, u), returns_along(rec, v)
        if xu[0] - xv[0] != (-1)**n * (u - v) * product(xu[:-1]) * product(xv[:-1]):
            raise ValueError('two-closure identity fails')

    # T2 a cell is positive exactly while x_j x_{j+1} < 1: R2's class ends at the horizon
    for xa, xb in ((F(1, 2), F(3, 2)), (F(9, 10), F(11, 10)), (F(3, 2), F(2))):
        if (rb.cells_of_profile([xa, xb])[0] > 0) != (xa*xb < 1):
            raise ValueError('cell sign law')

    # T3 toward a horizon along the tower: memory stays; bounded by x_0^2; the centre forgets
    rows = {}
    for x0, lv in ((F(1, 3), 6), (F(1, 10), 8), (F(1, 100), 11)):
        xs = tower_profile(x0, lv)
        sv = sensitivity(xs)
        if not (0 < sv < x0*x0) or sv - sensitivity(xs + [rb.tower(xs[-1])]) > sv/10**9:
            raise ValueError('tower sensitivity')
        rows[str(x0)] = float(sv)
    res['tower_sensitivity_by_starting_x'] = rows

    # T4 the verdict depends on the spacing: x_j = 1 - 1/(j+2) reaches the horizon too, and forgets
    for m in (5, 50, 500):
        slow = [1 - F(1, j + 2) for j in range(m + 1)]
        if min(rb.cells_of_profile(slow)) <= 0 or sensitivity(slow) != F(1, (m + 1)**2):
            raise ValueError('harmonic spacing')
    res['harmonic_spacing_sensitivity'] = '1/(n+1)^2 after n cells'
    flat = [F(1, j + 3) for j in range(40)]
    if sensitivity(flat) > F(1, 10**90):
        raise ValueError('flat end sensitivity')

    # T5 past the horizon: literal continuation leaves the class and its memory factor is the reciprocal of the mirror's
    beyond = [1/x for x in tower[::-1]]
    if max(rb.cells_of_profile(beyond)) >= 0:
        raise ValueError('cells beyond the horizon should be negative')
    if returns_along(recips_of_profile(beyond), beyond[-1])[0] != beyond[0]:
        raise ValueError('literal continuation is not a chain return')
    mirror = tower[::-1]
    sb, sm = sensitivity(beyond), sensitivity(mirror)
    if sb*sm != 1 or not (sm < 1 < sb) or min(rb.cells_of_profile(mirror)) <= 0:
        raise ValueError('mirror relation')
    res['past_horizon'] = {'literal': float(sb), 'mirror': float(sm)}
    return res


if __name__ == '__main__':
    out = run()
    with open('HM1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
