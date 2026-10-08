"""TS1: a second sheet with its own flat end removes the closure freedom of the horizon end.

Stdlib only, exact rationals.  Builds on RB1 (same line) and on
  Recognition-Kernel-Framework research/recognition_return  R1 RC1, RC2 and r2 RI1, RI2, SY2 (unchanged).
"""
from fractions import Fraction as F
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / 'rb1'))
import rb1_return_boost_boundary_memory as rb  # noqa: E402


def inward(x0, levels):
    xs = [F(x0)]
    for _ in range(levels):
        xs.append(rb.tower(xs[-1]))
    return xs


def flat_end(n, start=4):
    return [F(1, j + start) for j in range(n)]


def two_sheet_profile(x0, levels, n_out):
    xs = inward(x0, levels)
    return xs + xs[::-1] + flat_end(n_out)


def run():
    res = {}
    x0 = F(1, 3)
    # T1 the tower does not distinguish x from 1/x; in radius: r -> r_s^2/r (m = x^2 -> 1/m)
    for x in (F(1, 3), F(2, 7), F(9, 10)):
        if rb.tower(x) != rb.tower(1/x):
            raise ValueError('tower is not mirror symmetric')
        m = x*x
        if (1 + m)**2/(4*m) != (1 + 1/m)**2/(4/m):      # LT1 radius map (r + r_s)^2/4r in units r_s = 1, r = 1/m
            raise ValueError('radius map is not inversion symmetric')

    # T2 one sheet: freedom stays; two sheets with a flat end: it goes
    xs = inward(x0, 6)
    one = rb.interval(rb.cells_of_profile(xs))[2]
    rows = []
    for n_out in (0, 5, 20, 80):
        cells = rb.cells_of_profile(two_sheet_profile(x0, 6, n_out))
        if min(cells) <= 0:
            raise ValueError('cell not positive')
        e0, einf, w, _, _ = rb.interval(cells)
        if not (min(e0, einf) <= x0 <= max(e0, einf)):
            raise ValueError('point-mass value outside the interval')
        rows.append((n_out, w))
    if not all(rows[i+1][1] < rows[i][1] for i in range(len(rows) - 1)):
        raise ValueError('second sheet does not close the interval')
    if not (one > F(12, 100) and rows[0][1] < F(3, 100) and rows[-1][1] < F(1, 10**200)):
        raise ValueError('closure verdict fails')
    res['width_one_sheet'] = float(one)
    res['width_two_sheets'] = {str(n): float(w) for n, w in rows}

    # T3 the throat cell is R1's uniform coupling; its return is the selected cut up to d/(4q)
    xt = xs[-1]
    cells = rb.cells_of_profile(xs + xs[::-1])
    q = cells[len(xs) - 1]
    if q != xt/(1 - xt*xt):
        raise ValueError('throat cell is not the R1 coupling')
    Q = rb.el(F(1, 2), xt/2)
    defect = rb.add(rb.mul(Q, Q), rb.scal(-1, Q))
    if defect != rb.scal(-xt/(4*q), rb.ONE):
        raise ValueError('idempotent defect law fails')
    res['throat'] = {'coupling_digits': len(str(q.numerator // q.denominator)),
                     'idempotent_defect_below': float(xt/(4*q))}

    # T4 RC2 on the chain: strong throat alone (no depth behind it) does not select the cut
    alone = rb.chain_return([q], F(0))          # one strong cell, wall behind: returns q itself, not 1
    behind = rb.chain_return([q] + rb.cells_of_profile(xs[::-1] + flat_end(40)), F(0))
    if not (alone == q and abs(behind - xt) < F(1, 10**30)):
        raise ValueError('joint-limit statement fails')
    res['strong_cell_alone'] = 'returns q (unbounded)'
    res['strong_cell_with_sheet_behind'] = 'returns the throat value'
    return res


if __name__ == '__main__':
    out = run()
    with open('TS1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
