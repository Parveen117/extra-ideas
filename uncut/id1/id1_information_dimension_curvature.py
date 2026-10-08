"""ID1: information, dimension and curvature.

The owner's statement: gravity may be a coupling of information and curvature.  On a straight line the
dimension is one and the information is unbounded; as the dimension grows the information is distributed
among the phases (directions); in three dimensions it is distributed further.

Sources read (unchanged): extra-ideas physics mc1 T1 (least-cost spreading: r^-(d-2), log r, r), cv1 (what no frame change removes:
  r_s/r^3 x (1, -1/2, -1/2); contracted curvature zero iff m = r_s/r), mo1 (fall speed^2 = memory m), in1 T2 (a single reading:
  n^2 = sum r_i^2), gb1; uncut ci1 (exp(-2 I) = 1 - m), up1 T6 (the far end is the flat special case).
No measurement.  sympy, exact.  Python 3.12.
"""
import json

import sympy as sp

r = sp.Symbol('r', positive=True)


def coords(d):
    return sp.symbols(f'x1:{d + 1}', real=True)


def radius(xs):
    return sp.sqrt(sum(x*x for x in xs))


def spread(d, rad):
    """the memory of a point source in d space dimensions (MC1-T1), unit strength"""
    if d == 1:
        return rad
    if d == 2:
        return sp.log(rad)
    return rad**(-(d - 2))


def hessian_on_axis(f, xs):
    """Hessian of f at the point (r, 0, ..., 0)"""
    H = sp.hessian(f, xs)
    at = {xs[0]: r, **{x: 0 for x in xs[1:]}}
    return H.subs(at).applyfunc(sp.simplify)


def run():
    res = {}

    # T1 one cut among D holds on average 1/D of a single reading
    for D in (1, 2, 3, 4):
        rs = sp.symbols(f'r1:{D + 1}', real=True)
        n2 = sum(x*x for x in rs)
        if sp.simplify(sum(x*x/n2 for x in rs)/D - sp.Rational(1, D)) != 0:
            raise ValueError('share of one cut')
    res['share_of_one_cut'] = {'1': '1', '2': '1/2', '3': '1/3', 'D': '1/D'}

    # T2 spreading: flux through every shell the same; the far value is finite only from d = 3 on
    table = {}
    for d in (1, 2, 3, 4, 5):
        m = spread(d, r)
        flux = sp.simplify(r**(d - 1)*sp.diff(m, r))
        if flux.has(r):
            raise ValueError('flux is not constant')
        lap = sp.simplify(sp.diff(m, r, 2) + (d - 1)/r*sp.diff(m, r))
        if lap != 0:
            raise ValueError('not harmonic')
        far = sp.limit(m, r, sp.oo)
        table[d] = {'law': str(m), 'flux': str(flux), 'far_value': str(far)}
    if not (table[1]['far_value'] == 'oo' and table[2]['far_value'] == 'oo' and table[3]['far_value'] == '0'):
        raise ValueError('far values')
    res['spreading'] = table

    # T3 curvature = how the spread is shared among directions: one radial direction against d - 1 transverse ones
    tid = {}
    for d in (1, 2, 3, 4, 5):
        xs = coords(d)
        H = hessian_on_axis(spread(d, radius(xs)), xs)
        if not H.is_diagonal():
            raise ValueError('axis frame')
        diag = [sp.simplify(H[j, j]) for j in range(d)]
        if sp.simplify(sum(diag)) != 0:
            raise ValueError('trace')
        if d == 1:
            if diag[0] != 0:
                raise ValueError('line')
            tid[d] = 'none: a line has no direction to share with'
            continue
        if any(sp.simplify(q - diag[1]) != 0 for q in diag[1:]):
            raise ValueError('transverse')
        ratio = sp.simplify(diag[1]/diag[0])
        if ratio != -sp.Rational(1, d - 1):
            raise ValueError('ratio')
        tid[d] = f'radial : each transverse = 1 : {ratio}'
    res['sharing_among_directions'] = tid
    rs_ = sp.Symbol('r_s', positive=True)
    xs = coords(3)
    H3 = hessian_on_axis(rs_/radius(xs), xs)/2
    if [sp.simplify(H3[j, j]*r**3/rs_) for j in range(3)] != [1, -sp.Rational(1, 2), -sp.Rational(1, 2)]:
        raise ValueError('CV1 tidal part')
    res['cv1'] = 'half the Hessian of m = r_s/r is r_s/r^3 x (1, -1/2, -1/2): the part of the curvature no frame change removes'

    # T4 in terms of I (exp(-2I) = 1 - m): the free law is  Laplacian(I) = 2 |grad I|^2
    xs = coords(3)
    rad = radius(xs)
    I = -sp.log(1 - rs_/rad)/2
    lapI = sum(sp.diff(I, x, 2) for x in xs)
    gradsq = sum(sp.diff(I, x)**2 for x in xs)
    if sp.simplify(lapI - 2*gradsq) != 0:
        raise ValueError('law for I')
    Ig = sp.Function('I')(*xs)
    mg = 1 - sp.exp(-2*Ig)
    for a_ in range(3):
        for b_ in range(3):
            lhs = sp.diff(mg, xs[a_], xs[b_])
            rhs = 2*sp.exp(-2*Ig)*(sp.diff(Ig, xs[a_], xs[b_]) - 2*sp.diff(Ig, xs[a_])*sp.diff(Ig, xs[b_]))
            if sp.simplify(lhs - rhs) != 0:
                raise ValueError('Hessian in terms of I')
    res['law_for_I'] = {'free': 'Laplacian(I) = 2 |grad I|^2', 'curvature': '(1/2) Hessian(m) = exp(-2I) ( Hessian(I) - 2 grad I grad I )',
                        'far': 'Hessian(I), with Laplacian(I) = 0 to first order'}
    return res


if __name__ == '__main__':
    out = run()
    with open('ID1_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
