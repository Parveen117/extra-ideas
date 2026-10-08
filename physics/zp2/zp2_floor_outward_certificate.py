"""ZP2: exact finite certificates for the two floor constants of ZP1, and the hypothesis they need.

Stdlib only, exact rationals.  Sources read (unchanged):
  Recognition-Kernel-Framework theorum/28  section 4 (finite packet + declared correction + vanishing tail)
                                           section 9 (a top without an explicit e_n is a shadow certificate)
  this line: ZP1 (Z3 one cut, Z4 three cuts), PS1-T4, LC1-C4 / WQ1-W3 (floor of a mode), RB1 (edge memory).
Cut-off: each mode's floor is weighted by eta(w/W), eta(t) = (1-t)^m on [0,1], y = W L / (pi c) an integer.
"""
from fractions import Fraction as F
import json


def line_sum(m, y):
    """one cut: sum over n of n (1 - n/y)^m"""
    return sum(n * (1 - F(n, y))**m for n in range(y + 1))


def g(m, u):
    """integral from u to 1 of t^2 (1-t)^m dt, exact (binomial expansion)"""
    u = F(u)
    c, total = F(1), F(0)
    for k in range(m + 1):            # (1-t)^m = sum (-1)^k C(m,k) t^k
        total += (-1)**k * c * (1 - u**(k + 3)) / (k + 3)
        c = c * (m - k) / (k + 1)
    return total


def plate_sum(m, y):
    """three cuts: g(0)/2 + sum over n >= 1 of g(n/y)"""
    return g(m, 0)/2 + sum(g(m, F(n, y)) for n in range(1, y + 1))


def interpolate(points):
    """exact coefficients (low to high) of the polynomial through the points"""
    n = len(points)
    coef = [F(0)]*n
    for i, (xi, yi) in enumerate(points):
        basis, den = [F(1)], F(1)
        for j, (xj, _) in enumerate(points):
            if j != i:
                basis = [F(0)] + basis
                for k in range(len(basis) - 1):
                    basis[k] -= xj*basis[k + 1]
                den *= xi - xj
        for k in range(n):
            coef[k] += yi*basis[k]/den
    return coef


def laurent(fn, m, shift, degree, extra=6):
    """y^shift * fn(m, y) is a polynomial of the given degree; return {power of y: coefficient}, checked on extra points."""
    ys = list(range(1, degree + 2))
    coef = interpolate([(F(y), fn(m, y)*F(y)**shift) for y in ys])
    for y in range(degree + 2, degree + 2 + extra):
        if sum(c*F(y)**k for k, c in enumerate(coef)) != fn(m, y)*F(y)**shift:
            raise ValueError('not a polynomial of the declared degree')
    return {k - shift: c for k, c in enumerate(coef) if c}


def run():
    res = {'line': {}, 'plates': {}}
    # one cut
    for m in range(0, 8):
        L = laurent(line_sum, m, m, m + 2)
        res['line'][m] = {str(k): str(v) for k, v in sorted(L.items(), reverse=True)}
        const = L.get(0, F(0))
        want = {0: F(0), 1: F(-1, 6)}.get(m, F(-1, 12))
        if const != want:
            raise ValueError('line constant')
        if L.get(2) != F(1, (m + 1)*(m + 2)):
            raise ValueError('bulk term')
        if m >= 1 and L.get(1, 0) != 0:
            raise ValueError('unexpected length-independent term')
    if laurent(line_sum, 2, 2, 4) != {2: F(1, 12), 0: F(-1, 12)}:
        raise ValueError('m = 2 is not exact')
    # explicit tail for m >= 2: |F - bulk + 1/12| <= e_m / y^2
    tails = {}
    for m in range(2, 8):
        L = laurent(line_sum, m, m, m + 2)
        e = sum(abs(v) for k, v in L.items() if k < 0)
        for y in (1, 2, 5, 40):
            if abs(line_sum(m, y) - L[2]*y*y + F(1, 12)) > e/F(y*y):
                raise ValueError('tail bound fails')
        tails[m] = e
    res['line_tail_e_m'] = {str(m): str(e) for m, e in tails.items()}

    # three cuts
    for m in range(0, 8):
        P = laurent(plate_sum, m, m + 3, m + 4)
        res['plates'][m] = {str(k): str(v) for k, v in sorted(P.items(), reverse=True)}
        c3 = P.get(-3, F(0))
        if m >= 3 and c3 != F(-1, 360):
            raise ValueError('plate constant')
        if m < 3 and c3 == F(-1, 360):
            raise ValueError('constant should fail below m = 3')
        if m >= 1 and any(k in P for k in (0, -1, -2)):
            raise ValueError('unexpected term between bulk and constant')
    ptails = {}
    for m in range(3, 8):
        P = laurent(plate_sum, m, m + 3, m + 4)
        e = sum(abs(v) for k, v in P.items() if k < -3)
        for y in (1, 3, 20):
            if abs(plate_sum(m, y) - P[1]*y + F(1, 360)/F(y)**3) > e/F(y)**5:
                raise ValueError('plate tail bound fails')
        ptails[m] = e
    res['plates_tail_e_m'] = {str(m): str(e) for m, e in ptails.items()}

    # a mixed cut-off closing with zero slope: constant fixed by eta(0) = 1 alone
    mix = {2: F(3, 2), 3: F(-2), 5: F(3, 2)}                 # weights sum to 1
    y = 30
    val = sum(c*line_sum(m, y) for m, c in mix.items())
    bulk = sum(c*F(1, (m + 1)*(m + 2)) for m, c in mix.items())*y*y
    e = sum(abs(c)*tails[m] for m, c in mix.items())
    if abs(val - bulk + F(1, 12)) > e/F(y*y):
        raise ValueError('mixed cut-off fails')
    res['mixed_cutoff'] = {'deviation_from_-1/12': float(val - bulk + F(1, 12)), 'bound': float(e/F(y*y))}
    return res


if __name__ == '__main__':
    out = run()
    with open('ZP2_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
