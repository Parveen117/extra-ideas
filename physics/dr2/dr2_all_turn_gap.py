"""DR2: closed radial floors and an explicit physical core gap for every d >= 3.

The note supplies the analytical proof (nodal form restriction, min-max,
row-parity inclusion, Gaussian integration, ground-state symmetry).
This stdlib certificate checks exact polynomial identities in symbolic d,
replays DR1's frozen small-d sign certificates, and rounds only outwards.
No float search, fitted dimension formula, or approximate deflation is used.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parent
spec = importlib.util.spec_from_file_location('dr2_dr1', PHYSICS/'dr1/dr1_core_by_rows.py')
dr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dr)


# Q[d], stored in ascending powers. Empty tuple is zero.
def poly(*xs):
    xs = list(map(F, xs))
    while xs and not xs[-1]:
        xs.pop()
    return tuple(xs)


def add(p, q):
    return poly(*( (p[i] if i < len(p) else 0)+(q[i] if i < len(q) else 0)
                   for i in range(max(len(p), len(q))) ))


def scale(p, a):
    return poly(*(a*x for x in p))


def mul(p, q):
    if not p or not q:
        return ()
    out = [F(0)]*(len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return poly(*out)


def evaluate(p, d):
    out = F(0)
    for x in reversed(p):
        out = d*out+x
    return out


def shift(p, start):
    return poly(*(sum(p[k]*comb(k, j)*start**(k-j) for k in range(j, len(p)))
                  for j in range(len(p))))


D, N = poly(0, 1), poly(-1, 1)
ZERO, E1, E2 = (0, 0, 0), (1, 0, 0), (0, 1, 0)


# Q[d][e1,e2,e3]. The degree lowering generator is CR1's exact L = Delta/2.
def degree(e):
    return e[0]+2*e[1]+3*e[2]


def p_add(p, q):
    out = dict(p)
    for e, v in q.items():
        out[e] = add(out.get(e, ()), v)
    return {e: v for e, v in out.items() if v}


def p_scale(p, q):
    return {e: v for e, c in p.items() if (v := mul(c, q))}


def p_mul(p, q):
    out = {}
    for a, x in p.items():
        for b, y in q.items():
            e = tuple(i+j for i, j in zip(a, b))
            out[e] = add(out.get(e, ()), mul(x, y))
    return {e: v for e, v in out.items() if v}


def gen_mono(e):
    a, b, c = e
    entries = [((a-1, b, c), poly(2*a*(a-1)+8*a*b+12*a*c, 3*a)),
               ((a+1, b-1, c), poly(-2*b+2*b*(b-1)+8*b*c, 2*b)),
               ((a, b+1, c-1), poly(-2*c+2*c*(c-1), c)),
               ((a, b-2, c+1), poly(6*b*(b-1)))]
    return {k: v for k, v in entries if v}


@lru_cache(None)
def mono_mean(e):
    if e == ZERO:
        return poly(1)
    out = ()
    for k, v in gen_mono(e).items():
        out = add(out, mul(v, mono_mean(k)))
    return scale(out, F(1, 2*degree(e)))


def mean(p):
    out = ()
    for e, v in p.items():
        out = add(out, mul(v, mono_mean(e)))
    return out


def scaled_action(p):
    """A = (d-1) hbar acting on polynomial times exp(-e1/2)."""
    out = {}
    for e, v in p.items():
        out = p_add(out, {k: scale(mul(N, mul(v, c)), -1) for k, c in gen_mono(e).items()})
        out = p_add(out, {e: mul(N, mul(v, poly(2*degree(e), F(3, 2))))})
    out = p_add(out, p_mul(p, {E1: scale(N, -F(1, 2)), E2: poly(F(1, 2))}))
    return out


def source_numerator():
    """(d-1) r; r is the actual Gaussian residual, not an inferred hidden mode."""
    return {ZERO: scale(mul(D, N), F(3, 8)), E1: scale(N, -F(1, 2)), E2: poly(F(1, 2))}


def integer_at_least(x, low, name):
    if isinstance(x, bool) or not isinstance(x, int) or x < low:
        raise ValueError(f'{name} must be an integer >= {low}')


def laguerre(m, j):
    integer_at_least(m, 3, 'm')
    integer_at_least(j, 0, 'j')
    alpha = F(2*m-4, 3)
    cs = [F(1)]
    for k in range(j):
        cs.append(-(j-k)*cs[-1]/((k+1)*(alpha+k+1)))
    return tuple(cs)


def laguerre_residual(cs, alpha, j):
    """Coefficients of x q'' + (alpha+1-x) q' + j q."""
    return poly(*((k+1)*(k+alpha+1)*(cs[k+1] if k+1 < len(cs) else 0)+(j-k)*cs[k]
                  for k in range(len(cs))))


def radial_floor_cube(m, j):
    integer_at_least(m, 3, 'm')
    integer_at_least(j, 0, 'j')
    return F(27, 64)*(2*m-1+6*j)**2


def root_floor(x, digits=12):
    x = F(x)
    integer_at_least(digits, 0, 'digits')
    if x < 0:
        raise ValueError('nonnegative radicand required')
    unit = 10**digits
    target = x.numerator*unit**3//x.denominator
    lo, hi = 0, 1
    while hi**3 <= target:
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**3 <= target:
            lo = mid
        else:
            hi = mid
    return F(lo, unit)


def root_upper(x, digits=12):
    lo = root_floor(x, digits)
    return lo if lo**3 == x else lo+F(1, 10**digits)


def down(x, digits=10):
    unit = 10**digits
    return F((x*unit).numerator//(x*unit).denominator, unit)


def up(x, digits=10):
    return -down(-x, digits)


def gain(d):
    integer_at_least(d, 3, 'd')
    return F(3*d*(23*d-35), 2*(d-1)*(521*d-512))


def gap_coefficient(d):
    integer_at_least(d, 9, 'd')
    return gain(d)-F(17, 32*(d-1))


def ground_upper_scaled(d):
    return F(9*d, 8)-gain(d)


def excitation_upper_scaled(d):
    integer_at_least(d, 3, 'd')
    return F(9*d, 8)+2-F(1, 2*(d-1))


def gap_upper_coefficient(d):
    integer_at_least(d, 3, 'd')
    return F(11, 4)-F(15, 32*(d-1))


def frozen_floors():
    result = json.loads((PHYSICS/'dr1/DR1_RESULT.json').read_text(), parse_float=F)
    rows = result['numbers']['floors of one row: m, first level, second level (l = 0)']
    return {m: (F(a), F(b) if b is not None else None) for m, a, b in rows}


def small_count_scaled(d, floors):
    integer_at_least(d, 3, 'd')
    if d > 10:
        raise ValueError('frozen sign floors cover only d <= 10')
    a, b = floors[d]
    value = min(2*a+min(b, floors[d+4][0]), 3*floors[d+2][0])
    return root_floor(F(d-1, 32))*value


def run():
    checks = {}
    source = source_numerator()
    eta = scale(D, F(9, 8))
    bracket = p_add(scaled_action(source), p_scale(source, scale(mul(N, eta), -1)))
    checks['actual Gaussian source: n(hbar-eta)1 = N'] = (
        p_add(scaled_action({ZERO: poly(1)}), {ZERO: scale(mul(N, eta), -1)}) == source)
    checks['symbolic Gaussian mean of the source is zero'] = mean(source) == ()
    checks['symbolic source variance numerator is 9d(d-1)/32'] = (
        mean(p_mul(source, source)) == scale(mul(D, N), F(9, 32)))
    checks['symbolic source energy numerator is 3d(d-1)(25d-13)/64'] = (
        mean(p_mul(source, bracket)) == scale(mul(mul(D, N), poly(-13, 25)), F(3, 64)))
    # The exact corrected Rayleigh gain, clearing its positive denominators.
    vnum = scale(D, F(9, 32))       # v = vnum/n
    bnum = scale(mul(D, poly(-13, 25)), F(3, 64))  # b = bnum/n^2
    gainnum = add(scale(mul(vnum, N), F(1, 2)), scale(bnum, -F(1, 16)))
    gainden = mul(N, add(N, scale(vnum, F(1, 16))))
    gnum, gden = scale(mul(D, poly(-35, 23)), 3), scale(mul(N, poly(-512, 521)), 2)
    checks['one-quarter actual-source correction has the claimed rational gain'] = (
        mul(gainnum, gden) == mul(gnum, gainden))
    checks['gain numerator and denominator are positive for every d >= 3'] = all(
        all(c > 0 for c in shift(p, 3)) for p in (gnum, gden))

    cases = [(m, j) for m in (3, 4, 7, 14, 101) for j in range(13)]
    checks['Laguerre coefficient recurrence solves the exact differential equation'] = all(
        not laguerre_residual(laguerre(m, j), F(2*m-4, 3), j) for m, j in cases)
    checks['Laguerre change-of-variable derivative terms cancel exactly'] = all(
        F(9, 2)*(F(2*m-4, 3)+1) == 6*F(m-1, 2)+F(3, 2) for m, _ in cases)
    checks['minimum of A/sqrt(s)+s/2 has the claimed cubed floor'] = all(
        F(27, 8)*F(9, 16)*F(2, 9)*(2*m-1+6*j)**2 == radial_floor_cube(m, j)
        for m, j in cases)
    # Taylor lower polynomial 1 + 2x/3 - x^2/9; third derivative is positive.
    # Listed coefficients are exact after x=a/n substitution, not numerical fits.
    radial_correction = 2*F(1, 2)**2+F(7, 2)**2
    checks['radial Taylor loss is exactly 17/(32n)'] = F(3, 8)*radial_correction/9 == F(17, 32)
    checks['odd-branch Taylor reserve is positive for all n >= 2'] = F(3, 4)-F(25, 64) > 0
    numerator = poly(8704, -10537, 1104)
    den = scale(mul(N, poly(-512, 521)), 32)
    checks['gap coefficient equals its stated rational polynomial'] = (
        add(scale(gnum, 16), scale(poly(-512, 521), -17)) == numerator and den == scale(gden, 16))
    checks['gap coefficient is positive for all d >= 9 by shifted coefficients'] = shift(numerator, 9) == poly(3295, 9335, 1104)
    uniform = add(scale(numerator, 2), scale(mul(N, poly(-512, 521)), -1))
    checks['coefficient exceeds 1/64 for all d >= 11 by shifted coefficients'] = shift(uniform, 11) == poly(572, 17073, 1687)

    floors = frozen_floors()
    checks['frozen d=3 Airy first and second floors replay with exact signs'] = (
        dr.gc1.level_one(floors[3][0]) and dr.gc1.level_two(floors[3][1])[0])
    checks['all frozen first-row floors m=4..14 replay with exact signs'] = all(
        dr.level_one(a, m-1, int(a)+7) for m, (a, _) in floors.items() if m >= 4)
    checks['all frozen second-row floors d=4..10 replay with exact signs'] = all(
        dr.level_two(b, m-1, int(b)+7) for m, (_, b) in floors.items() if m >= 4 and b is not None)
    small = {d: small_count_scaled(d, floors)-ground_upper_scaled(d) for d in range(3, 11)}
    checks['each small-d normalized gap exceeds 1/64 with the explicit source trial'] = all(v > F(1, 64) for v in small.values())

    qnorm = F(2)*(F(15, 4)-F(9, 4))
    xq2 = scale(add(scale(poly(-3, 2), F(45, 4)),
                   scale(mul(poly(-2, 1), poly(-3, 1)), F(27, 8))), F(2, 3))
    checks['turning trial has norm 3 and exact symbolic potential energy'] = (
        qnorm == 3 and xq2 == scale(poly(-12, 5, 3), F(3, 4)))
    # n*(T + V/(2n)) = n*(9d/8+2)-1/2.
    checks['turning trial energy is 9d/8+2-1/(2n) in w units'] = (
        add(mul(N, poly(1, F(3, 4))), scale(xq2, F(1, 6))) ==
        add(mul(N, poly(2, F(9, 8))), poly(-F(1, 2))))
    checks['upper gap coefficient follows from the ground Taylor floor'] = (
        F(2)+F(3, 4) == F(11, 4) and -F(1, 2)+F(1, 32) == -F(15, 32))
    checks['large-d lower coefficient is exactly 69/1042 before the 2^(1/3) factor'] = F(1104, 32*521) == F(69, 1042)

    rows = []
    for d in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 20, 100, 1000):
        lower = max(F(1, 64), gap_coefficient(d) if d >= 9 else F(0), small.get(d, F(0)))
        wl, wu = root_floor(2*(d-1)), root_upper(2*(d-1))
        rows.append({'d': d, 'gap_lower': str(down(wl*lower)),
                     'gap_upper': str(up(wu*gap_upper_coefficient(d))),
                     'source_trial_E0_upper': str(up(wu*ground_upper_scaled(d)))})
    return {'stage': 'DR2', 'checks': checks, 'all_pass': all(checks.values()),
            'universal_gap_over_w_lower': '1/64',
            'small_d_gap_over_w_lower': {str(d): str(down(v)) for d, v in small.items()},
            'bounds_using_DR2_trials': rows,
            'liminf_gap_over_d_cuberoot_lower': str(down(root_floor(2)*F(69, 1042))),
            'limsup_gap_over_d_cuberoot_upper': str(up(root_upper(2)*F(11, 4))),
            'claim_boundary': {'carrier': 'L2(R^(3xd),dC), SO(3) row-gauge singlets',
                'dimension': 'integer d >= 3; column count, not volume',
                'all_column_direction_sectors_in_lower_bound': True,
                'excited_trial_orthogonality': 'exact column-swap parity against the unique positive ground',
                'lowest_excited_sector_identified': False, 'exact_gap_coefficient': False,
                'volume_uniform_YM_gap': False, 'continuum_mass_gap': False,
                'proof_status': 'written analytic proof with exact arithmetic replay, not formal proof assistant verification'}}


def source_pins():
    paths = ['dr1/DR1_CORE_BY_ROWS.md', 'dr1/dr1_core_by_rows.py', 'dr1/DR1_RESULT.json',
             'gc1/gc1_core_gap_count.py', 'cm1/cm1_core_sector_memory.py', 'cr1/cr1_core_rates.py',
             'tc1/tc1_core_of_turns.py', 'dr2/dr2_all_turn_gap.py', 'dr2/DR2_ALL_TURN_GAP.md', 'dr2/test_dr2.py']
    return {'physics/'+p: hashlib.sha256((PHYSICS/p).read_bytes()).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if not result['all_pass']:
        raise SystemExit('FAIL: '+', '.join(k for k, v in result['checks'].items() if not v))
    result['source_sha256'] = source_pins()
    path = HERE/'DR2_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale DR2 result or source pin')
    print('PASS', len(result['checks']), 'exact/outward checks')
    print('Every integer d >= 3: w/64 <= gap <= w*(11/4-15/(32*(d-1))), w^3=2*(d-1).')
    print('Asymptotic lower / upper:', result['liminf_gap_over_d_cuberoot_lower'], '/', result['limsup_gap_over_d_cuberoot_upper'])


if __name__ == '__main__':
    main()
