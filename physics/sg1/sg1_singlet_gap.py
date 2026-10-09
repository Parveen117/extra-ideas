"""SG1: outward radial counts and a gap on the physical gauge-singlet core.

Proof arithmetic is Fraction only. The shooting equation is y''=(r-E)y,
y(0)=0, y'(0)=1; no Airy values or fitted floating constants are inputs.
"""
from fractions import Fraction as F
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parent


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cr = load(PHYSICS/'cr1/cr1_core_rates.py', 'sg1_cr1')


def iv(x):
    return F(x), F(x)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def scale(a, c):
    return mul(a, iv(c))


def sign(a):
    if a[0] > 0:
        return 1
    if a[1] < 0:
        return -1
    return 0


def ceil(x):
    return -((-x.numerator)//x.denominator)


def round_out(a, bits):
    grid = 2**bits
    return F((a[0]*grid).numerator//(a[0]*grid).denominator, grid), F(ceil(a[1]*grid), grid)


def horner(coefficients, t):
    out = iv(0)
    for c in reversed(coefficients):
        out = add(c, mul(t, out))
    return out


def taylor_cell(a, energy, initial, step=F(1, 8), order=32, bits=160):
    """Enclose endpoint and whole cell, using a complex-disk Cauchy tail.

    On |t|<=1, ||[[0,1],[a+t-E,0]]||_infty <= M=max(1,|a-E|+1).
    Thus ||U(t)|| <= exp(M)||U(0)|| < 4^ceil(M)||U(0)|| = B.
    Each omitted vector coefficient is <=B; tail <=B*h^(N+1)/(1-h).
    """
    if not (0 < step < 1) or order < 2:
        raise ValueError('Taylor disk requires 0 < step < 1 and order >= 2')
    ys = [initial[0], initial[1]]
    for n in range(order):
        prior = ys[n-1] if n else iv(0)
        ys.append(scale(add(scale(ys[n], a-energy), prior), F(1, (n+2)*(n+1))))
    coeffs = [ys[:order+1], [scale(ys[n+1], n+1) for n in range(order+1)]]
    norm = max(abs(x) for component in initial for x in component)
    majorant = norm*4**ceil(max(F(1), abs(a-energy)+1))
    error = majorant*step**(order+1)/(1-step)
    remainder = (-error, error)
    endpoint = tuple(round_out(add(horner(p, iv(step)), remainder), bits) for p in coeffs)
    ranges = tuple(add(horner(p, (F(0), step)), remainder) for p in coeffs)
    return endpoint, ranges, error


@lru_cache(None)
def radial_count(energy, radius=F(12), step=F(1, 8), order=32, bits=160):
    """Certified Dirichlet-at-zero/Neumann-at-radius eigenvalue count below E.

    Sturm oscillation: count = interior zeros + [y(R)y'(R)<0].
    Every cell either excludes zero or has a strictly signed derivative.
    The initial Dirichlet zero is not counted. Unresolved signs fail closed.
    """
    if radius <= energy:
        raise ValueError('tail potential floor must exceed the test energy')
    cells = radius/step
    if cells.denominator != 1:
        raise ValueError('radius must contain an integer number of cells')
    current = (iv(0), iv(1))
    zeros = []
    methods = {'zero_excluded': 0, 'monotone': 0}
    biggest_error = F(0)
    for cell in range(int(cells)):
        a = cell*step
        nxt, ranges, error = taylor_cell(a, energy, current, step, order, bits)
        biggest_error = max(biggest_error, error)
        left, right = sign(current[0]), sign(nxt[0])
        if not right or (cell and not left):
            raise ValueError(f'unresolved endpoint sign at cell {cell}')
        if sign(ranges[0]):
            methods['zero_excluded'] += 1
            if left != right:
                raise ValueError('whole-cell sign contradicts endpoints')
        elif sign(ranges[1]):
            methods['monotone'] += 1
            direction = sign(ranges[1])
            if cell == 0:
                if current[0] != iv(0) or right != direction:
                    raise ValueError('initial Dirichlet endpoint not resolved')
            elif left != right:
                if right != direction:
                    raise ValueError('crossing contradicts derivative direction')
                zeros.append([str(a), str(a+step)])
        else:
            raise ValueError(f'unresolved zero count in cell {cell}')
        current = nxt
    end_signs = [sign(component) for component in current]
    if not all(end_signs):
        raise ValueError('endpoint is an unresolved eigenvalue')
    count = len(zeros)+int(end_signs[0] != end_signs[1])
    return {'energy': str(energy), 'radius': str(radius), 'step': str(step),
            'order': order, 'rounding_bits': bits, 'cells': int(cells),
            'interior_zero_cells': zeros, 'endpoint_signs': end_signs,
            'endpoint_intervals': [[str(v) for v in component] for component in current],
            'cell_certification': methods, 'largest_local_remainder': str(biggest_error),
            'negative_count': count, 'dirichlet_negative_count': len(zeros),
            'tail_floor': str(radius)}


def coupled_spins(labels):
    """Exact SO(3) tensor-product support, for integer orbital spins."""
    support = {0}
    for ell in labels:
        support = {j for old in support for j in range(abs(old-ell), old+ell+1)}
    return support


def run():
    checks = {}
    a0, a1, b = F(23381, 10000), F(40879, 10000), F(649, 200)
    radial_upper = F(11691, 5000)
    radial = [radial_count(e) for e in (a0, a1, radial_upper)]
    checks['radial ground has no eigenvalue below 2.3381'] = radial[0]['negative_count'] == 0
    checks['radial second has at most one eigenvalue below 4.0879'] = radial[1]['negative_count'] == 1
    checks['Dirichlet trial certifies radial ground below 2.3382'] = radial[2]['dirichlet_negative_count'] >= 1
    checks['all radial cells resolved including the Neumann endpoint'] = all(
        sum(row['cell_certification'].values()) == row['cells'] for row in radial)
    checks['l>=1 positive-test-function bound exceeds 3.245'] = b**3 < F(2187, 64)
    checks['l>=1 lower is above l=0 ground lower'] = b > a0
    checks['one nonzero orbital spin contains no gauge singlet'] = all(
        0 not in coupled_spins([ell, 0, 0]) for ell in range(1, 9))
    checks['two spin-one factors admit a singlet'] = 0 in coupled_spins([1, 1, 0])
    numbers = {}
    for d, c, upper, lower in [(3, F(15749, 25000), F(12967, 2500), F(5521, 1000)),
                               (4, F(180281, 250000), F(8003, 1000), F(4003, 500))]:
        checks[f'outward cube-root scaling d={d}'] = c**3 < F(d-1, 8)
        radial_excited = c*((d-1)*a0+a1)
        angular_excited = c*((d-2)*a0+2*b)
        bound = min(radial_excited, angular_excited)
        checks[f'all singlet alternatives exceed second-rate lower d={d}'] = bound > lower
        keys, s, h = cr.ladder(d, cr.RATE[d], 10)
        count = cr.below(s, h, upper)
        checks[f'CR1 actual scalar trial below ground upper d={d}'] = count >= 1
        checks[f'positive certified physical-core gap d={d}'] = lower > upper
        c_upper = c+F(1, 1000000)
        trial_a = F(27, 25)
        ell1_upper = trial_a**2+F(5, 2)/trial_a
        false_low_mode = c_upper*((d-1)*radial_upper+ell1_upper)
        checks[f'ungauged one-spin comparison excitation lies below ground upper d={d}'] = (
            c_upper**3 > F(d-1, 8) and false_low_mode < upper)
        numbers[str(d)] = {'scale_lower': str(c), 'radial_ground_lower': str(a0),
            'radial_second_lower': str(a1), 'nonzero_spin_ground_lower': str(b),
            'core_ground_lower': str(c*d*a0), 'radial_excited_comparison_lower': str(radial_excited),
            'angular_excited_comparison_lower': str(angular_excited),
            'core_second_lower': str(lower), 'core_ground_upper': str(upper),
            'core_gap_lower': str(lower-upper), 'upper_trial_basis_size': len(keys),
            'upper_trial_negative_count': count,
            'excluded_one_spin_comparison_upper': str(false_low_mode)}
    return {'schema': 'sg1_singlet_gap_v1', 'arithmetic': 'Fraction with outward dyadic rounding',
            'checks': checks, 'radial_certificates': radial, 'numbers': numbers,
            'claim_boundary': {'free_bosonic_SO3_gauge_singlet_core_gap': 'PROVED_WITH_DECLARED_ANALYTIC_INPUTS',
                'all_direction_sectors': 'COVERED_BY_GAUGE_SINGLET_COMPARISON',
                'lowest_excited_sector_identity': 'OPEN', 'continuum_mass_gap': 'OPEN',
                'TVSP_to_VTSP_physical_identification': 'OWNER_RESEARCH_INTERPRETATION_NOT_A_GAP_ASSUMPTION'}}


def source_pins():
    paths = [PHYSICS/'cr1/cr1_core_rates.py', PHYSICS/'cr1/CR1_CORE_RATES.md',
             PHYSICS/'tc1/tc1_core_of_turns.py', PHYSICS/'cm1/CM1_CORE_SECTOR_MEMORY.md',
             HERE/'sg1_singlet_gap.py', HERE/'SG1_GAUGE_SINGLET_GAP.md', HERE/'test_sg1.py']
    return {str(p.relative_to(PHYSICS.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if not all(result['checks'].values()):
        raise SystemExit('FAIL: '+', '.join(k for k, v in result['checks'].items() if not v))
    result['source_sha256'] = source_pins()
    path = HERE/'SG1_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale SG1 result or source pin')
    print('PASS', len(result['checks']), 'exact/outward checks')
    print(json.dumps(result['numbers'], indent=2))


if __name__ == '__main__':
    main()
