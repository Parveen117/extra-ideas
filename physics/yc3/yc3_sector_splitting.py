"""YC3: sector-resolved compact energies with complete hidden-space errors.

All replay arithmetic is rational. Twelve variables are the four quaternion
coordinates of each of three links. The note supplies the harmonic completeness,
closed-form Schur bound and finite-window eigenvalue ordering arguments.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ZERO = (0,)*12


def padd(a, b, scale=F(1)):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, F(0))+scale*c
        if not out[e]:
            del out[e]
    return out


def pscale(a, c):
    return {e: c*x for e, x in a.items() if c*x}


def pmul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            g = tuple(x+y for x, y in zip(e, f))
            out[g] = out.get(g, F(0))+c*d
    return {e: c for e, c in out.items() if c}


def mono(*indices):
    e = [0]*12
    for i in indices:
        e[i] += 1
    return {tuple(e): F(1)}


@lru_cache(None)
def moment(e):
    if len(e) != 12 or any(not isinstance(x, int) or x < 0 for x in e):
        raise ValueError('twelve nonnegative integer powers required')
    if any(x % 2 for x in e):
        return F(0)
    answer = F(1)
    for link in range(3):
        ks = [x//2 for x in e[4*link:4*link+4]]
        answer *= F(prod(prod(range(1, 2*k, 2)) for k in ks),
                    prod(4+2*j for j in range(sum(ks))))
    return answer


def mean(p):
    return sum((c*moment(e) for e, c in p.items()), F(0))


def inner(p, q):
    return mean(pmul(p, q))


def h0(p):
    """Sum of minus spherical Laplacians, acting on ambient polynomials."""
    out = {}
    for e, c in p.items():
        for link in range(3):
            degree = sum(e[4*link:4*link+4])
            out = padd(out, {e: c*degree*(degree+2)})
            for i in range(4*link, 4*link+4):
                if e[i] >= 2:
                    f = list(e)
                    f[i] -= 2
                    out = padd(out, {tuple(f): -c*e[i]*(e[i]-1)})
    return out


def dot(i, j):
    out = {}
    for a in (1, 2, 3):
        out = padd(out, mono(4*i+a, 4*j+a))
    return out


def cross2(i, j):
    return padd(pmul(dot(i, i), dot(j, j)), pmul(dot(i, j), dot(i, j)), -1)


def potential():
    out = {}
    for i, j in combinations(range(3), 2):
        out = padd(out, cross2(i, j), 2)
    return out


def triple():
    out = {}
    for a, b, c in permutations((1, 2, 3)):
        inversions = int(a > b)+int(a > c)+int(b > c)
        out = padd(out, mono(a, 4+b, 8+c), (-1)**inversions)
    return out


def validate_k(k):
    if isinstance(k, bool) or not isinstance(k, int) or k not in range(4):
        raise ValueError('k is the number of odd links: 0, 1, 2 or 3')


def lowest_basis(k):
    validate_k(k)
    if k == 0:
        return [{ZERO: F(1)}]
    if k == 1:
        return [mono(0)]
    if k == 2:
        return [mono(0, 4), dot(0, 1)]
    return [mono(0, 4, 8), pmul(mono(0), dot(1, 2)),
            pmul(mono(4), dot(0, 2)), pmul(mono(8), dot(0, 1)), triple()]


def gauge_generator(p, a, b):
    """Simultaneous row rotation sum_i(u_ia d_ib-u_ib d_ia)."""
    out = {}
    for e, c in p.items():
        for link in range(3):
            ia, ib = 4*link+a, 4*link+b
            for i, j, sign in ((ia, ib, 1), (ib, ia, -1)):
                if e[j]:
                    f = list(e)
                    f[i] += 1
                    f[j] -= 1
                    out = padd(out, {tuple(f): sign*c*e[j]})
    return out


@lru_cache(None)
def sector_data(k):
    validate_k(k)
    basis, v = lowest_basis(k), potential()
    gram = [[inner(f, g) for g in basis] for f in basis]
    images = [pmul(v, f) for f in basis]
    matrix = [[inner(f, vg) for vg in images] for f in basis]
    slopes = [matrix[i][i]/gram[i][i] for i in range(len(basis))]
    residuals = [padd(vf, f, -a) for vf, f, a in zip(images, basis, slopes)]
    residual_gram = [[inner(f, g) for g in residuals] for f in residuals]
    residual_squares = [residual_gram[i][i]/gram[i][i] for i in range(len(basis))]
    return {'k': k, 'basis': basis, 'gram': gram, 'matrix': matrix,
            'slopes': slopes, 'residuals': residuals, 'residual_gram': residual_gram,
            'residual_squares': residual_squares, 'e': 3*k,
            'delta': 12 if k == 3 else 8,
            'a': min(slopes), 'b2': max(residual_squares)}


def lower_upper(k, theta):
    """Entire-sector floor, trial ceiling. Nonnegative theta, below denominator gate."""
    d = sector_data(k)
    theta = F(theta)
    if theta < 0 or theta*d['a'] >= d['delta']:
        raise ValueError('requires 0 <= theta and theta*a < delta')
    error = theta**2*d['b2']/(d['delta']-theta*d['a'])
    upper = d['e']+theta*d['a']
    return max(F(d['e']), upper-error), upper


def splitting(k, theta):
    if k not in (1, 2, 3):
        raise ValueError('a non-vacuum centre sector is required')
    l, u = lower_upper(k, theta)
    l0, u0 = lower_upper(0, theta)
    return l-u0, u-l0


@lru_cache(None)
def source_components(k):
    """Exact harmonic resolution of the lowest scalar-product trial's residual.

    V raises the degree of at most two links by two. Each odd link has
    increments 0,12; each even link 0,8. Polynomial spectral projectors on
    this finite residual support give the FULL second-order coefficient.
    """
    d = sector_data(k)
    steps = [12 if i < k else 8 for i in range(3)]
    shifts = sorted({0, *steps, *(steps[i]+steps[j] for i, j in combinations(range(3), 2))})
    residual, norm = d['residuals'][0], d['gram'][0][0]
    parts = {}
    for shift in shifts:
        p = residual
        for other in shifts:
            if other != shift:
                p = pscale(padd(h0(p), p, -(d['e']+other)), F(1, shift-other))
        parts[shift] = p
    weights = {s: inner(p, p)/norm for s, p in parts.items()}
    defects = {s: inner(padd(h0(p), p, -(d['e']+s)),
                       padd(h0(p), p, -(d['e']+s))) for s, p in parts.items()}
    total = {}
    for p in parts.values():
        total = padd(total, p)
    difference = padd(total, residual, -1)
    return {'weights': weights, 'eigen_defects': defects,
            'reconstruction_defect': inner(difference, difference),
            'coefficient': sum((w/s for s, w in weights.items() if s), F(0))}


def cubic_remainder(k, theta):
    """Full-space bound for the scalar branch, sector-bottom on [0,2]."""
    d, theta = sector_data(k), F(theta)
    if not 0 <= theta <= 2:
        raise ValueError('sector-bottom branch certified only on 0 <= theta <= 2')
    return 6*d['residual_squares'][0]*theta**3/(d['delta']*(d['delta']-d['a']*theta))


def second_order_window(k, theta):
    theta = F(theta)
    d = sector_data(k)
    remainder = cubic_remainder(k, theta)
    value = d['e']+d['a']*theta-source_components(k)['coefficient']*theta**2
    lower, upper = lower_upper(k, theta)
    return max(lower, value-remainder), min(upper, value+remainder)


def reflected_signature(p):
    signatures = {tuple(e[4*i] % 2 for i in range(3)) for e in p}
    return next(iter(signatures)) if len(signatures) == 1 else None


def is_diagonal(matrix):
    return all(x == 0 for i, row in enumerate(matrix) for j, x in enumerate(row) if i != j)


def singlet_counts(k):
    """Multiplicity of spin zero in (spin-0 + spin-1)^tensor k."""
    counts = {0: 1}
    for _ in range(k):
        new = dict(counts)  # scalar coordinate
        for spin, count in counts.items():
            for j in range(abs(spin-1), spin+2):
                new[j] = new.get(j, 0)+count
        counts = new
    return counts[0]


def run():
    checks, records = {}, []
    data = [sector_data(k) for k in range(4)]
    expected_slopes = [[F(9, 4)], [F(7, 4)], [F(4, 3), F(20, 9)],
                       [F(1), F(5, 3), F(5, 3), F(5, 3), F(10, 3)]]
    expected_variances = [[F(19, 16)], [F(43, 48)], [F(47, 72), F(715, 648)],
                          [F(11, 24), F(175, 216), F(175, 216), F(175, 216), F(25, 24)]]
    expected_c = [F(31, 256), F(311, 3840), F(451, 8640), F(19, 576)]

    checks['complete lowest gauge multiplets have dimensions 1,1,2,5'] = (
        [singlet_counts(k) for k in range(4)] == [len(d['basis']) for d in data] == [1, 1, 2, 5])
    checks['every retained vector is a simultaneous-conjugation singlet'] = all(
        not gauge_generator(p, a, b) for d in data for p in d['basis']
        for a, b in ((1, 2), (1, 3), (2, 3)))
    checks['retained vectors have precisely k odd centre links'] = all(
        all(tuple(sum(e[4*i:4*i+4]) % 2 for i in range(3)) == tuple(int(i < k) for i in range(3))
            for e in p) for k, d in enumerate(data) for p in d['basis'])
    checks['free spherical generator gives energy 3k on the entire retained multiplet'] = all(
        not padd(h0(p), p, -d['e']) for d in data for p in d['basis'])
    checks['next parity-allowed free shell is at least 8,8,8,12 above 3k'] = all(
        min(sum(n*(n+2) for n in ns)-3*k
            for ns in product(*(range(1, 6, 2) if i < k else range(0, 5, 2) for i in range(3)))
            if sum(n*(n+2) for n in ns) > 3*k) == d['delta']
        for k, d in enumerate(data))
    checks['Haar Gram and potential compressions are diagonal in the declared bases'] = all(
        is_diagonal(d['gram']) and is_diagonal(d['matrix']) for d in data)
    checks['all compressed first-order branches have exact rational slopes'] = (
        [d['slopes'] for d in data] == expected_slopes)
    checks['full residual Gram is diagonal with the exact rational squared norms'] = (
        all(is_diagonal(d['residual_gram']) for d in data)
        and [d['residual_squares'] for d in data] == expected_variances)
    checks['residuals are orthogonal to every retained vector, not only their own trial'] = all(
        inner(r, p) == 0 for d in data for r in d['residuals'] for p in d['basis'])
    checks['scalar-coordinate reflections distinguish the retained branches'] = all(
        len({reflected_signature(p) for p in d['basis']}) == len(d['basis']) for d in data)
    checks['potential respects scalar-coordinate reflections and centre parities'] = all(
        all(e[4*i] == 0 and sum(e[4*i:4*i+4]) % 2 == 0 for i in range(3))
        for e in potential())

    # Each expression (a_j-a_0) - theta b_j^2/(delta-theta a_j)
    # decreases on [0,2]. Positivity at 2 establishes the entire interval.
    margins = [[a-d['a']-2*b/(d['delta']-2*a)
                for a, b in zip(d['slopes'][1:], d['residual_squares'][1:])] for d in data]
    checks['scalar-product branch is the sector bottom for every 0<theta<=2'] = (
        all(d['delta']-2*a > 0 for d in data for a in d['slopes'])
        and all(m > 0 for row in margins for m in row))

    c0 = source_components(0)['coefficient']
    for k, d in enumerate(data):
        parts = source_components(k)
        checks[f'k={k}: entire residual resolves into exact free shells with no zero-energy part'] = (
            all(v == 0 for v in parts['eigen_defects'].values())
            and parts['reconstruction_defect'] == 0 and parts['weights'][0] == 0
            and sum(parts['weights'].values()) == d['residual_squares'][0]
            and parts['coefficient'] == expected_c[k])
        records.append({
            'odd_links': k, 'number_of_centre_sectors': (1, 3, 3, 1)[k],
            'free_bottom': d['e'], 'free_bottom_multiplicity_in_one_sector': len(d['basis']),
            'full_complement_free_increment': d['delta'],
            'compression_slopes': list(map(str, d['slopes'])),
            'residual_norm_squares': list(map(str, d['residual_squares'])),
            'a': str(d['a']), 'b2': str(d['b2']),
            'second_order_coefficient_c': str(parts['coefficient']),
            'splitting_second_order_coefficient': str(c0-parts['coefficient']),
            'exact_residual_spectral_weights': {str(s): str(w) for s, w in parts['weights'].items() if w},
            'branch_selection_margins_at_theta_2': list(map(str, margins[k])),
            'energy_enclosure_at_theta_2': list(map(str, lower_upper(k, 2))),
        })

    # Entire ordering interval: subtract the linear E_(k=1) upper.
    # Differences for k=2,3 decrease on [0,2], and the vacuum excited
    # level is >=8 from H>=H0. The first k=1 level is below its second (>=11).
    checks['single-odd sectors precede k=2,k=3 and vacuum excitations throughout [0,2]'] = (
        8 > lower_upper(1, 2)[1]
        and all(lower_upper(k, 2)[0] > lower_upper(1, 2)[1] for k in (2, 3))
        and data[1]['e']+data[1]['delta'] > lower_upper(1, 2)[1])
    checks['all-sector gap has the uniform lower 65/54 on [0,2]'] = (
        splitting(1, 2)[0] == F(65, 54) > 0)
    checks['gap curvature coefficient is exactly 77/1920'] = (
        c0-source_components(1)['coefficient'] == F(77, 1920))
    checks['YC2 Haar variance becomes the vacuum residual norm with factor four'] = (
        data[0]['b2'] == 4*F(19, 64))
    checks = {k: bool(v) for k, v in checks.items()}
    if not all(checks.values()):
        raise AssertionError([k for k, v in checks.items() if not v])
    return {'stage': 'YC3', 'checks': checks, 'all_pass': True,
            'reviewed_merge_head': '4aad5a0ec0fbdac727894ad0c315ed822cc0350e',
            'sectors': records,
            'enclosure': 'max(3k,3k+a theta-b2 theta^2/(delta-a theta)) <= E_k <= 3k+a theta',
            'second_order': 'E_k=3k+a theta-c theta^2+R; |R|<=6 b_scalar^2 theta^3/[delta(delta-a theta)] on [0,2]',
            'all_sector_gap': {
                'window': '0 <= theta <= 2', 'first_excitation': 'three single-centre-odd sector bottoms',
                'exact_multiplicity': 3, 'uniform_lower': '65/54',
                'gap_at_theta_1': list(map(str, splitting(1, 1))),
                'gap_at_theta_2': list(map(str, splitting(1, 2))),
                'expansion': '3-theta/2+(77/1920)theta^2+R_gap',
                'remainder': '|R_gap|<=theta^3[(43/8)/(8(8-7theta/4))+(57/8)/(8(8-9theta/4))]'},
            'claim_boundary': {
                'actual_quantum_operator': 'H_theta=-sum Delta_S3+2 theta sum cross-squares; gauge invariant one-site SU2^3',
                'static_Gibbs_compass_substituted_for_quantum_resolvent': False,
                'complete_hidden_space_floor_used': True,
                'weak_limit_centre_splitting_rate': False,
                'larger_lattice_or_continuum_mass_gap': False,
                'DR2_column_count_identified_with_volume': False,
                'twisted_trace_is_unique_possible_sector_observer': False,
                'status': 'written analytic proof with rational replay, not formal proof-assistant verification'},
            'source_sha256': source_pins()}


def source_pins():
    paths = ['physics/ol1/OL1_ONE_SITE_LATTICE.md',
             'physics/cb1/CB1_COMPACT_BRIDGE.md',
             'physics/cm2/CM2_COMPACT_CORE_MATCHING.md',
             'physics/cm2/CM2_RESULT.json',
             'physics/dr2/DR2_ALL_TURN_GAP.md', 'physics/dr2/DR2_RESULT.json',
             'physics/yc2/YC2_COMPACT_CENTRE_RESPONSE.md',
             'physics/yc3/YC3_SECTOR_SPLITTING.md',
             'physics/yc3/yc3_sector_splitting.py', 'physics/yc3/test_yc3.py']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC3_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale result or source pin')
    print('PASS', len(result['checks']), 'exact checks')
    print('All-sector gap on [0,2]: >=65/54; first excitation has exactly three centre characters.')
    print('Gap = 3-theta/2+(77/1920)theta^2 with the full cubic remainder bound.')
