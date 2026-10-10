"""YC12: actual cube kinetics and a weakly joined correlated-cube lattice.

Exact rational certificate arithmetic. Analytic completeness, inverse order,
Bochner/min-max and the all-volume cluster proof are in the accompanying note.
No finite harmonic space is substituted for the full hidden compression.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cz = load('yc12_cz1', 'physics/cz1/cz1_centre_record_of_a_closed_surface.py')
y8 = load('yc12_yc8', 'physics/yc8/yc8_local_vacuum_dressing.py')
y10 = load('yc12_yc10', 'physics/yc10/yc10_closed_return_and_scale.py')
# Reuse the already tested symmetric rational congruence, including 2x2 pivots.
y4 = load('yc12_yc4', 'physics/yc4/yc4_harmonic_return.py')

GRAD_F2 = Q(77, 7488)
HESS_F2 = Q(193, 1248)


@lru_cache(None)
def cube():
    words = cz.box_surface(1, 1, 1)
    edges = tuple(sorted({e for word in words for e, sign in word}))
    number = {e: i for i, e in enumerate(edges)}
    faces = tuple(frozenset(number[e] for e, sign in word) for word in words)
    adjacent = tuple((p, q) for p, q in combinations(range(6), 2) if faces[p] & faces[q])
    opposite = tuple((p, q) for p, q in combinations(range(6), 2) if not faces[p] & faces[q])
    return edges, faces, adjacent, opposite


def cycle_supports(length):
    edges = cube()[0]
    found = []
    for support in combinations(range(12), length):
        degree = Counter(v for i in support for v in edges[i])
        if set(degree.values()) != {2}:
            continue
        visited = {next(iter(degree))}
        while True:
            new = visited | {v for i in support for v in edges[i] if any(w in visited for w in edges[i])}
            if new == visited:
                break
            visited = new
        if len(visited) == length:
            found.append(frozenset(support))
    return tuple(found)


def centre_flip():
    faces = cube()[1]
    return next(mask for mask in range(1 << 12)
                if all(sum((mask >> e) & 1 for e in face) % 2 for face in faces))


@lru_cache(None)
def source_grams():
    """Build Gram matrices from independent electric-harmonic features."""
    _, faces, adjacent, opposite = cube()
    features = {18: [(pair, Q(-1, 2), Q(1, 4)) for pair in adjacent],
                24: [(pair, Q(-2), Q(1, 16)) for pair in opposite],
                26: [(pair, Q(-2), Q(3, 64)) for pair in adjacent],
                32: [((p,), Q(-1, 2), Q(1)) for p in range(6)]}
    result = {}
    for energy, rows in features.items():
        matrix = [[Q(0) for _ in range(7)] for _ in range(7)]
        for indices, coefficient, norm in rows:
            for p in indices:
                for q in indices:
                    matrix[p+1][q+1] += coefficient**2 * norm
        result[energy] = matrix
    return result


def return_matrix(theta, zeta, upper):
    """Centered full-hidden parity envelope. The source is J-even."""
    theta, zeta = Q(theta), Q(zeta)
    if theta < 0 or zeta >= 18:
        raise ValueError('theta >=0 and zeta <18 required')
    ratio = 6*theta/(18-zeta)
    if ratio >= 1:
        raise ValueError('full hidden Neumann reserve must be positive')
    factor = 1/(1-ratio**2) if upper else Q(1)
    return [[factor*sum(matrix[i][j]/(energy-zeta)
                       for energy, matrix in source_grams().items())
             for j in range(7)] for i in range(7)]


def pencil(theta, zeta, lower):
    theta, zeta = Q(theta), Q(zeta)
    # More return means a smaller effective restoring operator.
    sigma = return_matrix(theta, zeta, upper=lower)
    result = [[-theta**2*sigma[i][j] for j in range(7)] for i in range(7)]
    for i in range(7):
        result[i][i] += (12 if i else 0)-zeta
    for i in range(1, 7):
        result[0][i] = result[i][0] = -theta/2
    return result


def spectral_enclosure(theta, index, bits=22):
    """Full gauge-carrier eigenvalue enclosure via exact Schur inertia.

Return centered energies; index 0 is the ground, index 1 the first
excitation. Rational pencil roots only bracket the actual eigenvalues.
"""
    theta = Q(theta)
    if not 0 <= theta <= 1 or index not in (0, 1):
        raise ValueError('0<=theta<=1 and index in {0,1} required')
    if theta == 0:
        return (Q(12*index), Q(12*index))
    left, right = Q(-1), Q(0) if index == 0 else Q(12)-Q(1, 10**8)
    endpoints = []
    for lower in (True, False):
        lo, hi = left, right
        if y4.inertia(pencil(theta, lo, lower))[0] > index:
            raise AssertionError('lower search endpoint does not enclose the pencil root')
        if y4.inertia(pencil(theta, hi, lower))[0] <= index:
            raise AssertionError('upper search endpoint does not enclose the pencil root')
        for _ in range(bits):
            mid = (lo+hi)/2
            if y4.inertia(pencil(theta, mid, lower))[0] <= index:
                lo = mid
            else:
                hi = mid
        endpoints.append(lo if lower else hi)
    return tuple(endpoints)


def cube_bounds(theta):
    theta = Q(theta)
    if not 0 <= theta <= 1:
        raise ValueError('certified cube window is 0<=theta<=1')
    rho = 2-Q(4, 3)*theta-2*HESS_F2*theta**2
    osc = 8*GRAD_F2*theta**3+12*GRAD_F2**2*theta**4
    return {'comparison_gap': rho, 'residual_oscillation': osc,
            'actual_full_gap_lower': rho-osc}


def tiling(shape):
    """Disjoint 12-edge cube factors and single bridge-edge factors."""
    if len(shape) != 3 or any(type(n) is not int or n < 4 or n % 2 for n in shape):
        raise ValueError('three even integer side lengths >=4 required')
    points = tuple(product(*(range(n) for n in shape)))

    def shift(x, axis):
        out = list(x)
        out[axis] = (out[axis]+1) % shape[axis]
        return tuple(out)

    def block(x):
        return tuple(a//2 for a in x)

    factors = {}
    for x in points:
        for axis in range(3):
            edge = (x, axis)
            other = shift(x, axis)
            factors[edge] = ('cube', block(x)) if block(x) == block(other) else ('bridge', edge)
    internal, external = [], []
    for x in points:
        for i, j in combinations(range(3), 2):
            face = ((x, i), (shift(x, i), j), (shift(x, j), i), (x, j))
            support = frozenset(factors[e] for e in face)
            (internal if len(support) == 1 else external).append(support)
    incidence = Counter(factor for support in external for factor in support)
    return {'edge_factors': factors, 'internal': internal, 'external': external,
            'incidence': incidence, 'shape': shape}


def join_bounds(eta):
    eta = Q(eta)
    if not 0 <= eta <= Q(1, 16384):
        raise ValueError('certified interface window is 0<=eta<=1/16384')
    g, k, s, r, M, exp_upper = Q(1, 4), 4, 1, Q(1, 8), 4, Q(3)
    beta, a = 24*eta, Q(2*k, s)
    mapping = M*beta/g*(1+2*r)*exp_upper
    contraction = M*beta/g*exp_upper*(2+a*(1+2*r))
    relative = 2*M*beta/g*exp_upper
    return dict(g=g, k=k, s=s, r=r, M=M, beta=beta, exponent=a*r,
                exp_upper=exp_upper, mapping=mapping, contraction=contraction,
                relative_return=relative, actual_full_gap_lower=g*(1-relative))


def encode(obj):
    if isinstance(obj, Q):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): encode(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode(x) for x in obj]
    return obj


def source_pins():
    paths = ['physics/cz1/cz1_centre_record_of_a_closed_surface.py',
             'physics/yc3/yc3_sector_splitting.py',
             'physics/yc4/yc4_harmonic_return.py',
             'physics/yc8/YC8_LOCAL_VACUUM_DRESSING.md',
             'physics/yc8/yc8_local_vacuum_dressing.py',
             'physics/yc9/YC9_SHARED_PLANE_GAP.md',
             'physics/yc9/yc9_shared_plane_gap.py',
             'physics/yc10/YC10_CLOSED_RETURN_AND_SCALE.md',
             'physics/yc10/yc10_closed_return_and_scale.py',
             'physics/yc11/YC11_CLOSED_CUBE_RESPONSE.md',
             'physics/yc12/YC12_CORRELATED_CUBE_GAP.md',
             'physics/yc12/yc12_correlated_cube_gap.py', 'physics/yc12/test_yc12.py']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    checks = {f'YC8 identity: {k}': bool(v) for k, v in y8.symbolic_checks().items()}
    edges, faces, adjacent, opposite = cube()
    checks['ordinary cube has twelve links and six faces'] = len(edges) == 12 and len(faces) == 6
    checks['six complete free energy-12 cycles'] = set(cycle_supports(4)) == set(faces)
    checks['full energy-18 shell includes four skew cycles'] = len(cycle_supports(6)) == 16
    checks['adjacent and opposite face incidence'] = len(adjacent) == 12 and len(opposite) == 3
    incidence = [(sum(e in f for f in faces),
                  sum(e in faces[p] | faces[q] for p, q in adjacent),
                  sum(e in faces[p] & faces[q] for p, q in adjacent)) for e in range(12)]
    checks['every edge sees two faces, seven pairs, one common pair'] = set(incidence) == {(2, 7, 1)}
    grad = Q(8, 4608)+Q(6, 1404)+Q(8, 1872)
    hess = Q(112, 4608)+Q(36, 1404)+Q(196, 1872)
    checks['cube derivative constants from incidence'] = grad == GRAD_F2 and hess == HESS_F2
    bounds = cube_bounds(1)
    checks['comparison curvature remains positive through one'] = bounds['comparison_gap'] == Q(223, 624)
    checks['finite residual oscillation is retained'] = bounds['residual_oscillation'] == Q(390313, 4672512)
    checks['actual cube gap includes all boundary charges'] = bounds['actual_full_gap_lower'] == Q(1279511, 4672512) > Q(1, 4)
    mask = centre_flip()
    checks['unitary link sign choice flips all six faces'] = all(sum(mask >> e & 1 for e in f) % 2 for f in faces)
    grams = source_grams()
    checks['source has no vacuum column'] = all(not any(m[0]) for m in grams.values())
    checks['full source Gram retains all cross columns'] = all(
        sum(m[i+1][j+1] for m in grams.values()) == (Q(3, 2) if i == j else Q(1, 4))
        for i in range(6) for j in range(6))
    checks['each electric source Gram is positive'] = all(y4.inertia(m)[0] == 0 for m in grams.values())
    # Octahedral face irreducibles: uniform, opposite-even traceless, opposite-odd.
    vectors = [[1]*6, [1, 1, -1, -1, 0, 0], [1, -1, 0, 0, 0, 0]]
    eigenvalues = {18: (Q(1, 2), Q(1, 8), Q(1, 4)),
                   24: (Q(1, 2), Q(1, 2), Q(0)),
                   26: (Q(3, 2), Q(3, 8), Q(3, 4)),
                   32: (Q(1, 4), Q(1, 4), Q(1, 4))}
    checks['resolved source respects all three face symmetry channels'] = all(
        sum(grams[e][i+1][j+1]*v[j] for j in range(6)) == eigenvalues[e][a]*v[i]
        for e in grams for a, v in enumerate(vectors) for i in range(6))
    samples = []
    for theta in (Q(1, 4), Q(1, 2), Q(1)):
        ground, first = spectral_enclosure(theta, 0), spectral_enclosure(theta, 1)
        gap = (first[0]-ground[1], first[1]-ground[0])
        checks[f'theta={theta}: full Schur spectral brackets ordered'] = ground[0] <= ground[1] < first[0] <= first[1]
        for index, enclosure in enumerate((ground, first)):
            nl = y4.inertia(pencil(theta, enclosure[0], True))
            nu = y4.inertia(pencil(theta, enclosure[1], False))
            checks[f'theta={theta}: level {index} exact outward inertias'] = nl[0] <= index < nu[0] and nl[1] == nu[1] == 0
        samples.append(dict(theta=theta, centered_ground=ground, centered_first=first,
                            actual_ground=[6*theta+x for x in ground],
                            actual_first=[6*theta+x for x in first], gauge_gap=gap))
    geometries = []
    for shape in ((4, 4, 4), (4, 6, 8), (6, 6, 6)):
        t = tiling(shape)
        sizes = Counter(t['edge_factors'].values())
        blocks = sum(f[0] == 'cube' for f in sizes)
        prefix = f'{shape}: '
        checks[prefix+'disjoint complete factors have twelve or one links'] = all(n == (12 if f[0] == 'cube' else 1) for f, n in sizes.items())
        checks[prefix+'each cube has its six internal plaquettes'] = len(t['internal']) == 6*blocks
        kinds = Counter(tuple(sorted(Counter(f[0] for f in x).items())) for x in t['external'])
        checks[prefix+'external plaquettes use exactly four factors'] = all(len(x) == 4 for x in t['external'])
        checks[prefix+'cube and bridge interface degrees are twenty-four and four'] = all(n == (24 if f[0] == 'cube' else 4) for f, n in t['incidence'].items())
        checks[prefix+'the two interface types have the claimed multiplicities'] = kinds == {(('bridge', 2), ('cube', 2)): 12*blocks, (('bridge', 4),): 6*blocks}
        geometries.append(dict(shape=shape, cube_factors=blocks,
                               bridge_factors=sum(f[0] == 'bridge' for f in sizes),
                               internal_faces=len(t['internal']), external_faces=len(t['external'])))
    join = join_bounds(Q(1, 16384))
    checks['exponential at one has a strict rational upper bound three'] = y10.exponential_upper(Q(1), 2) < 3
    checks['correlated-block map preserves the ball'] = join['mapping'] == Q(45, 512) < join['r']
    checks['complete cluster return contracts'] = join['contraction'] == Q(27, 32) < 1
    checks['relative hidden-return reserve stays positive'] = join['relative_return'] == Q(9, 64) < 1
    checks['actual joined-block gap is volume uniform'] = join['actual_full_gap_lower'] == Q(55, 256)
    if not all(checks.values()):
        raise AssertionError([k for k, v in checks.items() if not v])
    return encode(dict(stage='YC12', base_commit='8f15a7d', checks=checks,
                       cube_constants=dict(gradient_f2=grad, hessian_f2=hess, endpoint=bounds),
                       centre_flip_link_mask=mask, source_energies=list(grams), source_grams=grams,
                       finite_gauge_spectral_certificates=samples,
                       geometry_controls=geometries, uniform_join_certificate=join,
                       claim_boundary=dict(actual_full_cube_gap=True,
                           all_boundary_charge_sectors_kept=True, complete_hidden_inverse_enclosed=True,
                           actual_correlated_blocks_joined=True, volume_uniform_weak_interface_gap=True,
                           stronger_isotropic_coupling_window=False, physical_continuum_mass=False,
                           ST_identified_with_retained_Hamiltonian=False,
                           source_lambda_identified_with_theta=False, formal_machine_verification=False),
                       source_pins=source_pins()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC12_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC12: {len(result['checks'])} exact checks pass.")
    print('Actual correlated-cube block gap and volume-uniform weak-interface gap; no continuum claim.')
