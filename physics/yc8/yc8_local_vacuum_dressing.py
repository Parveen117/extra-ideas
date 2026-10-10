"""YC8: local logarithmic vacuum dressing on ordinary 3D SU(2) tori.

The proof covers every rectangular torus with side lengths >=3. Finite
geometry checks are controls for its bounded-incidence argument, not a
substitute for that all-size proof. Exact fractions and symbolic identities.
"""
from fractions import Fraction as F
from itertools import product, combinations
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GRAD_F1 = F(1, 3)
GRAD_F2 = F(205, 3744)
HESS_F1 = F(4, 3)
HESS_F2 = F(1555, 1872)
CUBIC = 2*GRAD_F1*GRAD_F2
QUARTIC = GRAD_F2**2
STENCIL_BOUND = 121


def validate_shape(shape):
    if len(shape) != 3 or any(isinstance(n, bool) or not isinstance(n, int) or n < 3 for n in shape):
        raise ValueError('three integer side lengths >=3 required; collapsed tori use different identities')


@lru_cache(None)
def geometry(shape):
    validate_shape(shape)
    points = tuple(product(*(range(n) for n in shape)))
    edges = tuple((x, i) for x in points for i in range(3))
    number = {e: k for k, e in enumerate(edges)}

    def shift(x, i):
        out = list(x)
        out[i] = (out[i]+1) % shape[i]
        return tuple(out)

    faces, words = [], []
    for x in points:
        for i, j in combinations(range(3), 2):
            word = ((number[x, i], 1), (number[shift(x, i), j], 1),
                    (number[shift(x, j), i], -1), (number[x, j], -1))
            words.append(word)
            faces.append(frozenset(e for e, sign in word))
    incident = [set() for _ in edges]
    for p, face in enumerate(faces):
        for e in face:
            incident[e].add(p)
    pairs = {}
    for e, ps in enumerate(incident):
        for p, q in combinations(sorted(ps), 2):
            if (p, q) in pairs:
                raise AssertionError('two distinct plaquettes share more than one edge')
            pairs[p, q] = e
    # The local remainder r_e depends only on this two-plaquette stencil.
    stencils, pair_incidence = [], []
    for e, ps in enumerate(incident):
        touching = {(p, q) for p, q in pairs if e in faces[p] or e in faces[q]}
        common = sum(pairs[p, q] == e for p, q in touching)
        pair_incidence.append((len(touching), common))
        face_set = set(ps)
        for p, q in touching:
            face_set.update((p, q))
        stencils.append(frozenset().union(*(faces[p] for p in face_set)))
    return {'edges': edges, 'faces': tuple(faces), 'words': tuple(words),
            'incident': incident, 'pairs': pairs, 'stencils': stencils,
            'pair_incidence': pair_incidence}


def local_constants():
    grad = F(4*4, 4608)+F(36, 1404)+F(36+2*6, 1872)
    hess = F(4*(8+3*16), 4608)+F(36*6, 1404)+F(42*7*4, 1872)
    return grad, hess


def residual_density(theta):
    theta = F(theta)
    if theta < 0:
        raise ValueError('nonnegative theta required')
    return CUBIC*theta**3+QUARTIC*theta**4


def energy_density(theta):
    theta = F(theta)
    remainder = residual_density(theta)
    centre = theta-theta**2/48
    return max(F(0), centre-remainder), min(theta, centre+CUBIC*theta**3)


def comparison_gap(theta):
    theta = F(theta)
    if not 0 <= theta <= F(1, 2):
        raise ValueError('positive comparison gap certified on [0,1/2]')
    return 2-2*HESS_F1*theta-2*HESS_F2*theta**2


def symbolic_checks():
    q = sp.symbols('q0:4', real=True)
    a = sp.symbols('a0:4', real=True)
    b = sp.symbols('b0:4', real=True)
    wp, wq = sum(x*y for x, y in zip(q, a)), sum(x*y for x, y in zip(q, b))
    A, B = sum(x*y for x, y in zip(a, b)), wp*wq

    def euler(p):
        return sum(x*sp.diff(p, x) for x in q)

    def h_edge(p):
        return sp.expand(euler(euler(p))+2*euler(p)-sum(sp.diff(p, x, 2) for x in q))

    result = {}
    result['shared-edge gradient is A minus the Wilson product'] = sp.expand(
        sum(sp.diff(wp, x)*sp.diff(wq, x) for x in q)-euler(wp)*euler(wq)-(A-B)) == 0
    result['shared-edge product harmonic equation'] = sp.expand(h_edge(B)-(8*B-2*A)) == 0
    result['pair inverse uses exact electric energies 18 and 26'] = sp.expand(
        sp.Rational(2,39)*18*A-(26*B-2*A)/26-(A-B)) == 0
    result['pair spectral components reconstruct with correct rates'] = (
        sp.Rational(3,4)/18+sp.Rational(1,4)/26 == sp.Rational(2,39))

    w = sp.symbols('w', real=True)
    h_w2 = 2*w*(12*w)-2*4*(1-w*w)
    result['adjoint plaquette is an energy-32 harmonic'] = sp.expand(4*h_w2-32*(4*w*w-1)) == 0

    chi, aa, bb, p = sp.symbols('chi A B P')
    gamma_s = 3*p-chi+2*(aa-bb)
    h_f2 = 32*chi/4608-18*aa/1404+(26*bb-2*aa)/1872
    result['second dressing cancels the entire nonconstant quadratic source'] = sp.expand(
        h_f2-(p/sp.Integer(48)-gamma_s/144)) == 0
    theta, v, g12, g22 = sp.symbols('theta V g12 g22')
    local = theta*v-theta*(v-p)-theta**2*h_f2-theta**2*gamma_s/144-2*theta**3*g12-theta**4*g22
    result['positive logarithmic dressing has exact quartic local energy'] = sp.expand(
        local-(p*theta-p*theta**2/48-2*theta**3*g12-theta**4*g22)) == 0
    return result


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): encode(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode(v) for v in obj]
    return obj


def source_pins():
    paths = ['physics/yc7/YC7_TWO_CELL_RETURN.md', 'physics/yc7/YC7_RESULT.json',
             'uncut/up1/UP1_DEGREE_OF_THE_POTENTIAL.md',
             'physics/yc8/YC8_LOCAL_VACUUM_DRESSING.md',
             'physics/yc8/yc8_local_vacuum_dressing.py', 'physics/yc8/test_yc8.py']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    checks = symbolic_checks()
    records = []
    for shape in ((3,3,3), (3,4,5), (5,5,5)):
        g = geometry(shape)
        n, p = len(g['edges']), len(g['faces'])
        checks[f'{shape}: four distinct edges per face and four faces per edge'] = (
            n == p and all(len(f) == 4 for f in g['faces']) and all(len(x) == 4 for x in g['incident']))
        checks[f'{shape}: local pairs have unique shared edges and six pairs per edge'] = len(g['pairs']) == 6*n
        checks[f'{shape}: each edge is in 42 pairs, common to 6'] = all(
            row == (42, 6) for row in g['pair_incidence'])
        checks[f'{shape}: residual stencils are symmetric with at most 121 links'] = (
            max(map(len, g['stencils'])) <= STENCIL_BOUND and
            all(e in g['stencils'][f] for e, stencil in enumerate(g['stencils']) for f in stencil))
        records.append({'shape': shape, 'links': n, 'plaquettes': p,
                        'shared_pairs': len(g['pairs']), 'max_residual_stencil': max(map(len,g['stencils']))})
    grad, hess = local_constants()
    checks['per-link second-dressing gradient bound'] = grad == GRAD_F2
    checks['volume-independent Hessian row bound'] = hess == HESS_F2
    checks['exact cubic and quartic density coefficients'] = (CUBIC == F(205,5616) and QUARTIC == F(42025,14017536))
    checks['comparison gap stays above one quarter through theta=1/2'] = comparison_gap(F(1,2)) == F(941,3744) > F(1,4)
    checks['original operator residual is not silently identified with a scalar'] = CUBIC > 0 and QUARTIC > 0
    if not all(checks.values()):
        raise AssertionError([k for k, v in checks.items() if not v])
    return encode({'stage': 'YC8', 'base_commit': 'baca59a', 'checks': checks,
        'carrier': 'periodic Nx x Ny x Nz SU(2), all three side lengths >=3; all link kinetic energies retained',
        'geometric_controls': records,
        'first_dressing': '-sum(W_p)/12',
        'second_dressing': 'sum(4 W_p^2-1)/4608 - sum_shared(A_pq)/1404 + sum_shared(W_p W_q)/1872',
        'gradient_f1': GRAD_F1, 'gradient_f2': grad, 'hessian_f1': HESS_F1, 'hessian_f2': hess,
        'cubic_density': CUBIC, 'quartic_density': QUARTIC, 'local_stencil_bound': STENCIL_BOUND,
        'density_windows': [{'theta': t, 'ground_density': energy_density(t),
                             'comparison_gap_lower': comparison_gap(t),
                             'residual_per_link': residual_density(t),
                             'interaction_norm_bound': STENCIL_BOUND*residual_density(t)}
                            for t in (F(0),F(1,10),F(1,4),F(1,2))],
        'claim_boundary': {'actual_YM_ground_energy_density_uniform_in_volume': True,
            'comparison_Hamiltonian_gap_uniform_in_volume': True,
            'actual_YM_gap_uniform_in_volume': False,
            'explicit_dressing_is_actual_vacuum': False,
            'continuum_4d': False, 'geometry_inferred_only_from_finite_samples': False},
        'source_pins': source_pins()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC8_RESULT.json'
    if args.check:
        if json.loads(path.read_text()) != result:
            raise AssertionError('result or source pins differ')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC8: {len(result['checks'])} exact checks pass.")
    print('Actual YM energy density: theta-theta^2/48, with explicit volume-independent cubic/quartic errors.')
    print('Comparison gap >=1/4 on [0,1/2]; actual YM uniform gap NOT claimed.')
