"""YC5: full harmonic source measure and a two-sided return comparison."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'physics'/'yc4'))
import yc4_harmonic_return as c
y = c.y


def generator(p):
    """Same spherical generator as YC3, accumulated without dictionary copies."""
    out = {}
    for e, value in p.items():
        for link in range(3):
            degree = sum(e[4*link:4*link+4])
            out[e] = out.get(e, F(0))+value*degree*(degree+2)
        for i in range(12):
            if e[i] >= 2:
                f = list(e)
                f[i] -= 2
                f = tuple(f)
                out[f] = out.get(f, F(0))-value*e[i]*(e[i]-1)
    return {e: value for e, value in out.items() if value}


def possible_energies(p, floor):
    degrees = {tuple(sum(e[4*i:4*i+4]) for i in range(3)) for e in p}
    return sorted({sum(n*(n+2) for n in ns) for ds in degrees
                   for ns in product(*(range(d % 2, d+1, 2) for d in ds))
                   if sum(n*(n+2) for n in ns) >= floor})


@lru_cache(None)
def lagrange_coefficients(nodes, target):
    coefficients = [F(1)]
    for other in nodes:
        if other == target:
            continue
        out = [F(0)]*(len(coefficients)+1)
        for i, value in enumerate(coefficients):
            out[i] -= other*value/(target-other)
            out[i+1] += value/(target-other)
        coefficients = out
    return coefficients


def resolve(source, floor):
    nodes = tuple(possible_energies(source, floor))
    if not nodes:
        if c.inner(source, source):
            raise AssertionError('nonzero unresolved source below claimed floor')
        return {}, {'reconstruction': F(0), 'eigen_defect': F(0)}
    powers = [source]
    for _ in range(len(nodes)-1):
        powers.append(generator(powers[-1]))
    parts, summed, eigen_defect = {}, {}, F(0)
    for energy in nodes:
        p = {}
        for value, power in zip(lagrange_coefficients(nodes, energy), powers):
            p = y.padd(p, power, value)
        if c.inner(p, p):
            parts[energy] = p
            summed = y.padd(summed, p)
            defect = y.padd(generator(p), p, -energy)
            eigen_defect += c.inner(defect, defect)
    diff = y.padd(summed, source, -1)
    errors = {'reconstruction': c.inner(diff, diff), 'eigen_defect': eigen_defect}
    if any(errors.values()):
        raise AssertionError(('incomplete harmonic source resolution', errors))
    return parts, errors


@lru_cache(None)
def resolved_blocks(k):
    output = {}
    for sig, base in c.matrices(k).items():
        n = len(base['norms'])
        sources, resolutions, errors = [], [], []
        for i, (p, norm, _, _, _) in enumerate(base['vectors']):
            r = y.pmul(y.potential(), p)
            for j, vector in enumerate(base['vectors']):
                r = y.padd(r, vector[0], -base['V'][j][i]/base['norms'][j])
            parts, defects = resolve(r, base['floor'])
            sources.append(r)
            resolutions.append(parts)
            errors.append(defects)
        energies = sorted(set().union(*(p.keys() for p in resolutions)))
        weights = {}
        for energy in energies:
            matrix = [[F(0)]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1):
                    matrix[i][j] = matrix[j][i] = c.inner(resolutions[i].get(energy, {}),
                                                       resolutions[j].get(energy, {}))
            weights[energy] = matrix
        if any(sum(weights[e][i][j] for e in energies) != base['B'][i][j]
               for i in range(n) for j in range(n)):
            raise AssertionError('resolved matrices do not reconstruct full source Gram')
        output[sig] = {'base': base, 'sources': sources, 'parts': resolutions,
                       'weights': weights, 'defects': errors}
    return output


@lru_cache(None)
def hidden_interaction(k):
    """Actual QVQ matrix on the free-spectral source span of the scalar block."""
    if k not in (0, 1):
        raise ValueError('hidden-interaction refinement is built for k=0,1 only')
    sig = tuple(int(i < k) for i in range(3))
    block = resolved_blocks(k)[sig]
    basis = []
    for energy in block['weights']:
        at_energy = []
        for parts in block['parts']:
            p = parts.get(energy, {})
            for old, norm in at_energy:
                p = y.padd(p, old, -c.inner(p, old)/norm)
            norm = c.inner(p, p)
            if norm:
                at_energy.append((p, norm))
                basis.append((energy, p, norm))
    coords = [[c.inner(p, source)/norm for source in block['sources']] for _, p, norm in basis]
    images = [y.pmul(y.potential(), p) for _, p, _ in basis]
    size = len(basis)
    potential = [[F(0)]*size for _ in range(size)]
    for i in range(size):
        for j in range(i+1):
            potential[i][j] = potential[j][i] = c.inner(basis[i][1], images[j])
    by_energy = {e: [i for i, v in enumerate(basis) if v[0] == e] for e in block['weights']}
    rank = len(block['sources'])
    tensor = {}
    for e, ids in by_energy.items():
        for f, jds in by_energy.items():
            if f < e:
                continue
            applied = [[sum(potential[a][b]*coords[b][j] for b in jds)
                        for j in range(rank)] for a in ids]
            raw = [[sum(coords[a][i]*applied[t][j] for t, a in enumerate(ids))
                    for j in range(rank)] for i in range(rank)]
            tensor[(e, f)] = [[raw[i][j]+(raw[j][i] if e != f else 0)
                               for j in range(rank)] for i in range(rank)]
    return {'basis': basis, 'coords': coords, 'potential': potential, 'tensor': tensor}


def interaction_return(k, z):
    """< (D0-z)^-1 r_i, V (D0-z)^-1 r_j > in the true Haar metric."""
    h = hidden_interaction(k)
    rank = len(h['coords'][0])
    out = [[F(0)]*rank for _ in range(rank)]
    for (e, f), matrix in h['tensor'].items():
        denominator = (e-z)*(f-z)
        for i in range(rank):
            for j in range(rank):
                out[i][j] += matrix[i][j]/denominator
    return out


def pencil(block, theta, z, side, dynamic=False):
    theta, z = F(theta), F(z)
    base = block['base']
    if theta < 0 or z >= base['floor']:
        raise ValueError('requires theta>=0 and z below the full hidden floor')
    if side not in ('lower', 'upper'):
        raise ValueError('side must be lower or upper')
    out = c.pencil(base, base['k'], theta, z, False)
    n = len(out)
    refined = dynamic and base['k'] in (0, 1) and y.reflected_signature(base['vectors'][0][0]) == tuple(int(i < base['k']) for i in range(3))
    tangent = 2*theta/(base['floor']-z)
    for energy, weight in block['weights'].items():
        scale = theta**2/(energy+(6*theta if side == 'upper' else 0)-z)
        if refined and side == 'upper':
            scale = theta**2*(1+2*tangent)/((1+tangent)**2*(energy-z))
        for i in range(n):
            for j in range(n):
                out[i][j] -= scale*weight[i][j]
    if refined:
        interaction = interaction_return(base['k'], z)
        coefficient = theta**3/((1+tangent)**2 if side == 'upper' else 1+6*theta/(base['floor']-z))
        for i in range(n):
            for j in range(n):
                out[i][j] += coefficient*interaction[i][j]
    return out


def counts(k, theta, z, side, dynamic=False):
    answer = [0, 0, 0]
    for block in resolved_blocks(k).values():
        for i, value in enumerate(c.inertia(pencil(block, theta, z, side, dynamic))):
            answer[i] += value
    return tuple(answer)


def mesh():
    return [F(n, 8) for n in range(89)]+[F(176+n, 16) for n in range(1, 17)]


def validate_endpoint(record):
    theta = F(record['theta'])
    if not 0 <= theta <= 12:
        raise ValueError('the certificate window is [0,12]')
    if len(record['ground']) != 4 or any(len(pair) != 2 for pair in record['ground']) or len(record['second_lower']) != 2:
        raise ValueError('four ground pairs and two second-level floors required')
    for k, pair in enumerate(record['ground']):
        lower, upper = map(F, pair)
        if not 3*k <= lower <= upper:
            raise AssertionError(('invalid ground endpoints', k, theta))
        if counts(k, theta, lower, 'lower', True)[0] != 0:
            raise AssertionError(('failed ground lower', k, theta, lower))
        n, zero, _ = counts(k, theta, upper, 'upper', True)
        if n+zero < 1:
            raise AssertionError(('failed ground upper', k, theta, upper))
    for k, bound in enumerate(record['second_lower']):
        if counts(k, theta, F(bound), 'lower', True)[0] > 1:
            raise AssertionError(('failed second lower', k, theta, bound))


def validate_mesh(endpoints):
    if [F(r['theta']) for r in endpoints] != mesh():
        raise AssertionError('complete mesh required: eighths to 11, sixteenths to 12')


def coupling_cells(endpoints):
    validate_mesh(endpoints)
    cells = []
    for left, right in zip(endpoints, endpoints[1:]):
        upper1 = F(right['ground'][1][1])
        cells.append({
            'left': left['theta'], 'right': right['theta'],
            'gap_lower': str(F(left['ground'][1][0])-F(right['ground'][0][1])),
            'other_centre_margin': str(min(F(left['ground'][k][0])-upper1 for k in (2, 3))),
            'second_level_margin': str(min(F(z)-upper1 for z in left['second_lower']))})
    return cells


def source_pins():
    paths = ['physics/yc3/yc3_sector_splitting.py',
             'physics/yc4/YC4_HARMONIC_RETURN.md', 'physics/yc4/yc4_harmonic_return.py',
             'physics/yc4/YC4_RESULT.json', 'physics/yc5/YC5_RESOLVED_RETURN.md',
             'physics/yc5/yc5_resolved_return.py', 'physics/yc5/test_yc5.py',
             'physics/yc5/YC5_ENDPOINTS.json']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    endpoints = json.loads((HERE/'YC5_ENDPOINTS.json').read_text())['endpoints']
    validate_mesh(endpoints)
    predecessor = json.loads((ROOT/'physics/yc4/YC4_RESULT.json').read_text())
    checks = {'YC4 inputs match its frozen source pins': c.source_pins() == predecessor['source_sha256']}
    sectors, hidden = [], []
    for k in range(4):
        blocks = resolved_blocks(k)
        rank = sum(len(b['sources']) for b in blocks.values())
        checks[f'k={k}: complete retained cut unchanged'] = rank == (13, 21, 32, 23)[k]
        checks[f'k={k}: exact full source reconstruction and harmonic equations'] = all(
            not any(d.values()) for b in blocks.values() for d in b['defects'])
        checks[f'k={k}: resolved weights positive semidefinite'] = all(
            c.inertia(w)[0] == 0 for b in blocks.values() for w in b['weights'].values())
        checks[f'k={k}: resolved weights sum to full source Gram'] = all(
            sum(w[i][j] for w in b['weights'].values()) == b['base']['B'][i][j]
            for b in blocks.values() for i in range(len(b['sources'])) for j in range(len(b['sources'])))
        checks[f'k={k}: entire resolved source above full hidden floor'] = all(
            e >= 3*k+24 for b in blocks.values() for e in b['weights'])
        sectors.append({'odd_links': k, 'retained_rank': rank, 'full_hidden_floor': 3*k+24,
                        'source_energy_support': {''.join(map(str, sig)): list(b['weights'])
                                                  for sig, b in blocks.items()}})
    for k in (0, 1):
        h = hidden_interaction(k)
        basis, coords, potential = h['basis'], h['coords'], h['potential']
        sig = tuple(int(i < k) for i in range(3))
        block = resolved_blocks(k)[sig]
        n = len(basis)
        checks[f'k={k}: hidden representation has exact positive Haar metric'] = all(
            c.inner(p, q) == (norm if i == j else 0)
            for i, (_, p, norm) in enumerate(basis) for j, (_, q, _) in enumerate(basis)) and all(
                norm > 0 for _, _, norm in basis)
        checks[f'k={k}: hidden coordinates reconstruct every spectral weight'] = all(
            sum(norm*coords[a][i]*coords[a][j] for a, (energy, _, norm) in enumerate(basis) if energy == e) == w[i][j]
            for e, w in block['weights'].items() for i in range(len(w)) for j in range(len(w)))
        checks[f'k={k}: actual hidden potential form lies between 0 and 6'] = (
            c.inertia(potential)[0] == 0 and c.inertia([
                [(6*basis[i][2] if i == j else 0)-potential[i][j] for j in range(n)] for i in range(n)])[0] == 0)
        checks[f'k={k}: interaction tensor is symmetric'] = all(
            w[i][j] == w[j][i] for w in h['tensor'].values() for i in range(len(w)) for j in range(len(w)))
        hidden.append({'odd_links': k, 'reflection_signature': list(sig), 'source_span_rank': n,
                       'energy_ranks': {str(e): sum(v[0] == e for v in basis) for e in block['weights']}})
    for record in endpoints:
        validate_endpoint(record)
    checks['all 105 endpoints pass exact two-sided full-space counts'] = True
    cells = coupling_cells(endpoints)
    checks['all 104 coupling cells have gap at least 1/20'] = all(F(cell['gap_lower']) >= F(1, 20) for cell in cells)
    checks['all cells isolate precisely the three single-odd first excitations'] = all(
        F(cell['other_centre_margin']) > 0 and F(cell['second_level_margin']) > 0 for cell in cells)
    checks = {name: bool(ok) for name, ok in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    point_gaps = {r['theta']: [str(F(r['ground'][1][0])-F(r['ground'][0][1])),
                              str(F(r['ground'][1][1])-F(r['ground'][0][0]))]
                  for r in endpoints if F(r['theta']) in (2, 4, 6, 8, 10, 12)}
    return {'stage': 'YC5', 'all_pass': True, 'checks': checks, 'sectors': sectors,
            'hidden_interaction_representations': hidden,
            'coupling_window': '[0,12]', 'uniform_gap_lower': '1/20',
            'first_excitation_multiplicity': 3,
            'minimum_mesh_gap_lower': str(min(F(cell['gap_lower']) for cell in cells)),
            'minimum_other_centre_margin': str(min(F(cell['other_centre_margin']) for cell in cells)),
            'minimum_second_level_margin': str(min(F(cell['second_level_margin']) for cell in cells)),
            'point_gap_enclosures': point_gaps, 'coupling_cells': cells,
            'claim_boundary': {'one_site_only': True, 'all_gauge_and_centre_sectors_covered': True,
                               'uncomputed_complement_floor_proved': True,
                               'full_source_resolved': True, 'hidden_interaction_span_assumed_invariant': False,
                               'retained_cut_enlarged_relative_to_yc4': False,
                               'floating_point_used_to_certify': False,
                               'weak_limit_absolute_splitting_proved': False, 'volume_uniform_or_4d_gap': False,
                               'formal_proof_assistant_verification': False},
            'source_sha256': source_pins()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC5_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale result or source pins')
    print('PASS', len(result['checks']), 'exact checks; 105 endpoints, 104 complete coupling cells')
    print('All-sector gap >=1/20 on [0,12]; point gap windows:', result['point_gap_enclosures'])
