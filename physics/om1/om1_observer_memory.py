"""OM1: an actual hidden floor and source-return enclosure for the core.

Frozen rational trial vectors are replayed with Fraction arithmetic only.
Higher moments retain the full Gaussian metric; no hidden matrix truncation
is identified with the full hidden inverse.
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


cm = load(PHYSICS/'cm1/cm1_core_sector_memory.py', 'om1_cm1')
cr = cm.cr
SECOND_LOWER = {3: F(5521, 1000), 4: F(4003, 500)}


def down(x, places=9):
    scale = 10**places
    return F((x*scale).numerator//(x*scale).denominator, scale)


def up(x, places=9):
    return -down(-x, places)


def quad(a, matrix, b=None):
    b = a if b is None else b
    return sum(x*sum(v*y for v, y in zip(row, b)) for x, row in zip(a, matrix))


def cut_bounds(eta, variance, second_lower):
    """Temple floor, true Q_psi h Q_psi floor and vacuum-overlap lower."""
    if variance < 0 or eta >= second_lower:
        raise ValueError('nonnegative variance and eta < the actual second-level lower are required')
    separation = second_lower-eta
    temple = eta-variance/separation
    delta = second_lower-variance/separation
    fidelity = separation**2/(separation**2+variance)
    return temple, delta, fidelity


@lru_cache(None)
def trial_packet(d):
    record = json.loads((HERE/'OM1_TRIALS.json').read_text())['trials'][str(d)]
    degree, width = record['degree'], F(record['width'])
    basis = [tuple(e) for e in record['basis']]
    if basis != cm.inv_basis(degree) or len(basis) != len(record['coefficients']):
        raise ValueError('frozen trial basis mismatch')
    p = {e: F(v) for e, v in zip(basis, record['coefficients']) if F(v)}
    mean = cr.free_means(d, width)
    inner = lambda a, b: mean(cr.pmul(a, b))
    vectors = [p]
    for _ in range(3):
        vectors.append(cm.scalar_h(vectors[-1], d, width))
    norm = inner(p, p)
    if norm <= 0:
        raise ValueError('trial norm must be positive')
    raw = []
    for n in range(7):
        raw.append(inner(vectors[n//2], vectors[(n+1)//2])/norm)
    symmetry = raw[2] == inner(vectors[0], vectors[2])/norm
    symmetry &= raw[3] == inner(vectors[0], vectors[3])/norm
    symmetry &= raw[4] == inner(vectors[1], vectors[3])/norm
    eta, variance = raw[1], raw[2]-raw[1]**2
    z = SECOND_LOWER[d]
    if not eta < z or not variance > 0:
        raise ValueError('ground threshold and nonzero source are required')
    temple, delta_exact, fidelity = cut_bounds(eta, variance, z)
    delta = down(delta_exact, 8)
    if delta <= eta:
        raise ValueError('computed hidden floor does not exceed the trial rate')

    # D^j r is a polynomial in h applied to psi; Q subtracts its full psi component.
    source_polynomials = [[-eta, F(1)]]
    for _ in range(2):
        shifted = [F(0)]+source_polynomials[-1]
        shifted[0] = -sum(c*raw[i] for i, c in enumerate(shifted))
        source_polynomials.append(shifted)
    if any(sum(c*raw[i] for i, c in enumerate(a)) for a in source_polynomials):
        raise ValueError('a hidden source moment lost orthogonality')
    moment_inner = lambda a, b: sum(x*y*raw[i+j] for i, x in enumerate(a) for j, y in enumerate(b))
    mu = [moment_inner(source_polynomials[k//2], source_polynomials[(k+1)//2]) for k in range(5)]
    if mu[0] != variance or mu[2] != moment_inner(source_polynomials[0], source_polynomials[2]):
        raise ValueError('source moment reconstruction failed')
    return {'degree': degree, 'width': width, 'size': len(basis), 'eta': eta,
            'variance': variance, 'temple': temple, 'delta_exact': delta_exact, 'delta': delta,
            'overlap_lower': fidelity, 'raw': raw, 'mu': mu,
            'source_polynomials': source_polynomials, 'vectors': vectors, 'norm': norm,
            'symmetry': symmetry, 'inner': inner}


def blocks(mu, size):
    if size not in (1, 2) or len(mu) < 2*size+1:
        raise ValueError('this packet supports one or two complete source directions')
    return {'s': [[mu[i+j] for j in range(size)] for i in range(size)],
            'h': [[mu[i+j+1] for j in range(size)] for i in range(size)],
            'h2': [[mu[i+j+2] for j in range(size)] for i in range(size)],
            'b': mu[:size], 'c': mu[1:size+1], 'norm': mu[0]}


def return_from_blocks(block, energy, delta):
    if energy >= delta:
        raise ValueError('energy must lie below a proved full hidden floor')
    s, h, h2 = (block[k] for k in ('s', 'h', 'h2'))
    a = [[x-energy*y for x, y in zip(row, sr)] for row, sr in zip(h, s)]
    coefficients = [row[0] for row in cm.solve_spd(a, [[v] for v in block['b']])]
    lower = sum(x*y for x, y in zip(coefficients, block['b']))
    cross = [x-energy*y for x, y in zip(block['c'], block['b'])]
    square = [[t-2*energy*x+energy**2*y for t, x, y in zip(tr, hr, sr)]
              for tr, hr, sr in zip(h2, h, s)]
    residual = block['norm']-2*sum(x*y for x, y in zip(coefficients, cross))+quad(coefficients, square)
    if residual < 0 or lower < 0:
        raise ValueError('a source/residual Gram is not positive')
    trial_upper = lower+residual/(delta-energy)
    bare_upper = block['norm']/(delta-energy)
    upper = min(trial_upper, bare_upper)
    if upper < lower:
        raise ValueError('inconsistent source-return enclosure')
    return {'lower': lower, 'upper': upper, 'trial_upper': trial_upper,
            'bare_upper': bare_upper, 'residual': residual, 'coefficients': coefficients}


def return_bounds(packet, energy, size):
    return return_from_blocks(blocks(packet['mu'], size), energy, packet['delta'])


def ground_window(packet, size, places=9):
    """Each final endpoint is independently signed; no monotone error-model assumption."""
    eta = packet['eta']
    scale = 10**places
    left = int(down(packet['temple']-F(1, 1000000), places)*scale)
    right = int(up(eta, places)*scale)
    def signs(n):
        energy = F(n, scale)
        ret = return_bounds(packet, energy, size)
        return eta-energy-ret['upper'], eta-energy-ret['lower']
    if signs(left)[0] <= 0 or signs(right)[1] >= 0:
        raise ValueError('initial Schur bracket is not certified')
    low, high = left, right
    while high-low > 1:
        mid = (low+high)//2
        if signs(mid)[0] > 0:
            low = mid
        else:
            high = mid
    lower = F(low, scale)
    low, high = left, right
    while high-low > 1:
        mid = (low+high)//2
        if signs(mid)[1] < 0:
            high = mid
        else:
            low = mid
    upper = F(high, scale)
    if not lower < upper or signs(int(lower*scale))[0] <= 0 or signs(int(upper*scale))[1] >= 0:
        raise ValueError('final Schur endpoint signs failed')
    return lower, upper


def run():
    checks, numbers = {}, {}
    sg = json.loads((PHYSICS/'sg1/SG1_RESULT.json').read_text())
    for d in (3, 4):
        p = trial_packet(d)
        checks[f'second-rate floor matches the SG1 certificate d={d}'] = SECOND_LOWER[d] == F(sg['numbers'][str(d)]['core_second_lower'])
        checks[f'full-action Gaussian moment identities d={d}'] = p['symmetry']
        checks[f'actual hidden floor above retained rate d={d}'] = p['delta_exact'] >= p['delta'] > p['eta']
        checks[f'positive vacuum-overlap lower d={d}'] = 0 < p['overlap_lower'] < 1
        checks[f'Temple and rank-one floor use the same residual d={d}'] = p['delta_exact']-p['temple'] == SECOND_LOWER[d]-p['eta']
        energy = down(p['eta'], 6)
        rows = []
        for size in (1, 2):
            memory = return_bounds(p, energy, size)
            low, high = ground_window(p, size)
            checks[f'positive full residual and ordered return d={d} size={size}'] = memory['residual'] >= 0 and 0 <= memory['lower'] <= memory['upper']
            checks[f'returned ground upper improves original reading d={d} size={size}'] = high < p['eta']
            rows.append({'hidden_trial_size': size, 'ground_lower': str(low), 'ground_upper': str(high),
                         'memory_lower_at_test_energy': str(down(memory['lower'], 12)),
                         'memory_upper_at_test_energy': str(up(memory['upper'], 12)),
                         'residual_norm_squared_upper': str(up(memory['residual'], 12)),
                         'bare_memory_upper': str(up(memory['bare_upper'], 12))})
        checks[f'two hidden directions improve the ground lower d={d}'] = F(rows[-1]['ground_lower']) > p['temple']
        checks[f'nested hidden trials increase the memory lower d={d}'] = return_bounds(p, energy, 2)['lower'] >= return_bounds(p, energy, 1)['lower']
        lo, hi = F(rows[-1]['ground_lower']), F(rows[-1]['ground_upper'])
        turning = F(75787, 10000) if d == 3 else F(54869, 5000)
        numbers[str(d)] = {'retained_trial_degree': p['degree'], 'retained_trial_coefficients': p['size'],
            'retained_rank': 1, 'source_rank': 1, 'eta_lower': str(down(p['eta'], 12)),
            'eta_upper': str(up(p['eta'], 12)), 'variance_upper': str(up(p['variance'], 12)),
            'temple_ground_lower': str(down(p['temple'], 9)), 'hidden_floor': str(p['delta']),
            'vacuum_overlap_lower': str(down(p['overlap_lower'], 9)),
            'test_energy': str(energy), 'return_refinements': rows,
            'gap_lower': str(down(SECOND_LOWER[d]-hi, 9)),
            'gap_upper': str(up(turning-lo, 9)),
            'normalized_raw_moments': [str(v) for v in p['raw']],
            'normalized_source_moments': [str(v) for v in p['mu']]}
    return {'schema': 'om1_observer_memory_v1', 'arithmetic': 'Frozen Fractions; exact means, solves and residual Grams',
            'checks': checks, 'numbers': numbers,
            'claim_boundary': {'full_hidden_floor_for_computed_rank_one_cut': 'PROVED_ON_FREE_GAUGE_SINGLET_CORE',
                'returned_memory': 'ENCLOSED_WITH_FULL_RESIDUAL_NOT_EXACTLY_REPLACED',
                'computed_reading_is_exact_vacuum': False,
                'CM1_67_reading_source_rank_26': 'UNCHANGED_DIFFERENT_CUT',
                'nonconstant_mode_matching': 'OPEN', 'continuum_mass_gap': 'OPEN'}}


def source_pins():
    paths = [PHYSICS/'cm1/cm1_core_sector_memory.py', PHYSICS/'cr1/cr1_core_rates.py',
             PHYSICS/'tc1/tc1_core_of_turns.py', PHYSICS/'sg1/sg1_singlet_gap.py',
             PHYSICS/'sg1/SG1_RESULT.json', PHYSICS/'sg1/SG1_GAUGE_SINGLET_GAP.md',
             HERE/'OM1_TRIALS.json', HERE/'om1_observer_memory.py',
             HERE/'OM1_OBSERVER_RETURN.md', HERE/'test_om1.py']
    return {str(path.relative_to(PHYSICS.parent)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if not all(result['checks'].values()):
        raise SystemExit('FAIL: '+', '.join(k for k, v in result['checks'].items() if not v))
    result['source_sha256'] = source_pins()
    path = HERE/'OM1_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale OM1 result or source pin')
    print('PASS', len(result['checks']), 'exact checks')
    for d, row in result['numbers'].items():
        print(d, 'hidden floor', row['hidden_floor'], 'ground window',
              row['return_refinements'][-1]['ground_lower'], row['return_refinements'][-1]['ground_upper'])


if __name__ == '__main__':
    main()
