"""YC4 complete compact harmonic cuts. Exact rational certificate arithmetic."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import lcm, prod
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'physics'/'yc3'))
import yc3_sector_splitting as y


@lru_cache(None)
def packed_moment(encoded):
    return y.moment(tuple((encoded >> (5*i)) & 31 for i in range(12)))


@lru_cache(None)
def prepare(items):
    denom = lcm(*(c.denominator for _, c in items))
    groups, degrees = {}, [0, 0, 0]
    for e, c in items:
        if any(x >= 16 for x in e):
            raise ValueError('packed moment requires each input power <16')
        key = sum((x % 2) << i for i, x in enumerate(e))
        groups.setdefault(key, []).append((sum(x << (5*i) for i, x in enumerate(e)),
                                          c.numerator*(denom//c.denominator)))
        degrees = [max(degrees[i], sum(e[4*i:4*i+4])) for i in range(3)]
    return denom, groups, degrees


def inner(p, q):
    """Haar pairing, integer accumulation with parity-pruned monomial pairs."""
    dp, gp, np = prepare(tuple(sorted(p.items())))
    dq, gq, nq = prepare(tuple(sorted(q.items())))
    common = prod(prod(4+2*j for j in range((a+b)//2)) for a, b in zip(np, nq))
    total = 0
    for parity in gp.keys() & gq.keys():
        for e, c in gp[parity]:
            for f, d in gq[parity]:
                w = packed_moment(e+f)
                total += c*d*w.numerator*(common//w.denominator)
    return F(total, dp*dq*common)


def hi(p, link):
    out = {}
    for e, c in p.items():
        n = sum(e[4*link:4*link+4])
        out = y.padd(out, {e: c*n*(n+2)})
        for j in range(4*link, 4*link+4):
            if e[j] >= 2:
                f = list(e)
                f[j] -= 2
                out = y.padd(out, {tuple(f): -c*e[j]*(e[j]-1)})
    return out


def harmonic(p, ns):
    for i, n in enumerate(ns):
        for lower in range(n % 2, n, 2):
            lam = lower*(lower+2)
            p = y.pscale(y.padd(hi(p, i), p, -lam), F(1, n*(n+2)-lam))
    return p


def ppow(p, n):
    out = {y.ZERO: F(1)}
    for _ in range(n):
        out = y.pmul(out, p)
    return out


def spin_count(ns):
    return sum(abs(a-b) <= c <= a+b for a, b, c in product(*(range(n+1) for n in ns)))


@lru_cache(None)
def shell(ns):
    """Full invariant harmonic shell; exact Gram-Schmidt after SO3 contractions."""
    vectors = []
    for t in (0, 1):
        for b01, b02, b12 in product(range(max(ns)+1), repeat=3):
            ds = (b01+b02+t, b01+b12+t, b02+b12+t)
            if any(d > n for d, n in zip(ds, ns)):
                continue
            p = ppow(y.triple(), t)
            for pair, power in zip(((0, 1), (0, 2), (1, 2)), (b01, b02, b12)):
                p = y.pmul(p, ppow(y.dot(*pair), power))
            for i, (n, d) in enumerate(zip(ns, ds)):
                p = y.pmul(p, y.mono(*([4*i]*(n-d))))
            p = harmonic(p, ns)
            sig = y.reflected_signature(p)
            for old, norm, oldsig in vectors:
                if sig == oldsig:
                    p = y.padd(p, old, -inner(p, old)/norm)
            norm = inner(p, p)
            if norm:
                vectors.append((p, norm, sig))
    if len(vectors) != spin_count(ns):
        raise AssertionError(('incomplete shell', ns, len(vectors), spin_count(ns)))
    return vectors


def shells(k, increment=24):
    y.validate_k(k)
    cutoff = 3*k+increment
    if isinstance(increment, bool) or not isinstance(increment, int) or increment <= 0:
        raise ValueError('positive integer kinetic increment required')
    choices = [range(int(i < k), cutoff, 2) for i in range(3)]
    return sorted((ns for ns in product(*choices) if sum(n*(n+2) for n in ns) < cutoff),
                  key=lambda ns: (sum(n*(n+2) for n in ns), ns))


@lru_cache(None)
def matrices(k, increment=24):
    vectors = []
    for ns in shells(k, increment):
        energy = sum(n*(n+2) for n in ns)
        vectors.extend((p, norm, sig, energy, ns) for p, norm, sig in shell(ns))
    blocks = {}
    for sig in sorted({v[2] for v in vectors}):
        vs = [v for v in vectors if v[2] == sig]
        images = [y.pmul(y.potential(), v[0]) for v in vs]
        size = len(vs)
        norms = [v[1] for v in vs]
        vm = [[F(0)]*size for _ in vs]
        wm = [[F(0)]*size for _ in vs]
        for i in range(size):
            for j in range(i+1):
                vm[i][j] = vm[j][i] = inner(vs[i][0], images[j])
                wm[i][j] = wm[j][i] = inner(images[i], images[j])
        b = [[wm[i][j]-sum(vm[i][l]*vm[l][j]/norms[l] for l in range(size))
              for j in range(size)] for i in range(size)]
        blocks[sig] = {'k': k, 'floor': 3*k+increment,
                       'vectors': vs, 'norms': norms, 'V': vm, 'W': wm, 'B': b}
    return blocks


def inertia(matrix):
    """Exact symmetric congruence: (negative, zero, positive) including 2x2 pivots."""
    a = [list(map(F, row)) for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a) or any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('symmetric square matrix required')
    neg = zero = pos = 0
    while a:
        n = len(a)
        pivot = next((i for i in range(n) if a[i][i]), None)
        if pivot is not None:
            order = [pivot]+[i for i in range(n) if i != pivot]
            a = [[a[i][j] for j in order] for i in order]
            d = a[0][0]
            neg += int(d < 0)
            pos += int(d > 0)
            a = [[a[i][j]-a[i][0]*a[0][j]/d for j in range(1, n)] for i in range(1, n)]
        else:
            pair = next(((i, j) for i in range(n) for j in range(i) if a[i][j]), None)
            if pair is None:
                zero += n
                break
            order = list(pair)+[i for i in range(n) if i not in pair]
            a = [[a[i][j] for j in order] for i in order]
            b = a[0][1]
            neg += 1
            pos += 1
            a = [[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/b
                  for j in range(2, n)] for i in range(2, n)]
    return neg, zero, pos


def pencil(block, k, theta, z, returned):
    theta, z = F(theta), F(z)
    if block['k'] != k:
        raise ValueError('sector label must match the harmonic block')
    if theta < 0:
        raise ValueError('nonnegative coupling required')
    if returned and z >= block['floor']:
        raise ValueError('Schur threshold must lie strictly below full hidden floor')
    n = len(block['norms'])
    out = [[theta*block['V'][i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        out[i][i] += (block['vectors'][i][3]-z)*block['norms'][i]
    if returned:
        for i in range(n):
            for j in range(n):
                out[i][j] -= theta**2*block['B'][i][j]/(block['floor']-z)
    return out


def counts(k, theta, z, returned):
    totals = [0, 0, 0]
    for block in matrices(k).values():
        for i, v in enumerate(inertia(pencil(block, k, theta, z, returned))):
            totals[i] += v
    return tuple(totals)


def validate_endpoint(record):
    theta = F(record['theta'])
    if not 0 <= theta <= 8:
        raise ValueError('packet interval is [0,8]')
    if len(record['ground']) != 4 or any(len(pair) != 2 for pair in record['ground']) or len(record['second_lower']) != 2:
        raise ValueError('four ground pairs and both second-level floors are required')
    for k, pair in enumerate(record['ground']):
        lower, upper = map(F, pair)
        if not 3*k <= lower <= upper:
            raise AssertionError(('invalid ground endpoints', k, theta))
        if counts(k, theta, lower, True)[0] != 0:
            raise AssertionError(('failed full-space ground lower', k, theta, lower))
        n, zero, _ = counts(k, theta, upper, False)
        if n+zero < 1:
            raise AssertionError(('failed retained ground upper', k, theta, upper))
    for k, lower in enumerate(record['second_lower']):
        if counts(k, theta, F(lower), True)[0] > 1:
            raise AssertionError(('failed full-space second lower', k, theta, lower))


def source_pins():
    paths = ['physics/yc3/YC3_SECTOR_SPLITTING.md', 'physics/yc3/yc3_sector_splitting.py',
             'physics/yc3/YC3_RESULT.json', 'physics/cm2/CM2_RESULT.json',
             'physics/dr2/DR2_RESULT.json', 'physics/yc4/YC4_HARMONIC_RETURN.md',
             'physics/yc4/yc4_harmonic_return.py', 'physics/yc4/test_yc4.py',
             'physics/yc4/YC4_ENDPOINTS.json']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    endpoints = json.loads((HERE/'YC4_ENDPOINTS.json').read_text())['endpoints']
    if [F(r['theta']) for r in endpoints] != [F(i, 8) for i in range(65)]:
        raise AssertionError('the complete 1/8 mesh from 0 to 8 is required')
    checks, sectors = {}, []
    for k in range(4):
        data = matrices(k)
        vectors = [v for b in data.values() for v in b['vectors']]
        checks[f'k={k}: complete singlet shell rank'] = len(vectors) == sum(spin_count(ns) for ns in shells(k))
        checks[f'k={k}: exact Haar orthogonality and positive metric'] = all(
            inner(v[0], w[0]) == (v[1] if i == j else 0)
            for b in data.values() for i, v in enumerate(b['vectors']) for j, w in enumerate(b['vectors'])) and all(v[1] > 0 for v in vectors)
        checks[f'k={k}: free harmonic eigen-equations and centre characters'] = all(
            inner(y.padd(y.h0(p), p, -e), y.padd(y.h0(p), p, -e)) == 0
            and all(tuple(sum(power[4*i:4*i+4]) % 2 for i in range(3)) == tuple(int(i < k) for i in range(3))
                    for power in p) for p, norm, sig, e, ns in vectors)
        checks[f'k={k}: gauge invariance and reflection signatures'] = all(
            not y.gauge_generator(p, a, b) and y.reflected_signature(p) == sig
            for p, norm, sig, e, ns in vectors for a, b in ((1, 2), (1, 3), (2, 3)))
        checks[f'k={k}: full residual Gram is positive semidefinite'] = all(inertia(b['B'])[0] == 0 for b in data.values())
        sectors.append({'odd_links': k, 'retained_rank': len(vectors), 'full_hidden_floor': 3*k+24,
                        'shells': [{'degrees': list(ns), 'energy': sum(n*(n+2) for n in ns),
                                    'rank': spin_count(ns)} for ns in shells(k)],
                        'reflection_block_ranks': {''.join(map(str, s)): len(b['norms']) for s, b in data.items()},
                        'residual_rank': sum(inertia(b['B'])[2] for b in data.values())})
    for r in endpoints:
        validate_endpoint(r)
    checks['all 65 rational endpoints pass full-space lower counts and Ritz upper counts'] = True
    cells = []
    for left, right in zip(endpoints, endpoints[1:]):
        g = F(left['ground'][1][0])-F(right['ground'][0][1])
        upper1 = F(right['ground'][1][1])
        ordering = min(F(left['ground'][k][0])-upper1 for k in (2, 3))
        isolated = min(F(z)-upper1 for z in left['second_lower'])
        cells.append({'left': left['theta'], 'right': right['theta'],
                      'gap_lower': str(g), 'other_centre_margin': str(ordering),
                      'second_level_margin': str(isolated)})
    checks['all 64 complete coupling cells have gap at least 1/5'] = all(F(c['gap_lower']) >= F(1, 5) for c in cells)
    checks['all cells identify the three single-odd bottoms as the first excitation'] = all(
        F(c['other_centre_margin']) > 0 and F(c['second_level_margin']) > 0 for c in cells)
    checks = {name: bool(ok) for name, ok in checks.items()}
    if not all(checks.values()):
        raise AssertionError([name for name, ok in checks.items() if not ok])
    point_gaps = {r['theta']: [str(F(r['ground'][1][0])-F(r['ground'][0][1])),
                              str(F(r['ground'][1][1])-F(r['ground'][0][0]))]
                  for r in endpoints if F(r['theta']) in (2, 4, 6, 8)}
    return {'stage': 'YC4', 'all_pass': True, 'checks': checks, 'sectors': sectors,
            'coupling_window': '[0,8]', 'first_excitation_multiplicity': 3,
            'uniform_gap_lower': '1/5', 'minimum_mesh_gap_lower': str(min(F(c['gap_lower']) for c in cells)),
            'point_gap_enclosures': point_gaps, 'coupling_cells': cells,
            'claim_boundary': {'one_site_only': True, 'all_gauge_and_centre_sectors_covered': True,
                               'uncomputed_complement_floor_proved': True, 'floating_point_used_to_certify': False,
                               'weak_limit_absolute_splitting_proved': False, 'volume_uniform_or_4d_gap': False,
                               'formal_proof_assistant_verification': False},
            'source_sha256': source_pins()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC4_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale result or source pins')
    print('PASS', len(result['checks']), 'exact checks; 65 endpoint certificates, 64 complete coupling cells')
    print('All-sector gap >=1/5 on [0,8]; point gap windows:', result['point_gap_enclosures'])
