"""YC7: exact Haar/source calculation on the periodic 2 x 1 x 1 SU(2) lattice.

Six links: A0,A1,B0,C0,B1,C1. No gauge fixing or kinetic truncation.
The seven retained vectors exhaust vacuum electric energy below 12.
"""
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
sys.path.insert(0, str(ROOT/'physics'/'yc4'))
from yc4_harmonic_return import inertia

LINKS = 6
DIM = 4*LINKS
ZERO = (0,)*DIM
ONE = {ZERO: F(1)}
EDGES = ((0, 1), (1, 0), (0, 0), (0, 0), (1, 1), (1, 1))
ENERGIES = (0, 6, 6, 8, 8, 8, 8)
FLOOR = F(12)


def add(p, q, scale=F(1)):
    out = dict(p)
    for e, c in q.items():
        out[e] = out.get(e, F(0))+scale*c
        if not out[e]:
            del out[e]
    return out


def scale(p, c):
    return {e: c*v for e, v in p.items() if c*v}


def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            g = tuple(a+b for a, b in zip(e, f))
            out[g] = out.get(g, F(0))+c*d
    return {e: v for e, v in out.items() if v}


def mono(*indices):
    e = [0]*DIM
    for i in indices:
        e[i] += 1
    return {tuple(e): F(1)}


@lru_cache(None)
def moment(e):
    if len(e) != DIM or any(not isinstance(a, int) or a < 0 for a in e):
        raise ValueError('24 nonnegative integer powers required')
    if any(a % 2 for a in e):
        return F(0)
    ans = F(1)
    for i in range(LINKS):
        ks = [a//2 for a in e[4*i:4*i+4]]
        ans *= F(prod(prod(range(1, 2*k, 2)) for k in ks),
                 prod(4+2*j for j in range(sum(ks))))
    return ans


@lru_cache(None)
def packed_moment(e):
    return moment(tuple((e >> (5*i)) & 31 for i in range(DIM)))


@lru_cache(None)
def prepare(items):
    denom = lcm(*(c.denominator for _, c in items))
    groups, degrees = {}, [0]*LINKS
    for e, c in items:
        if any(a >= 16 for a in e):
            raise ValueError('packed pairing input power must be below 16')
        parity = sum((a % 2) << i for i, a in enumerate(e))
        packed = sum(a << (5*i) for i, a in enumerate(e))
        groups.setdefault(parity, []).append((packed, c.numerator*(denom//c.denominator)))
        degrees = [max(n, sum(e[4*i:4*i+4])) for i, n in enumerate(degrees)]
    return denom, groups, degrees


def inner(p, q):
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


def generator(p):
    out = {}
    for e, c in p.items():
        degree = sum(sum(e[4*i:4*i+4])*(sum(e[4*i:4*i+4])+2) for i in range(LINKS))
        out[e] = out.get(e, F(0))+c*degree
        for i in range(DIM):
            if e[i] >= 2:
                f = list(e)
                f[i] -= 2
                f = tuple(f)
                out[f] = out.get(f, F(0))-c*e[i]*(e[i]-1)
    return {e: c for e, c in out.items() if c}


def quat(link):
    return [mono(4*link+i) for i in range(4)]


def conjugate(q):
    return [q[0], *[scale(p, -1) for p in q[1:]]]


def qmul(p, q):
    out = [add(mul(p[0], q[0]), sum_polys(mul(p[i], q[i]) for i in range(1, 4)), -1)]
    for i, j, k in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        out.append(sum_polys([mul(p[0], q[i]), mul(p[i], q[0]),
                              mul(p[j], q[k]), scale(mul(p[k], q[j]), -1)]))
    return out


def sum_polys(ps):
    out = {}
    for p in ps:
        out = add(out, p)
    return out


def dot(i, j):
    return sum_polys(mono(4*i+a, 4*j+a) for a in (1, 2, 3))


@lru_cache(None)
def plaquettes():
    out = []
    for a, source, target in ((0, 2, 4), (1, 4, 2), (0, 3, 5), (1, 5, 3)):
        w = qmul(qmul(qmul(quat(a), quat(target)), conjugate(quat(a))),
                 conjugate(quat(source)))[0]
        out.append(add(ONE, w, -1))
    for i, j in ((2, 3), (4, 5)):
        out.append(scale(add(mul(dot(i, i), dot(j, j)), mul(dot(i, j), dot(i, j)), -1), 2))
    return tuple(out)


@lru_cache(None)
def potential_parts():
    ps = plaquettes()
    return sum_polys(ps[4:]), sum_polys(ps[:4])


@lru_cache(None)
def basis():
    return [ONE, scale(mono(8, 16), 4), scale(mono(12, 20), 4),
            *[add(scale(mono(4*i, 4*i), 4), ONE, -1) for i in (2, 3, 4, 5)]]


def gauge_generator(p, vertex, axis):
    """Infinitesimal q -> g_source q g_target^-1, quaternion convention."""
    unit = [{}, {}, {}, {}]
    unit[axis] = ONE
    coefficients = []
    for i, (source, target) in enumerate(EDGES):
        left, right = qmul(unit, quat(i)), qmul(quat(i), unit)
        coefficients.extend([add(scale(left[a], int(source == vertex)),
                                 right[a], -int(target == vertex)) for a in range(4)])
    out = {}
    for e, c in p.items():
        for i, n in enumerate(e):
            if n:
                f = list(e)
                f[i] -= 1
                out = add(out, mul({tuple(f): c*n}, coefficients[i]))
    return out


def centre_even(p):
    flips = ((0,), (2, 4), (3, 5))
    return all(sum(sum(e[4*i:4*i+4]) for i in flip) % 2 == 0
               for e in p for flip in flips)


def singlet_multiplicity(spins):
    """Tensor SU(2) singlet multiplicity; spins are twice physical spin."""
    counts = {0: 1}
    for spin in spins:
        nxt = {}
        for old, count in counts.items():
            for new in range(abs(old-spin), old+spin+1, 2):
                nxt[new] = nxt.get(new, 0)+count
        counts = nxt
    return counts.get(0, 0)


def electric_shells(cutoff=12):
    choices = [n for n in range(cutoff+1) if n*(n+2) < cutoff]
    out = []
    for ns in product(choices, repeat=LINKS):
        energy = sum(n*(n+2) for n in ns)
        if energy >= cutoff or ns[0] % 2 or (ns[2]+ns[4]) % 2 or (ns[3]+ns[5]) % 2:
            continue
        dims = [singlet_multiplicity([ns[0], ns[1], ns[2], ns[2], ns[3], ns[3]]),
                singlet_multiplicity([ns[0], ns[1], ns[4], ns[4], ns[5], ns[5]])]
        if prod(dims):
            out.append((ns, energy, prod(dims)))
    return sorted(out, key=lambda row: (row[1], row[0]))


def possible_energies(p):
    degrees = {tuple(sum(e[4*i:4*i+4]) for i in range(LINKS)) for e in p}
    return tuple(sorted({sum(n*(n+2) for n in ns) for ds in degrees
                        for ns in product(*(range(d % 2, d+1, 2) for d in ds))}))


def resolve(p):
    # Include every ambient harmonic possibility, even below the hidden floor.
    # Their zero norms are checked, rather than assumed during interpolation.
    nodes = possible_energies(p)
    powers = [p]
    for _ in range(len(nodes)-1):
        powers.append(generator(powers[-1]))
    parts = {}
    for lam in nodes:
        coefficients = [F(1)]
        for other in nodes:
            if other == lam:
                continue
            nxt = [F(0)]*(len(coefficients)+1)
            for i, value in enumerate(coefficients):
                nxt[i] -= other*value/(lam-other)
                nxt[i+1] += value/(lam-other)
            coefficients = nxt
        part = sum_polys(scale(power, c) for power, c in zip(powers, coefficients))
        if inner(part, part):
            defect = add(generator(part), part, -lam)
            if inner(defect, defect) or lam < FLOOR:
                raise AssertionError(('invalid hidden harmonic', lam))
            parts[lam] = part
    defect = add(sum_polys(parts.values()), p, -1)
    if inner(defect, defect):
        raise AssertionError('incomplete full-source reconstruction')
    return parts


@lru_cache(None)
def data():
    bs = basis()
    local, interface = potential_parts()
    v = add(local, interface)
    images = [mul(v, p) for p in bs]
    m = [[inner(p, q) for q in images] for p in bs]
    ml = [[inner(p, mul(local, q)) for q in bs] for p in bs]
    mi = [[m[i][j]-ml[i][j] for j in range(7)] for i in range(7)]
    sources = [add(q, sum_polys(scale(p, m[i][j]) for i, p in enumerate(bs)), -1)
               for j, q in enumerate(images)]
    b = [[inner(p, q) for q in sources] for p in sources]
    parts = [resolve(p) for p in sources]
    energies = sorted(set().union(*(d.keys() for d in parts)))
    weights = {lam: [[inner(p.get(lam, {}), q.get(lam, {})) for q in parts] for p in parts]
               for lam in energies}
    rl = [add(mul(local, q), sum_polys(scale(p, ml[i][j]) for i, p in enumerate(bs)), -1)
          for j, q in enumerate(bs)]
    ri = [add(sources[i], rl[i], -1) for i in range(7)]
    cross = [[inner(rl[i], ri[j])+inner(ri[i], rl[j]) for j in range(7)] for i in range(7)]
    return {'M': m, 'M_local': ml, 'M_interface': mi, 'B': b,
            'weights': weights, 'sources': sources, 'parts': parts,
            'source_local': rl, 'source_interface': ri, 'B_cross': cross}


def pencil(theta, z, which='lower'):
    theta, z = F(theta), F(z)
    if theta < 0 or which not in ('lower', 'upper', 'ritz'):
        raise ValueError('nonnegative theta and lower/upper/ritz required')
    if which != 'ritz' and z >= FLOOR:
        raise ValueError('Schur threshold must be below 12')
    d = data()
    out = [[theta*d['M'][i][j]+(ENERGIES[i]-z if i == j else 0)
            for j in range(7)] for i in range(7)]
    if which != 'ritz':
        for lam, matrix in d['weights'].items():
            denominator = F(lam)-z+(12*theta if which == 'upper' else 0)
            for i in range(7):
                for j in range(7):
                    out[i][j] -= theta**2*matrix[i][j]/denominator
    return out


def counts(theta, z, which='lower'):
    return inertia(pencil(theta, z, which))


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): encode(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode(v) for v in obj]
    return obj


def validate_endpoint(record):
    theta = F(record['theta'])
    if not 0 <= theta <= 2:
        raise ValueError('coupling must lie in [0,2]')
    for index, name in enumerate(('ground', 'second')):
        lower, upper = map(F, record[name])
        if not 0 <= lower <= upper or lower >= FLOOR:
            raise AssertionError('invalid energy window')
        if counts(theta, lower, 'lower')[0] > index:
            raise AssertionError(('failed lower', theta, index, lower))
        method = record[name+'_upper_method']
        if method not in ('upper', 'ritz'):
            raise ValueError('upper comparison must be Schur upper or Ritz')
        n, zero, _ = counts(theta, upper, method)
        if n+zero < index+1:
            raise AssertionError(('failed upper', theta, index, upper))


def source_pins():
    paths = ['physics/yc4/yc4_harmonic_return.py',
             'physics/yc6/YC6_VACUUM_STEP.md', 'physics/yc6/YC6_RESULT.json',
             'physics/cb1/CB1_COMPACT_BRIDGE.md',
             'physics/SIGNED_COMPASS_YM_HANDOFF.md',
             'physics/yc7/YC7_TWO_CELL_RETURN.md',
             'physics/yc7/yc7_two_cell_return.py', 'physics/yc7/test_yc7.py',
             'physics/yc7/YC7_ENDPOINTS.json']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    d = data()
    bs, ps = basis(), plaquettes()
    checks = {}
    shells = electric_shells()
    checks['complete vacuum cut has seven readings at energies 0,6,6,8,8,8,8'] = (
        sorted(e for _, e, multiplicity in shells for _ in range(multiplicity)) == list(ENERGIES))
    checks['the next allowed free shell is exactly 12'] = min(
        e for _, e, _ in electric_shells(13) if e >= 12) == 12
    checks['Haar Gram is identity'] = all(inner(p, q) == int(i == j)
        for i, p in enumerate(bs) for j, q in enumerate(bs))
    checks['retained free eigen-equations'] = all(not add(generator(p), p, -e)
        for p, e in zip(bs, ENERGIES))
    checks['local Gauss generators annihilate every plaquette and retained reading'] = all(
        not gauge_generator(p, v, a) for p in (*ps, *bs) for v in (0, 1) for a in (1, 2, 3))
    checks['three global centre flips preserve the retained readings and every plaquette'] = all(
        centre_even(p) for p in (*ps, *bs))
    checks['six plaquette Haar means are 1,1,1,1,3/4,3/4'] = (
        [inner(ONE, p) for p in ps] == [F(1)]*4+[F(3, 4)]*2)
    checks['all retained-source projections vanish'] = all(inner(p, q) == 0
        for p in bs for q in d['sources'])
    checks['full source Gram reconstructs from all electric energies'] = all(
        sum(m[i][j] for m in d['weights'].values()) == d['B'][i][j]
        for i in range(7) for j in range(7))
    checks['each resolved source measure is positive'] = all(
        inertia(m)[0] == 0 for m in d['weights'].values())
    checks['ground hidden source has weights 3/4 at 14 and 7/24 at 16'] = (
        {lam: m[0][0] for lam, m in d['weights'].items() if m[0][0]} == {14: F(3, 4), 16: F(7, 24)})
    checks['ground full variance and hidden variance retain their different cuts'] = (
        inner(add(sum_polys(ps), ONE, -F(11, 2)), add(sum_polys(ps), ONE, -F(11, 2))) == F(43, 24)
        and d['B'][0][0] == F(25, 24))
    coefficient = sum(d['M'][i][0]**2/F(e) for i, e in enumerate(ENERGIES) if e)
    coefficient += sum(m[0][0]/lam for lam, m in d['weights'].items())
    checks['full second-order ground coefficient is 167/896'] = coefficient == F(167, 896)
    checks['local/interface cross source is nonzero and retained'] = d['B_cross'][1][3] == F(1, 4)
    checks['cross source reconstructs the full Gram with both source families'] = all(
        inner(d['source_local'][i], d['source_local'][j])+
        inner(d['source_interface'][i], d['source_interface'][j])+d['B_cross'][i][j] == d['B'][i][j]
        for i in range(7) for j in range(7))
    checks['retained potential obeys 0 <= M <= 12'] = (
        inertia(d['M'])[0] == 0 and inertia([[12*int(i == j)-d['M'][i][j] for j in range(7)] for i in range(7)])[0] == 0)

    endpoints = json.loads((HERE/'YC7_ENDPOINTS.json').read_text())['endpoints']
    if [F(row['theta']) for row in endpoints] != [F(i, 8) for i in range(17)]:
        raise AssertionError('complete eighth mesh from 0 to 2 is required')
    for row in endpoints:
        validate_endpoint(row)
    checks['17 endpoint enclosures pass exact full-space inertia'] = True
    reserves = [F(a['second'][0])-F(b['ground'][1]) for a, b in zip(endpoints, endpoints[1:])]
    checks['all 16 coupling cells have vacuum gap at least 6/5'] = min(reserves) >= F(6, 5)
    if not all(checks.values()):
        raise AssertionError([k for k, ok in checks.items() if not ok])
    points = [{**row, 'gap': [F(row['second'][0])-F(row['ground'][1]),
                               F(row['second'][1])-F(row['ground'][0])]} for row in endpoints]
    return encode({'stage': 'YC7', 'base_commit': 'e7282d5', 'checks': checks,
        'carrier': 'periodic 2x1x1, six SU(2) links; two local Gauss laws; three even global centre characters',
        'rank': 7, 'hidden_floor': FLOOR, 'electric_shells': shells,
        'matrices': {k: d[k] for k in ('M', 'M_local', 'M_interface', 'B', 'B_cross', 'weights')},
        'ground_second_order': coefficient, 'uniform_gap': F(6, 5),
        'minimum_cell_reserve': min(reserves), 'point_windows': points,
        'claim_boundary': {'full_vacuum_sector': True, 'hidden_space_truncated': False,
            'volume_uniform': False, 'continuum_4d': False,
            'physical_mass_scale_calibrated': False, 'lowest_excited_symmetry_identified': False},
        'source_pins': source_pins()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'YC7_RESULT.json'
    if args.check:
        if result != json.loads(path.read_text()):
            raise AssertionError('result or source pins differ from the frozen packet')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(f"YC7: {len(result['checks'])} exact checks pass; gap >= 6/5 on [0,2].")
    print('Minimum coupling-cell reserve:', result['minimum_cell_reserve'])
