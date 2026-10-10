"""Floating proposals only; all certificate signs use rational arithmetic."""
from functools import lru_cache
from fractions import Fraction as F
from math import floor, ceil
import argparse
import json
import numpy as np
from scipy.linalg import eigvalsh
import yc5_resolved_return as m


@lru_cache(None)
def numeric(k):
    data = []
    for block in m.resolved_blocks(k).values():
        b = block['base']
        norms = np.array(b['norms'], float)
        scale = np.sqrt(norms[:, None]*norms[None, :])
        data.append((np.diag([v[3] for v in b['vectors']]), np.array(b['V'], float)/scale,
                     {e: np.array(w, float)/scale for e, w in block['weights'].items()}))
    return data


@lru_cache(None)
def numeric_hidden(k):
    h = m.hidden_interaction(k)
    sig = tuple(int(i < k) for i in range(3))
    norms = np.array(m.resolved_blocks(k)[sig]['base']['norms'], float)
    return (np.array([p[0] for p in h['basis']]),
            np.array(h['coords'], float)/np.sqrt(norms)[None, :], np.array(h['potential'], float))


def values(k, theta, z, side, dynamic=True):
    vals = []
    for sig, (kinetic, potential, weights) in zip(m.resolved_blocks(k), numeric(k)):
        matrix = kinetic+theta*potential-z*np.eye(len(kinetic))
        refined = dynamic and k < 2 and sig == tuple(int(i < k) for i in range(3))
        tangent = 2*theta/(3*k+24-z)
        if side != 'ritz':
            for energy, weight in weights.items():
                scale = theta**2/(energy+(6*theta if side == 'upper' else 0)-z)
                if refined and side == 'upper':
                    scale = theta**2*(1+2*tangent)/((1+tangent)**2*(energy-z))
                matrix -= scale*weight
            if refined:
                energies, coords, interaction = numeric_hidden(k)
                coeff = coords/(energies[:, None]-z)
                factor = theta**3/((1+tangent)**2 if side == 'upper' else 1+6*theta/(3*k+24-z))
                matrix += factor*(coeff.T@interaction@coeff)
        vals.extend(eigvalsh(matrix))
    return sorted(vals)


def root(k, theta, side, index=0):
    lo, hi = float(3*k), min(values(k, theta, 0, 'ritz')[index], 3*k+24-1e-7)
    if values(k, theta, hi, side)[index] > 1e-7:
        return None
    for _ in range(40):
        mid = (lo+hi)/2
        if values(k, theta, mid, side)[index] > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-proposals', action='store_true')
    args = parser.parse_args()
    for k in range(4):
        numeric(k)
        print('resolved', k, flush=True)
    if args.write_proposals:
        records = []
        for index, theta in enumerate(m.mesh()):
            ground = [[str(F(floor(root(k, float(theta), 'lower')*10000), 10000)),
                       str(F(ceil(root(k, float(theta), 'upper')*10000), 10000))] for k in range(4)]
            if theta == 0:
                ground = [[str(3*k)]*2 for k in range(4)]
            seconds = [str(F(floor(root(k, float(theta), 'lower', 1)*10000), 10000)) for k in (0, 1)]
            record = {'theta': str(theta), 'ground': ground, 'second_lower': seconds}
            m.validate_endpoint(record)
            records.append(record)
            if index % 8 == 0:
                print('exact signs validated through', theta, flush=True)
        cells = [F(a['ground'][1][0])-F(b['ground'][0][1]) for a, b in zip(records, records[1:])]
        print('minimum cell gap', min(cells), float(min(cells)), flush=True)
        (m.HERE/'YC5_ENDPOINTS.json').write_text(json.dumps({
            'provenance': 'Floating proposals rounded to rationals and checked by exact resolved-source, tangent/secant Schur counts. No floating arithmetic is used in replay.',
            'mesh': 'step 1/8 from 0 to 11; step 1/16 from 11 to 12', 'endpoints': records}, indent=2)+'\n')
        print('hidden ranks', [len(m.hidden_interaction(k)['basis']) for k in (0, 1)], flush=True)
    for theta in (2, 4, 6, 8, 10, 12, 14, 16):
        pairs = [(root(k, theta, 'lower'), root(k, theta, 'upper')) for k in range(4)]
        print(theta, pairs, 'seconds', [root(k, theta, 'lower', 1) for k in (0, 1)],
              'gap', None if any(v is None for v in pairs[0]+pairs[1]) else pairs[1][0]-pairs[0][1], flush=True)
