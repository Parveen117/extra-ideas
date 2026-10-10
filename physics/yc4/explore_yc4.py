"""Non-certificate float proposals; replay never imports NumPy or this script."""
import numpy as np
from scipy.linalg import eigvalsh
from fractions import Fraction as F
from math import floor, ceil
import argparse
import json
import yc4_harmonic_return as m


def eigenvalues(k, theta, z=None):
    out = []
    for block in m.matrices(k).values():
        norms = np.array([float(x) for x in block['norms']])
        scale = np.sqrt(norms[:, None]*norms[None, :])
        mat = np.diag([v[3] for v in block['vectors']])+theta*np.array(block['V'], float)/scale
        if z is not None:
            mat -= np.eye(len(norms))*z+theta**2*np.array(block['B'], float)/scale/(3*k+24-z)
        out.extend(eigvalsh(mat))
    return sorted(out)


def lower(k, theta, index=0):
    a, b = float(3*k), min(eigenvalues(k, theta)[index], 3*k+24-1e-8)
    for _ in range(35):
        mid = (a+b)/2
        if eigenvalues(k, theta, mid)[index] > 0:
            a = mid
        else:
            b = mid
    return a


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-proposals', action='store_true')
    args = parser.parse_args()
    for k in range(4):
        m.matrices(k)
        print('built', k, flush=True)
    if args.write_proposals:
        records = []
        for n in range(65):
            theta = F(n, 8)
            ground = [[str(F(floor(lower(k, float(theta))*10000), 10000)),
                       str(F(ceil(eigenvalues(k, float(theta))[0]*10000), 10000))] for k in range(4)]
            if n == 0:
                ground = [[str(3*k)]*2 for k in range(4)]
            seconds = [str(F(floor(lower(k, float(theta), 1)*10000), 10000)) for k in (0, 1)]
            record = {'theta': str(theta), 'ground': ground, 'second_lower': seconds}
            m.validate_endpoint(record)
            records.append(record)
            print('validated', theta, flush=True)
        (m.HERE/'YC4_ENDPOINTS.json').write_text(json.dumps({
            'provenance': 'Float eigenvalue proposals rounded outward to rationals; every recorded bound passes exact full-space/retained inertia checks. Replay does not run this generator.',
            'mesh_step': '1/8', 'endpoints': records}, indent=2)+'\n')
    for t in (2, 3, 4, 6, 8, 10, 12):
        windows = [(lower(k, t), eigenvalues(k, t)[0]) for k in range(4)]
        print(t, windows, 'vac2 lower', lower(0, t, 1), 'gap', windows[1][0]-windows[0][1], flush=True)
