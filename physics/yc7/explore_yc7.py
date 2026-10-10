"""Propose endpoint rationals with float roots; exact replay decides every sign."""
from fractions import Fraction as F
import json
import math
import numpy as np
from scipy.optimize import brentq
import yc7_two_cell_return as y


def build():
    d = y.data()
    matrix = np.array(d['M'], float)
    weights = {lam: np.array(m, float) for lam, m in d['weights'].items()}

    def eig(theta, z, method, index):
        a = np.diag(np.array(y.ENERGIES, float)-z)+theta*matrix
        if method != 'ritz':
            for lam, m in weights.items():
                a -= theta**2*m/(lam-z+(12*theta if method == 'upper' else 0))
        return np.linalg.eigvalsh(a)[index]

    def root(theta, index, method):
        hi = 12-1e-7 if method != 'ritz' else 9+12*theta
        if eig(theta, hi, method, index) >= 0:
            return None
        return brentq(lambda z: eig(theta, z, method, index), 0, hi, xtol=1e-11)

    output = []
    for k in range(17):
        theta = F(k, 8)
        row = {'theta': str(theta)}
        for index, name in enumerate(('ground', 'second')):
            if not theta:
                row[name] = [str(y.ENERGIES[index])]*2
                row[name+'_upper_method'] = 'ritz'
                continue
            lower_root = root(float(theta), index, 'lower')
            lower = F(math.floor(lower_root*10000), 10000)
            while y.counts(theta, lower, 'lower')[0] > index:
                lower -= F(1, 10000)
            proposals = [(root(float(theta), index, m), m) for m in ('upper', 'ritz')]
            upper_root, method = min((r, m) for r, m in proposals if r is not None)
            upper = F(math.ceil(upper_root*10000), 10000)
            while sum(y.counts(theta, upper, method)[:2]) <= index:
                upper += F(1, 10000)
            row[name] = [str(lower), str(upper)]
            row[name+'_upper_method'] = method
        y.validate_endpoint(row)
        output.append(row)
    (y.HERE/'YC7_ENDPOINTS.json').write_text(json.dumps({'endpoints': output}, indent=2)+'\n')
    for row in output[::4]:
        l = F(row['second'][0])-F(row['ground'][1])
        u = F(row['second'][1])-F(row['ground'][0])
        print(row['theta'], 'gap', str(l), str(u), f'[{float(l):.4f},{float(u):.4f}]')
    print('minimum cell reserve', min(F(a['second'][0])-F(b['ground'][1]) for a,b in zip(output,output[1:])))


if __name__ == '__main__':
    build()
