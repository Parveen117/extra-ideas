"""Optional floating proposals; acceptance and replay use exact rational signs."""
from fractions import Fraction as F
from math import floor, ceil
import json
import yc6_vacuum_step as v
import explore_yc5 as e


def rounded(x, up=False):
    return F((ceil if up else floor)(x*10000), 10000)


if __name__ == '__main__':
    records = []
    for i, theta in enumerate(v.mesh()):
        t = float(theta)
        ground = [rounded(e.root(0, t, 'lower')), rounded(e.root(0, t, 'upper'), True)]
        lower1 = rounded(e.root(0, t, 'lower', 1))
        schur = e.root(0, t, 'upper', 1)
        ritz = e.values(0, t, 0, 'ritz')[1]
        method = 'schur' if schur is not None and schur <= ritz else 'ritz'
        upper1 = rounded(schur if method == 'schur' else ritz, True)
        if not theta:
            ground, lower1, upper1 = [F(0), F(0)], F(8), F(8)
        record = {'theta': str(theta), 'levels': [list(map(str, ground)), [str(lower1), str(upper1)]],
                  'excited_upper_method': method}
        v.validate_endpoint(record)
        records.append(record)
        if i % 16 == 0:
            print('exact vacuum endpoint signs through', theta, flush=True)
    rows = v.envelopes(records)
    cells = v.coupling_cells(records, rows)
    print('minimum full-cell reserve:', min(F(c['gap_lower']) for c in cells), flush=True)
    print('point gaps:', {r['theta']: r['gap'] for r in rows if F(r['theta']) in (2,4,6,8,10,12,14)}, flush=True)
    (v.HERE/'YC6_ENDPOINTS.json').write_text(json.dumps({
        'provenance': 'Floating proposals rounded to rationals; every spectral sign checked exactly. Replay imports no floating solver.',
        'mesh': 'eighths from 0 through 14', 'endpoints': records}, indent=2)+'\n')
