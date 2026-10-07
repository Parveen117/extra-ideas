"""PH2: the return angle is the rotation left by a cycle of pure strains.

Exact part (stdlib): two pure strains compose to a strain times a rotation;
the polar angle formula.  Numerical part (numpy + CoolProp): the chain of N
symmetric positive elements that follows the CO2 cycle of PH1 and the rotation
it leaves.  Python 3.12.
"""
from fractions import Fraction as F
import json
import math
import sys


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def polar_tangent(g):
    """(numerator, denominator) of tan(theta) for g = R(theta) U, U symmetric positive."""
    return g[1][0]-g[0][1], g[0][0]+g[1][1]


def exact_witness():
    c, s = F(3, 5), F(4, 5)
    q = [[c, -s], [s, c]]
    qt = [[c, s], [-s, c]]
    m1 = [[F(2), F(0)], [F(0), F(1, 2)]]
    m2 = mul(mul(q, [[F(3), F(0)], [F(0), F(1, 3)]]), qt)
    for m in (m1, m2):
        if m[0][1] != m[1][0] or m[0][0] <= 0 or m[0][0]*m[1][1]-m[0][1]**2 <= 0:
            raise ValueError('element is not a pure strain')
    g = mul(m2, m1)
    num, den = polar_tangent(g)
    # R^T g symmetric with (cos, sin) proportional to (den, num)
    left = den*g[0][1]+num*g[1][1]
    right = -num*g[0][0]+den*g[1][0]
    if left != right:
        raise ValueError('polar angle formula failed')
    if num == 0:
        raise ValueError('two non-aligned pure strains must leave a rotation')
    aligned = mul([[F(3), F(0)], [F(0), F(1, 3)]], m1)
    if polar_tangent(aligned)[0] != 0:
        raise ValueError('aligned strains must leave no rotation')
    return dict(product=[[str(x) for x in row] for row in g], tan_theta=str(F(num, den)),
                aligned_tan_theta='0')


def chain(Hs):
    """Elements M_k and the final F for a closed list of responses (numpy)."""
    import numpy as np
    f = np.linalg.cholesky(Hs[0]).T          # upper triangular, F^T F = H_0
    f0 = f.copy()
    elements = []
    for h in Hs[1:]+[Hs[0]]:
        fi = np.linalg.inv(f)
        a = fi.T@h@fi
        ev, o = np.linalg.eigh(a)
        m = o@np.diag(np.sqrt(ev))@o.T
        elements.append(m)
        f = m@f
    rot = f@np.linalg.inv(f0)
    theta = math.atan2(rot[1, 0]-rot[0, 1], rot[0, 0]+rot[1, 1])
    defect = float(np.linalg.norm(rot@rot.T-np.eye(2)))
    rows = []
    for m in elements:
        ev, o = np.linalg.eigh(m)
        rows.append(dict(weak_over_strong=float(ev[0]/ev[1]),
                         strong_axis_deg=float(math.degrees(math.atan2(o[1, 1], o[0, 1])) % 180)))
    return theta, defect, rows, elements


def fluid_cycle(n):
    sys.path.insert(0, '../ph1')
    import ph1_real_fluid_holonomy as p
    H, _, _ = p.real_fluid('CarbonDioxide')
    T, rho = p.loop(n, 310.0, 350.0, 200.0, 700.0)
    return [H(float(t), float(r)) for t, r in zip(T, rho)], p


def run():
    import numpy as np
    out = dict(exact=exact_witness(), fluid='CarbonDioxide', cycle='PH1 Table 1', rows=[])
    reference = None
    for n in (6, 12, 24, 48, 96, 384, 8000):
        Hs, p = fluid_cycle(n)
        theta, defect, rows, elements = chain(Hs)
        D = np.diag([37.0, 0.004])
        theta_u, _, rows_u, elements_u = chain([D@h@D for h in Hs])
        unit_shift = max(float(np.linalg.norm(a-b)) for a, b in zip(elements, elements_u))
        if n == 8000:
            reference = p.holonomy(Hs)['theta_transport']
        out['rows'].append(dict(elements=n, rotation_rad=theta, rotation_deg=math.degrees(theta),
                                orthogonality_defect=defect, unit_change_rotation=theta_u,
                                unit_change_element_shift=unit_shift,
                                weakest_element_ratio=min(r['weak_over_strong'] for r in rows),
                                chain_weak_axis_transmission=float(np.prod([r['weak_over_strong'] for r in rows]))))
        if n == 12:
            out['twelve_element_recipe'] = rows
    out['transport_return_angle_PH1'] = reference
    return out


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PH2_RESULT.json', 'w'), indent=1)
    print('exact', res['exact']['tan_theta'])
    for r in res['rows']:
        print('N=%5d rotation=%+.6f rad (%+.3f deg) defect=%.1e unit=%+.6f elem-shift=%.1e weakest=%.4f chain=%.3e' % (
            r['elements'], r['rotation_rad'], r['rotation_deg'], r['orthogonality_defect'],
            r['unit_change_rotation'], r['unit_change_element_shift'], r['weakest_element_ratio'],
            r['chain_weak_axis_transmission']))
    print('PH1 transport angle', res['transport_return_angle_PH1'])
    for k, r in enumerate(res['twelve_element_recipe'], 1):
        print('  element %2d  axis %7.2f deg   weak/strong %.4f' % (k, r['strong_axis_deg'], r['weak_over_strong']))
