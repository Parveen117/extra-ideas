"""YC6: vacuum-sector enclosures and an exact, typed coupling-step audit."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'physics/yc5'))
import yc5_resolved_return as m


def mesh():
    return [F(n, 8) for n in range(113)]


def validate_endpoint(record):
    theta = F(record['theta'])
    if not 0 <= theta <= 14:
        raise ValueError('vacuum certificate window is [0,14]')
    if len(record['levels']) != 2 or any(len(pair) != 2 for pair in record['levels']):
        raise ValueError('both vacuum levels need lower and upper endpoints')
    kind = record['excited_upper_method']
    if kind not in ('schur', 'ritz'):
        raise ValueError('unknown excited upper method')
    for index, pair in enumerate(record['levels']):
        lower, upper = map(F, pair)
        if not 8*index <= lower <= upper:
            raise AssertionError('invalid ordered vacuum endpoints')
        if m.counts(0, theta, lower, 'lower', True)[0] > index:
            raise AssertionError(('failed lower count', theta, index))
        if index == 1 and kind == 'ritz':
            negative, zero, _ = m.c.counts(0, theta, upper, False)
        else:
            negative, zero, _ = m.counts(0, theta, upper, 'upper', True)
        if negative+zero < index+1:
            raise AssertionError(('failed upper count', theta, index))


def envelopes(records):
    if [F(r['theta']) for r in records] != mesh():
        raise AssertionError('all 113 eighth-mesh nodes from 0 through 14 required')
    floor = F(8)
    source_theta = F(0)
    rows = []
    for r in records:
        raw = F(r['levels'][1][0])
        if raw > floor:
            floor, source_theta = raw, F(r['theta'])
        rows.append({'theta': r['theta'], 'second_floor': str(floor),
                     'floor_source_theta': str(source_theta),
                     'gap': [str(floor-F(r['levels'][0][1])),
                             str(F(r['levels'][1][1])-F(r['levels'][0][0]))]})
    return rows


def coupling_cells(records, rows):
    return [{'left': a['theta'], 'right': b['theta'],
             'gap_lower': str(F(row['second_floor'])-F(b['levels'][0][1]))}
            for a, b, row in zip(records, records[1:], rows)]


def step_rows(rows):
    by_theta = {F(row['theta']): list(map(F, row['gap'])) for row in rows}
    output = []
    for theta in (F(1, 2), F(1), F(3, 2), F(2), F(5, 2), F(3), F(7, 2)):
        low, high = by_theta[theta]
        next_low, next_high = by_theta[4*theta]
        if low <= 0:
            raise AssertionError('positive denominator gap enclosure required')
        a, b = next_low/high, next_high/low
        output.append({'theta': str(theta), 'next_theta': str(4*theta),
                       'gap_ratio': [str(a), str(b)],
                       'normalized_ratio_cubed': [str(a**3/4), str(b**3/4)]})
    return output


def transpose(a):
    return list(map(list, zip(*a)))


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def add(a, b, scale=F(1)):
    return [[F(x)+scale*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def inverse(a):
    n = len(a)
    rows = [[F(x) for x in row]+[F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        pivot = next((j for j in range(k, n) if rows[j][k]), None)
        if pivot is None:
            raise ValueError('singular matrix')
        rows[k], rows[pivot] = rows[pivot], rows[k]
        value = rows[k][k]
        rows[k] = [x/value for x in rows[k]]
        for i in range(n):
            if i != k:
                value = rows[i][k]
                rows[i] = [x-value*y for x, y in zip(rows[i], rows[k])]
    return [row[n:] for row in rows]


def noncommuting_step_control():
    d0, w = [[F(5), F(1)], [F(1), F(8)]], [[F(2), F(1)], [F(1), F(3)]]
    source, eye = [[F(1), F(2)], [F(-1), F(1)]], [[F(1), F(0)], [F(0), F(1)]]
    theta, next_theta, z = F(1), F(4), F(1)
    r = inverse(add(add(d0, w, theta), eye, -z))
    next_r = inverse(add(add(d0, w, next_theta), eye, -z))
    sigma = mm(transpose(source), mm(r, source))
    next_sigma = mm(transpose(source), mm(next_r, source))
    correction = mm(transpose(source), mm(next_r, mm(w, mm(r, source))))
    identity = next_sigma == add(sigma, correction, -(next_theta-theta))
    # At the base point, the Taylor coefficients are (-1)^n r*(R W)^n R r.
    power = r
    coefficients = []
    for n in range(5):
        coefficients.append(mm(transpose(source), mm(power, source)))
        power = mm(r, mm(w, power))
    return {'noncommuting': mm(d0, w) != mm(w, d0), 'resolvent_step_exact': identity,
            'return_decreases': m.c.inertia(add(sigma, next_sigma, -1))[0] == 0,
            'alternating_derivative_coefficients_positive': all(m.c.inertia(a)[0] == 0 for a in coefficients),
            'sigma_before': sigma, 'sigma_after': next_sigma, 'unsigned_taylor': coefficients}


def scale_controls():
    # TC1's C=g^(2/3)c makes both kinetic and quartic coefficients g^(2/3).
    dilation = F(2, 3)
    theta_growth = F(1, 3)  # CM2 fixed-sector gap exponent
    physical_g_power = 2-4*theta_growth
    return {'core_kinetic_power': 2-2*dilation,
            'core_quartic_power': -2+4*dilation,
            'compact_physical_g_power': physical_g_power,
            # Along a=exp(-t), g=t^-1/2: log(g^p/a)=t-(p/2)log(t).
            'one_site_log_power': -physical_g_power/2,
            'normalized_gap_counterexamples': [F(-1), F(0), F(1, 6)],
            'refinement_example': {'inverse_g2': F(5), 'increment': F(1),
                                   'theta_before': F(25), 'theta_after': F(36),
                                   'wrong_quadrupled_theta': F(100)}}


def source_pins():
    paths = ['physics/yc5/YC5_RESOLVED_RETURN.md', 'physics/yc5/yc5_resolved_return.py',
             'physics/yc5/YC5_RESULT.json', 'physics/yc4/yc4_harmonic_return.py',
             'physics/yc3/yc3_sector_splitting.py', 'physics/cm2/CM2_COMPACT_CORE_MATCHING.md',
             'physics/cm2/CM2_RESULT.json', 'physics/tc1/TC1_CORE_OF_TURNS.md',
             'physics/ol1/OL1_ONE_SITE_LATTICE.md', 'physics/yc6/YC6_VACUUM_STEP.md',
             'physics/yc6/yc6_vacuum_step.py', 'physics/yc6/test_yc6.py',
             'physics/yc6/YC6_ENDPOINTS.json']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def run():
    records = json.loads((HERE/'YC6_ENDPOINTS.json').read_text())['endpoints']
    rows = envelopes(records)
    prior = json.loads((ROOT/'physics/yc5/YC5_RESULT.json').read_text())
    checks = {'YC5 frozen source pins unchanged': m.source_pins() == prior['source_sha256']}
    for record in records:
        validate_endpoint(record)
    checks['all 113 endpoints pass both vacuum level counts'] = True
    cells = coupling_cells(records, rows)
    checks['all 112 complete coupling cells have vacuum gap at least one'] = all(F(c['gap_lower']) >= 1 for c in cells)
    checks['past certified floors only carried forward'] = all(F(r['floor_source_theta']) <= F(r['theta']) for r in rows)
    checks['retained vacuum cut remains complete with rank 13 and floor 24'] = (
        sum(len(b['sources']) for b in m.resolved_blocks(0).values()) == 13
        and all(b['base']['floor'] == 24 for b in m.resolved_blocks(0).values()))
    steps = step_rows(rows)
    checks['all seven fourfold-coupling ratios have positive rational endpoints'] = all(
        0 < F(r['gap_ratio'][0]) <= F(r['gap_ratio'][1]) for r in steps)
    control = noncommuting_step_control()
    for key in ('noncommuting', 'resolvent_step_exact', 'return_decreases', 'alternating_derivative_coefficients_positive'):
        checks['operator step control: '+key] = control[key]
    scale = scale_controls()
    checks['TC1 kinetic and potential have the same g power'] = scale['core_kinetic_power'] == scale['core_quartic_power'] == F(2, 3)
    checks['CM2 plus a declared lattice energy unit has g power 2/3'] = scale['compact_physical_g_power'] == F(2, 3)
    checks['normalized vanishing admits shrinking constant and growing absolute examples'] = all(p < F(1, 3) for p in scale['normalized_gap_counterexamples'])
    checks['fixed refinement is not a quadrupling of theta in the one-loop control'] = scale['refinement_example']['theta_after'] != scale['refinement_example']['wrong_quadrupled_theta']
    checks['TC1 one-site rate has inverse one-third logarithmic factor'] = scale['one_site_log_power'] == F(-1, 3)
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError([k for k, v in checks.items() if not v])
    selected = {'0', '2', '4', '6', '8', '10', '12', '14'}
    return {'stage': 'YC6', 'all_pass': True, 'checks': checks, 'coupling_window': '[0,14]',
            'sector': 'simultaneous gauge invariant, all three centre characters even',
            'uniform_vacuum_gap_lower': '1',
            'minimum_mesh_gap_lower': str(min(F(c['gap_lower']) for c in cells)),
            'point_gap_enclosures': {r['theta']: r['gap'] for r in rows if r['theta'] in selected},
            'coupling_cells': cells, 'monotone_second_floors': rows,
            'fourfold_coupling_steps': steps, 'step_control': control, 'scale_controls': scale,
            'claim_boundary': {'one_site_only': True, 'vacuum_sector_only': True,
                               'full_hidden_return_retained': True, 'ground_state_transform_analytic_identity': True,
                               'coupling_step_is_spatial_rg': False, 'all_sector_absolute_collapse_proved': False,
                               'information_chart_implies_volume_uniformity': False,
                               'logarithmic_running_derived_from_core': False,
                               'bare_coupling_identified_with_fixed_box_effective_coupling': False,
                               'constant_modes_excluded_from_full_theory': False,
                               'volume_uniform_or_4d_gap': False, 'floating_point_used_to_certify': False},
            'source_sha256': source_pins()}


def jsonable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = jsonable(run())
    path = HERE/'YC6_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale packet or source pins')
    print('PASS', len(result['checks']), 'exact checks; vacuum gap >=1 on [0,14]')
    print('113 endpoints, 112 coupling cells; point gaps:', result['point_gap_enclosures'])
