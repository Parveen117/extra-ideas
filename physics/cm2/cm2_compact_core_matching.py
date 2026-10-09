"""CM2: exact compact/core matching and weak internal gaps in all 8 centre sectors.

Analytical claims and min-max/domain arguments are in CM2_COMPACT_CORE_MATCHING.md.
This replay checks frozen PC1 source data, exact chart/IMS identities, and
outward uniform bounds. All verification arithmetic uses stdlib Fraction.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parent


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pc = load(PHYSICS/'pc1/pc1_pair_clusters.py', 'cm2_pc1')
cb = load(PHYSICS/'cb1/cb1_compact_bridge.py', 'cm2_cb1')
gc, cr, cm = pc.gc1, pc.cr1, pc.cm1
RHO = F(5692389, 10**6)             # downward-rounded PC1 full-space E1 floor
PI_UPPER = F(22, 7)
ZERO = (0, 0, 0)


@lru_cache(None)
def pair_data():
    record = json.loads((HERE/'CM2_PAIR_TRIAL.json').read_text())
    expected = sorted([(a, b, 0) for a in range(15) for b in range(8) if a+2*b <= 14],
                      key=lambda e: (cr.degree(e), e))
    basis = [tuple(e) for e in record['basis']]
    if record['d'] != 2 or record['degree'] != 14 or basis != expected or len(basis) != len(record['coefficients']):
        raise ValueError('invalid frozen pair trial basis')
    p = {e: F(v) for e, v in zip(basis, record['coefficients'])}
    width, shape = F(record['width']), F(record['shape'])
    if width <= 0:
        raise ValueError('positive Gaussian width required')
    mean = cr.free_means(2, width)
    square = cr.pmul(p, p)
    norm = mean(square)
    if norm <= 0:
        raise ValueError('positive trial norm required')
    hp = cm.scalar_h(p, 2, width)
    eta = mean(cr.pmul(p, hp))/norm
    second = mean(cr.pmul(hp, hp))/norm
    averaged_product = {ZERO: F(1), (1, 0, 0): shape,
                        (2, 0, 0): shape**2/8, (0, 1, 0): shape**2/2}
    norm_one = 1+3*shape/width+F(15, 4)*shape**2/width**2
    overlap = mean(cr.pmul(averaged_product, p))**2/(norm_one**2*norm)
    rho2 = (gc.E0+gc.E1)/2
    if not eta < rho2 or not 0 <= overlap <= 1:
        raise ValueError('pair trial not certified below the count')
    low = gc.floor_of_lowest(eta, second, rho2)
    share, true_overlap = pc.one_turn_share({'eta': eta, 'overlap2': overlap}, rho2, low)
    line = pc.cluster_line(3, rho2, low, share)
    return {'norm': norm, 'eta': eta, 'second': second, 'overlap': overlap,
            'low': low, 'share': share, 'true_overlap': true_overlap,
            'rho2': rho2, 'line3': line, 'p': p, 'hp': hp, 'width': width}


def character(parity, signs):
    if len(parity) != 3 or len(signs) != 3 or any(x not in (-1, 1) for x in parity+signs):
        raise ValueError('three signs in each character argument required')
    out = 1
    for sigma, sign in zip(parity, signs):
        if sign == -1:
            out *= sigma
    return out


def geometry(q, c):
    """Exact one-link inverse metric, half-density gradient, and scalar term."""
    q = F(q)
    y = sum(x*x for x in c)
    den = 1-q*y
    if q < 0 or len(c) != 3 or den <= 0:
        raise ValueError('point outside the compact coordinate chart')
    metric = [[F(i == j)-q*c[i]*c[j] for j in range(3)] for i in range(3)]
    a = [q*x/(2*den) for x in c]
    correction = 3*q/4+q*q*y/(8*den)
    return metric, a, correction


def localization(start, radius, width, level=RHO):
    """Constants for E_k(H) >= a*theta^(1/3)-b from an input core floor.

    The inner count is used only after verifying the exterior floor is above
    its threshold for the entire theta >= start interval.
    """
    start, radius, width, level = map(F, (start, radius, width, level))
    if min(start, radius, width, level) <= 0:
        raise ValueError('positive localization parameters required')
    reach = radius+width
    w_lower = cb.root_floor(start, 3)
    if w_lower <= reach**2:
        raise ValueError('localization reaches an equator')
    # Outside minus inner threshold, divided by 2w:
    # (R-level)+(A^2 level-21/4)/w. These conditions prove its sign for ALL w.
    if radius < level or reach**2*level < F(21, 4):
        raise ValueError('exterior has no certified uniform clearance')
    cost = PI_UPPER**2/(8*width**2)
    return {'start': start, 'radius': radius, 'width': width, 'reach': reach,
            'w_lower': w_lower, 'ims_scaled_upper': cost,
            'excitation_coefficient': 2*(level-cost),
            'subtracted_constant': 2*level*reach**2-F(9, 2)}


def weak_window(start, width, advertised, constant):
    local = localization(start, F(57, 10), width)
    # Same actual CB1 reading and complete upper error, with a smaller tail.
    lift = cb.compact_upper(start, F(6), F(1, 4), 28, F(1, 500000))
    coefficient = local['excitation_coefficient']-2*lift['upper']
    margin = advertised*local['w_lower']-constant
    if coefficient < advertised or local['subtracted_constant'] > constant or margin <= 0:
        raise ValueError('advertised interval gap is not certified')
    return {'theta_start': start, 'all_eight_fixed_centre_sectors': True,
            'advertised_coefficient': str(advertised), 'advertised_constant': str(constant),
            'threshold_gap_lower': str(cb.down(margin)),
            'available_coefficient_lower': str(cb.down(coefficient)),
            'actual_constant_upper': str(cb.up(local['subtracted_constant'])),
            'localization_radius': str(local['radius']), 'localization_width': str(width),
            'chart_reach': str(local['reach']), 'ims_scaled_upper': str(cb.up(local['ims_scaled_upper'])),
            'compact_ground_over_2_cuberoot_upper': str(cb.up(lift['upper'])),
            'trial_cutoff_radius': '6', 'trial_cutoff_width': '1/4', 'trial_tail_moment_order': 28,
            'trial_young_parameter': '1/500000', 'trial_tail_upper': str(cb.up(lift['tail'], 24))}


def run():
    checks = {}
    pair = pair_data()
    checks['PC1 pair input: exact frozen 64-reading trial has positive variance'] = (
        len(pair['p']) == 64 and pair['second'] > pair['eta']**2)
    checks['PC1 radial first and second signs replay with complete tails'] = gc.level_one(gc.E0) and gc.level_two(gc.E1)[0]
    checks['PC1 even angular and uneven odd floors replay with complete tails'] = all(
        pc.dr1.level_one(lam, m-1, int(lam)+7)
        for m, lam in ((5, F(33592, 10000)), (7, F(42461, 10000))))
    checks['full pair count: even radial floor is lower than every other required sector'] = (
        pair['rho2'] == F(3213, 1000) and F(33592, 10000) > pair['rho2']
        and (gc.E0+F(42461, 10000))/2 > pair['rho2'])
    checks['frozen pair ground and share meet the PC1 bounds'] = (
        pair['low'] > F(26592, 10000) and pair['eta'] < F(26594, 10000)
        and pair['share'] > F(9556, 10000))
    checks['PC1 full-space three-turn floor exceeds the retained 5.692389'] = pair['line3'] > RHO

    chart_ok, scalar_ok, metric_ok = True, True, True
    for q, c in ((F(1, 7), (F(1, 2), F(1, 3), F(2, 5))),
                 (F(1, 2), (F(1), F(1, 4), F(1, 5))), (F(0), (F(1), F(2), F(3)))):
        metric, a, correction = geometry(q, c)
        y = sum(x*x for x in c)
        ba = [sum(metric[i][j]*a[j] for j in range(3)) for i in range(3)]
        chart_ok &= ba == [q*x/2 for x in c]
        scalar_ok &= correction == 3*q/4+sum(a[i]*ba[i] for i in range(3))/2 and correction >= 3*q/4
        for v in ((F(1), F(-2), F(3)), c):
            form = sum(v[i]*metric[i][j]*v[j] for i in range(3) for j in range(3))
            metric_ok &= form >= (1-q*y)*sum(x*x for x in v)
    checks['exact half-density identity: B grad(log J)/2 = q c/2'] = chart_ok
    checks['exact positive geometric correction is 3q/4+q^2 r^2/(8(1-q r^2))'] = scalar_ok
    checks['inverse metric bound includes its negative radial correction'] = metric_ok

    parities = list(product((-1, 1), repeat=3))
    checks['all eight centre characters are orthogonal on the eight central patches'] = all(
        sum(character(s, e)*character(t, e) for e in parities) == (8 if s == t else 0)
        for s in parities for t in parities)
    checks['folded coordinate sign(a_i) u_i is unchanged by each centre flip'] = all(
        (-e)*(-x) == e*x for e in (-1, 1) for x in (F(-2, 7), F(3, 5)))
    checks['IMS quadratic partition has exact norm, cross cancellation, and unit derivative cost'] = all(
        ((1-t*t)/(1+t*t))**2+(2*t/(1+t*t))**2 == 1
        and ((1-t*t)/(1+t*t))*(-2*t/(1+t*t))+(2*t/(1+t*t))*((1-t*t)/(1+t*t)) == 0
        for t in (F(0), F(1, 3), F(2), F(7)))
    # x^4(1-x)^4 = (1+x^2)(x^6-4x^5+5x^4-4x^2+4)-4.
    dividend = {4: F(1), 5: F(-4), 6: F(6), 7: F(-4), 8: F(1)}
    quotient = {0: F(4), 2: F(-4), 4: F(5), 5: F(-4), 6: F(1)}
    recovered = {0: F(-4)}
    for k, v in quotient.items():
        for j in (0, 2):
            recovered[k+j] = recovered.get(k+j, F(0))+v
    checks['pi < 22/7 follows from the positive rational integral identity'] = (
        {k: v for k, v in recovered.items() if v} == dividend
        and sum(v/F(k+1) for k, v in quotient.items()) == PI_UPPER)
    checks['total vector radius has spherical gradient squared at most one'] = all(
        0 <= 1-sum(y*y for y in ys)/sum(ys) <= 1
        for ys in ((F(1, 9), F(1, 4), F(1, 5)), (F(0), F(1), F(0))))
    checks['quartic scaling gives epsilon^-2 = theta epsilon^4 = theta^(1/3)'] = all(
        (1/e**2) == (1/e**6)*e**4 for e in (F(1, 2), F(2, 9), F(1, 100)))
    p = cb.flat_trial()
    checks['CB1 retained ground reading has the exact full action and eta < 5.186744'] = p['action_check'] and p['eta'] < F(5186744, 10**6)
    limit = 2*(RHO-p['eta'])
    checks['fixed-sector gap limiting coefficient exceeds 1.0112'] = limit > F(10112, 10000)

    windows = [weak_window(2*10**9, F(31, 10), F(73, 100), F(878)),
               weak_window(10**10, F(37, 10), F(82, 100), F(1002)),
               weak_window(10**12, F(57, 10), F(93, 100), F(1476))]
    for row in windows:
        checks[f'uniform chart, exterior, IMS and trial reserve from theta={row["theta_start"]}'] = (
            F(row['threshold_gap_lower']) > 0
            and F(row['available_coefficient_lower']) >= F(row['advertised_coefficient']))
    return {'stage': 'CM2', 'checks': checks, 'all_pass': all(checks.values()),
            'frozen_pair': {'ground_lower': str(cb.down(pair['low'])),
                'ground_upper': str(cb.up(pair['eta'])), 'share_lower': str(cb.down(pair['share'])),
                'three_turn_line_lower': str(cb.down(pair['line3']))},
            'retained_full_core_E1_floor': str(RHO), 'flat_ground_trial_upper': str(cb.up(p['eta'])),
            'fixed_sector_gap_limit_lower': str(cb.down(limit)), 'weak_windows': windows,
            'claim_boundary': {'carrier': '(S3)^3 with Haar measure, simultaneous-conjugation invariance',
                'sectors': 'all 8 fixed independent centre characters',
                'spectral_limit': 'for each fixed k and sigma: E_k,sigma(theta)/(2 theta^(1/3)) -> e_k of the free physical core',
                'internal_sector_gap_limit': 'Delta_sigma(theta)/theta^(1/3) -> 2(e1-e0)',
                'gap_across_all_centre_sectors': 'Delta_all(theta)/theta^(1/3) -> 0; splitting rate not determined',
                'norm_resolvent_convergence': False, 'nonconstant_spatial_modes': False,
                'volume_uniformity': False, 'continuum_mass_gap': False,
                'proof_status': 'written analytic proof plus exact/outward replay; not proof-assistant verification'}}


def source_pins():
    paths = ['pc1/PC1_PAIR_CLUSTERS.md', 'pc1/pc1_pair_clusters.py', 'pc1/PC1_RESULT.json',
             'dr1/DR1_CORE_BY_ROWS.md', 'dr1/dr1_core_by_rows.py',
             'gc1/gc1_core_gap_count.py', 'cr1/cr1_core_rates.py', 'cm1/cm1_core_sector_memory.py',
             'tc1/tc1_core_of_turns.py', 'ol1/OL1_ONE_SITE_LATTICE.md', 'ol1/ol1_one_site_lattice.py',
             'cb1/CB1_COMPACT_BRIDGE.md', 'cb1/cb1_compact_bridge.py', 'cb1/CB1_TRIAL.json',
             'cm2/CM2_PAIR_TRIAL.json', 'cm2/CM2_COMPACT_CORE_MATCHING.md',
             'cm2/cm2_compact_core_matching.py', 'cm2/test_cm2.py']
    return {'physics/'+p: hashlib.sha256((PHYSICS/p).read_bytes()).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if not result['all_pass']:
        raise SystemExit('FAIL: '+', '.join(k for k, v in result['checks'].items() if not v))
    result['source_sha256'] = source_pins()
    path = HERE/'CM2_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale CM2 result or source pin')
    print('PASS', len(result['checks']), 'exact/outward checks')
    for row in result['weak_windows']:
        print('Every fixed centre sector, theta >=', row['theta_start'], ': gap >=',
              row['advertised_coefficient'], '* theta^(1/3) -', row['advertised_constant'])
    print('Common limiting internal gap coefficient >=', result['fixed_sector_gap_limit_lower'])


if __name__ == '__main__':
    main()
