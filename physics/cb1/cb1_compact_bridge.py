"""CB1: compact heat-return bounds and a correlated flat-to-compact trial.

Certificate replay is stdlib Fraction arithmetic. Bounded positive kernels,
the closed-form lift, and their spectral interfaces are proved in the note.
No floating trial search or finite hidden inverse enters the verifier.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, isqrt
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


ol = load(PHYSICS/'ol1/ol1_one_site_lattice.py', 'cb1_ol1')
gc = ol.gc1
cm = gc.cm1
cr = cm.cr


def down(x, digits=12):
    scale = 10**digits
    return F((x*scale).numerator//(x*scale).denominator, scale)


def up(x, digits=12):
    return -down(-x, digits)


def root_floor(x, power, digits=12):
    """Integer-only enclosure of a nonnegative rational root."""
    if x < 0 or power < 1:
        raise ValueError('nonnegative radicand and positive power required')
    scale = 10**digits
    target = x.numerator*scale**power//x.denominator
    if power == 2:
        return F(isqrt(target), scale)
    low, high = 0, 1
    while high**power <= target:
        high *= 2
    while high-low > 1:
        mid = (low+high)//2
        if mid**power <= target:
            low = mid
        else:
            high = mid
    return F(low, scale)


def root_upper(x, power, digits=12):
    low = root_floor(x, power, digits)
    return low if low**power == x else low+F(1, 10**digits)


def exp_minus_upper(x, terms=40):
    """exp(-x) <= reciprocal of a positive Taylor partial sum for exp(x)."""
    if x < 0 or terms < 1:
        raise ValueError('positive Taylor majorant needs x >= 0, terms >= 1')
    return 1/sum((x**k/F(factorial(k)) for k in range(terms+1)), F(0))


def heat_deviation_upper(even=False):
    """Full infinite character tail, not just a finite heat-kernel sum."""
    time, step, first = (F(1, 2), 2, 2) if even else (F(1), 1, 1)
    kept = [first+step*k for k in range(4)]
    term = lambda n: (n+1)**2*exp_minus_upper(time*n*(n+2))
    n = kept[-1]+step
    # For every later n the dimension ratio falls and eigenvalue increment rises.
    ratio = F((n+step+1)**2, (n+1)**2)*exp_minus_upper(time*step*(2*n+step+2))
    if not 0 < ratio < 1:
        raise ValueError('uncertified heat tail ratio')
    tail = term(n)/(1-ratio)
    total = sum((term(k) for k in kept), F(0))+tail
    return {'time': time, 'deviation_upper': total, 'tail_upper': tail, 'ratio_upper': ratio}


def return_reserve(links, deviation):
    if links < 1 or not 0 <= deviation < 1:
        raise ValueError('positive number of links and heat deviation below one required')
    return ((1-deviation)/(1+deviation))**links


def sin_ratio_lower(x):
    """sin(x)/x lower on 0 < x <= 1 by the degree-seven alternating sum."""
    if not 0 < x <= 1:
        raise ValueError('sine certificate is restricted to 0 < angle <= 1')
    return 1-x*x/6+x**4/120-x**6/5040


@lru_cache(None)
def flat_trial():
    record = json.loads((HERE/'CB1_TRIAL.json').read_text())['trial']
    basis = [tuple(e) for e in record['basis']]
    if basis != cm.inv_basis(record['degree']) or len(basis) != len(record['coefficients']):
        raise ValueError('frozen Gaussian basis mismatch')
    p = {e: F(c) for e, c in zip(basis, record['coefficients'])}
    width = F(record['width'])
    mean = cr.free_means(3, width)
    square = cr.pmul(p, p)
    norm = mean(square)
    if norm <= 0:
        raise ValueError('positive norm required')
    hp = cm.scalar_h(p, 3, width)
    kinetic_poly = cr.padd(hp, cr.pmul(cm.E2, p), -1)
    r2 = lambda poly: {(a+1, b, c): v for (a, b, c), v in poly.items()}
    potential = mean(cr.pmul(cm.E2, square))/norm
    kinetic = mean(cr.pmul(p, kinetic_poly))/norm
    eta = kinetic+potential
    # Integration by parts in R^9: weighted kinetic density = <psi,r^2 K psi>+9/2.
    weighted_kinetic = mean(r2(cr.pmul(p, kinetic_poly)))/norm+F(9, 2)
    weighted_potential = mean(r2(cr.pmul(cm.E2, square)))/norm
    return {'p': p, 'width': width, 'mean': mean, 'square': square, 'norm': norm,
            'eta': eta, 'kinetic': kinetic, 'potential': potential,
            'weighted_kinetic': weighted_kinetic,
            'weighted_energy': weighted_kinetic+weighted_potential,
            'action_check': eta == mean(cr.pmul(p, hp))/norm}


def radial_moment(order):
    if order < 0:
        raise ValueError('nonnegative moment order required')
    p = flat_trial()
    return p['mean']({(a+order, b, c): v for (a, b, c), v in p['square'].items()})/p['norm']


def compact_upper(start, radius, width, order, young):
    """E0(H_theta)/(2 theta^(1/3)) upper, uniformly for theta >= start.

    chi=1 up to radius, linear to zero at radius+width. All derivative, tail,
    Haar density and inverse-metric errors are kept as in CB1 equations (8)-(10).
    """
    if min(start, radius, width, young) <= 0 or order < 1:
        raise ValueError('positive lift parameters and moment order required')
    p = flat_trial()
    tail = radial_moment(order)/radius**(2*order)
    scale = root_upper(F(1, start), 3)
    reach = radius+width
    residual = 1-scale*reach**2
    if not 0 <= tail < 1 or residual <= 0:
        raise ValueError('cutoff or compact-chart denominator not certified')
    density = 1/root_floor(residual, 2)
    coefficient = density**3/2
    energy = (p['eta']+young*p['kinetic']
              +coefficient*scale*(p['weighted_energy']+young*p['weighted_kinetic'])
              +(1+1/young)*density*tail/(2*width**2))/(1-tail)
    return {'upper': energy, 'tail': tail, 'scale_upper': scale,
            'density_upper': density, 'metric_remainder_coefficient': coefficient,
            'radius': radius, 'width': width, 'moment_order': order, 'young': young}


def weak_row(start, angle, radius, width, order, young, advertised):
    sigma = sin_ratio_lower(angle)
    reach_cube = 4*root_floor(F(start), 2)*angle**3*sigma
    if reach_cube < 9**3:
        raise ValueError('compact comparison has not reached both sign certificates')
    count = min(2*gc.E0+gc.E1, gc.E0+2*gc.P0)
    line = count*root_floor(2*sigma*sigma, 3)
    lift = compact_upper(start, radius, width, order, young)
    coefficient = line-2*lift['upper']
    margin = advertised*root_floor(F(start), 3)-F(15, 2)
    if coefficient < advertised or margin <= 0:
        raise ValueError('advertised uniform gap is not certified')
    return {'start': start, 'angle': str(angle), 'sigma_lower': str(down(sigma)),
            'comparison_reach_cube_lower': str(down(reach_cube)),
            'comparison_coefficient_lower': str(down(line)),
            'compact_ground_coefficient_upper': str(up(2*lift['upper'])),
            'available_gap_coefficient_lower': str(down(coefficient)),
            'advertised_gap_coefficient': str(advertised), 'subtracted_constant': '15/2',
            'minimum_gap_lower': str(down(margin)),
            'cutoff_radius': str(radius), 'cutoff_width': str(width), 'moment_order': order,
            'young_parameter': str(young), 'tail_mass_upper': str(up(lift['tail'], 18)),
            'density_upper': str(up(lift['density_upper']))}


def run():
    checks = {}
    full, even = heat_deviation_upper(), heat_deviation_upper(True)
    checks['full SU2 heat kernel at t=1 stays within 1/4 of Haar, with infinite tail'] = full['deviation_upper'] < F(1, 4)
    checks['centre quotient heat kernel at t=1/2 stays within 1/5 of Haar, with infinite tail'] = even['deviation_upper'] < F(1, 5)
    checks['vacuum-sector return reserve coefficient is 8/27'] = return_reserve(3, F(1, 5)) == F(8, 27)
    checks['all-sector return reserve coefficient is 27/125'] = return_reserve(3, F(1, 4)) == F(27, 125)
    checks['vacuum-sector gap coefficient 2 beta/t is 32/27'] = 2*return_reserve(3, F(1, 5))/even['time'] == F(32, 27)
    checks['all-sector gap coefficient 2 beta/t is 54/125'] = 2*return_reserve(3, F(1, 4))/full['time'] == F(54, 125)
    checks['finite-volume return reserve loses a factor 3/5 per additional link'] = all(
        return_reserve(n+1, F(1, 4)) == F(3, 5)*return_reserve(n, F(1, 4)) for n in range(1, 20))
    # Algebraic sharp row-overlap bound (r is the square root of the cross-ratio bound).
    checks['likelihood-ratio contraction remainder is the square (r*a-1)^2'] = all(
        (r-1)/(r+1)-(r*r*a-1)*(1-a)/(a*(r*r-1)) == (r*a-1)**2/(a*(r*r-1))
        for r in (F(3, 2), F(2), F(7)) for a in (F(1)/(r*r), F(1)/r, F(1)))
    checks['GC1 first-level sign certificate replayed'] = gc.level_one(gc.E0)
    checks['GC1 second-level sign certificate replayed'] = gc.level_two(gc.E1)[0]
    checks['GC1 nonzero-angular floor and selected comparison channel'] = gc.P0**3 <= F(2187, 64) and 2*gc.E0+gc.E1 < gc.E0+2*gc.P0
    p = flat_trial()
    checks['frozen OM1 trial has 67 scalar coefficients and exact action identity'] = len(p['p']) == 67 and p['action_check']
    checks['kinetic and weighted energy densities are positive'] = min(p['kinetic'], p['potential'], p['weighted_kinetic'], p['weighted_energy']) > 0
    checks['correlated free-core reading stays below 5.186743367'] = p['eta'] < F(5186743367, 10**9)
    rows = [weak_row(2*10**6, F(513, 1000), F(19, 4), F(1, 4), 18, F(442, 10**6), F(11, 100)),
            weak_row(10**7, F(39, 100), F(5), F(1, 4), 20, F(1, 6000), F(36, 100))]
    for row in rows:
        checks[f'full compact-lift error below the weak gap reserve from theta={row["start"]}'] = F(row['available_gap_coefficient_lower']) > F(row['advertised_gap_coefficient']) and F(row['minimum_gap_lower']) > 0
    limit = (2*gc.E0+gc.E1)*root_floor(F(2), 3)-2*p['eta']
    checks['correlated compact lift doubles OL1 liminf lower from 0.327 to above 0.6685'] = limit > F(6685, 10000)
    # For total radius support, product(1-a_i) >= 1-sum(a_i), not an unproved flat Haar measure.
    checks['three-link compact Haar product domination on nonnegative squared radii'] = all(
        (1-a)*(1-b)*(1-c) >= 1-a-b-c
        for a, b, c in ((F(1, 10), F(1, 5), F(1, 4)), (F(0), F(1, 2), F(0)), (F(1, 7),)*3))
    return {'schema': 'cb1_compact_bridge_v1', 'arithmetic': 'Fraction only; frozen trial', 'checks': checks,
        'heat_kernel': {'full_deviation_upper': str(up(full['deviation_upper'])),
                        'full_tail_upper': str(up(full['tail_upper'], 24)),
                        'even_deviation_upper': str(up(even['deviation_upper'])),
                        'even_tail_upper': str(up(even['tail_upper'], 36))},
        'all_coupling_bounds': {'vacuum_sector_H': '(32/27)*exp(-3*theta)',
                               'all_gauge_sectors_H': '(54/125)*exp(-6*theta)',
                               'finite_lattice_H': '2*(3/5)^links*exp(-2*theta*plaquettes)',
                               'vacuum_sector_A': '(8/27)*exp(-12*theta_YM)',
                               'all_gauge_sectors_A': '(27/250)*exp(-24*theta_YM)'},
        'flat_trial': {'eta_upper': str(up(p['eta'])), 'kinetic_upper': str(up(p['kinetic'])),
                       'weighted_kinetic_upper': str(up(p['weighted_kinetic'])),
                       'weighted_energy_upper': str(up(p['weighted_energy']))},
        'weak_windows': rows, 'asymptotic_gap_liminf_lower': str(down(limit)),
        'claim_boundary': {'one_site_all_finite_couplings': 'EXPLICIT_POSITIVE_GAP',
            'flux_sector_bottom_above_vacuum': 'BOUNDED_FROM_BELOW_AT_EACH_FINITE_COUPLING',
            'within_each_flux_sector_excitation_gap': 'NOT_COMPUTED',
            'weak_vacuum_gap_scale': 'THETA_TO_ONE_THIRD',
            'finite_lattice_all_couplings': 'EXPLICIT_BUT_DECAYS_WITH_SIZE',
            'volume_uniformity': False, 'continuum_mass_gap': False,
            'flat_to_compact_full_spectral_convergence': 'NOT_CLAIMED_ONLY_TRIAL_LIFT'}}


def source_pins():
    paths = ['ol1/OL1_ONE_SITE_LATTICE.md', 'ol1/ol1_one_site_lattice.py', 'ol1/OL1_RESULT.json',
             'gc1/gc1_core_gap_count.py', 'cm1/cm1_core_sector_memory.py', 'cr1/cr1_core_rates.py',
             'tc1/tc1_core_of_turns.py', 'cb1/CB1_TRIAL.json', 'cb1/cb1_compact_bridge.py',
             'cb1/CB1_COMPACT_BRIDGE.md', 'cb1/test_cb1.py']
    return {'physics/'+p: hashlib.sha256((PHYSICS/p).read_bytes()).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if not all(result['checks'].values()):
        raise SystemExit('FAIL: '+', '.join(k for k, v in result['checks'].items() if not v))
    result['source_sha256'] = source_pins()
    path = HERE/'CB1_RESULT.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    if args.check and json.loads(path.read_text()) != result:
        raise SystemExit('FAIL: stale CB1 result or source pin')
    print('PASS', len(result['checks']), 'exact/outward checks')
    for row in result['weak_windows']:
        print('theta >=', row['start'], 'gap >=', row['advertised_gap_coefficient'],
              '* theta^(1/3) - 15/2; threshold margin >=', row['minimum_gap_lower'])
    print('asymptotic liminf >=', result['asymptotic_gap_liminf_lower'])


if __name__ == '__main__':
    main()
