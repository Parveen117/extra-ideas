"""MG1: is the gap of the Yang-Mills line a mass count of the physics line?

AC1 (later note): the one missing rule is what fixes the mass count n = r_s/2l of a kind of matter.
The Yang-Mills line of the framework certifies gaps.  This stage writes a gap as a mass count, lists
what the certified gaps give, and says exactly where the proton's small count would have to come from.
sympy for the algebra, mpmath for the numbers.  Python 3.12."""
import json, math, os
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def run():
    out, num = {}, {}
    z = lambda e: sp.simplify(e) == 0
    G, c, m, kap, gam, ac = sp.symbols('G c m kappa gamma a_c', positive=True)

    # M1: the lengths of a centre with the turn of one reading (J = kappa/2)
    rs = 2*G*m/c**2; a = (kap/2)/(m*c); l = sp.sqrt(G*kap/c**3); g = m*c/kap     # g: rest rate per length (BR1)
    n = rs/(2*l)
    out['M1_mass_length_times_turn_length_is_the_record_area'] = z(rs*a - l**2)
    out['M1_mass_count_is_rest_rate_times_record_length'] = z(n - g*l) and z(n - l/(2*a))

    # M2: a gap gamma per cell of length a_c is the rest rate of the lightest record: g = gamma/a_c, so n = (l/a_c) gamma
    out['M2_gap_as_a_mass_count'] = z((gam/ac)*l - (l/ac)*gam)

    # M3: what the Yang-Mills line has certified, in cell units
    mp.mp.dps = 30
    d_red = -mp.log(mp.besseli(2, 2)/mp.besseli(1, 2))                           # YM-1
    out['M3_YM1_reduced_gap'] = abs(d_red - mp.mpf('0.83672330623')) < mp.mpf('1e-10')
    certified = {'YM-1 theta graph (kappa = 2)': float(d_red), 'YM-9/22 free rate, SU(2)': 0.75, 'YM95/98 free rate, SU(3)': 4/3,
                 'YM93 window, SU(2)': 16/125, 'YM97 window, SU(3)': 2/5, 'YM98 window end, SU(2)': 0.37, 'YM98 window end, SU(3)': 0.78}
    num['certified_rates_in_cell_units'] = certified
    out['M3_all_certified_rates_are_of_order_one'] = all(0.1 < v < 2 for v in certified.values())

    # M4: the one exact rate known at every coupling: the chain, -Log(I2(k)/I1(k)) (YM-21).  It falls like 3/(2k).
    rate = lambda k: -mp.log(mp.besseli(2, k)/mp.besseli(1, k))
    out['M4_chain_rate_falls_as_a_power'] = all(abs(k*rate(k) - mp.mpf(3)/2) < mp.mpf(3)/k for k in (mp.mpf(10)**3, mp.mpf(10)**5, mp.mpf(10)**7))
    out['M4_chain_rate_is_decreasing'] = all(rate(k) > rate(2*k) for k in (mp.mpf(1)/8, 1, 8, 64, 512))

    # M5: numbers (constants recalled)
    hbar, cc, GG, me, mp_ = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.67262192369e-27
    mP = math.sqrt(hbar*cc/GG); lP = math.sqrt(hbar*GG/cc**3)
    n_p, n_e = mp_/mP, me/mP
    num['mass_count_proton'] = n_p; num['mass_count_electron'] = n_e
    num['log_inverse_count_proton'] = math.log(1/n_p); num['log_inverse_count_electron'] = math.log(1/n_e)
    num['chain_coupling_needed_for_the_proton_count'] = 1.5/n_p
    num['cell_needed_if_rate_is_4_over_3_metres'] = lP*(4/3)/n_p
    out['M5_order_one_rate_needs_a_cell_of_nuclear_size'] = 1e-16 < num['cell_needed_if_rate_is_4_over_3_metres'] < 1e-15
    out['M5_power_law_needs_a_large_number'] = num['chain_coupling_needed_for_the_proton_count'] > 1e19
    out['M5_an_exponential_law_needs_a_number_near_44'] = 43.5 < num['log_inverse_count_proton'] < 44.5 and 51 < num['log_inverse_count_electron'] < 52

    # known physics, not the line: leading-order running for three colours without matter, exponent 1/(2 b0 g^2), b0 = 11/(16 pi^2)
    b0 = 11/(16*math.pi**2); g2 = 1/(2*b0*num['log_inverse_count_proton'])
    num['known_physics_not_the_line'] = {'g_squared_at_the_record_length': g2, 'theta_of_the_YM_line': 6/g2**2}

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'MG1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k, v in o.items(): print(k, v)
    for k, v in n.items(): print(k, v)
