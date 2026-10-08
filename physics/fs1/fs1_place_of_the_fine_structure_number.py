"""FS1: where the number 1/137 sits in the line, exactly; and what does not give it.

GM1's centre has three constants: r_s (mass), a (displacement, turning), q2 (charge).  Memory + iota x swirl potential,
far away:   r_s/d + ( iota r_s a cos - q2 ) / d^2 + ...
Symbolic (sympy) and plain numbers.  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, num = {}, {}
    G, M, J, c, ke, e, hb = sp.symbols('G M J c k_e e hbar', positive=True)
    rs = 2*G*M/c**2
    a = J/(M*c)
    q2 = G*ke*e**2/c**4                                              # EG1 / GM1: the charge term of the memory
    # A1: the one pure number of the centre in which the unit of gravity cancels
    ratio = sp.simplify(q2/(rs*a))
    out['A1_gravity_unit_cancels'] = not ratio.has(G) and not ratio.has(M) and z(ratio - ke*e**2/(2*J*c))
    out['A1_is_alpha_for_half_unit_turning'] = z(ratio.subs(J, hb/2) - ke*e**2/(hb*c))
    # A2: far field of the complex memory of the charged displaced centre, on the axis: coefficient of 1/d^2 is r_s a (iota - alpha)
    d, al = sp.symbols('d alpha', positive=True)
    rs_, a_ = sp.symbols('r_s a', positive=True)
    R = d - sp.I*a_                                                  # on the axis, R = d - iota a
    Phi = rs_/R - al*rs_*a_/(R*sp.conjugate(R))                      # Re part of the first term is the memory; q2 = alpha r_s a
    ser = sp.series(Phi, d, sp.oo, 3).removeO()
    co2 = sp.simplify(sp.expand(ser).coeff(d, -2))
    out['A2_second_coefficient'] = z(co2 - rs_*a_*(sp.I - al))
    # numbers: electron (same constants as GM1)
    hbar, cc, Gn, me, ee, k_e = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.602176634e-19, 8.9875517923e9
    rs_e = 2*Gn*me/cc**2; a_e_len = (hbar/2)/(me*cc); q2_e = Gn*k_e*ee**2/cc**4
    alpha = q2_e/(rs_e*a_e_len)
    num['alpha_from_the_three_constants'] = alpha
    num['one_over_alpha'] = 1/alpha
    out['A1_number'] = abs(1/alpha - 137.036) < 0.001
    # A3: the left-over turn.  Moment ratio g: the direction of the moment goes round g/2 times per circuit in a magnetic field
    #     (relation recalled, not derived in the line).  g = 2 closes (GM1); measured g/2 - 1 leaves 2 pi (g/2 - 1) per circuit.
    anomaly = 0.00115965218                                         # recalled; not re-read at source
    left_over = 2*sp.pi*anomaly
    num['left_over_turn_per_circuit_rad'] = float(left_over)
    num['left_over_over_alpha'] = float(left_over/alpha)
    out['A3_left_over_turn_is_alpha_to_two_parts_in_1000'] = abs(float(left_over/alpha) - 1) < 2e-3
    # A4: what does not give it.
    #  (i) the diagonal on this centre: 'charge term = swirl term' would be alpha = 1
    #  (ii) ratios of the counts met so far (cuts 3, planes 3 or 6, d - 1, quarter turns) with one full turn: a search of small forms
    target = 1/alpha
    best = None
    for p in range(-3, 4):
        for q_ in range(-3, 4):
            for n in range(1, 13):
                for mden in range(1, 13):
                    val = (n/mden)*float(sp.pi)**p*2.0**(q_/2)
                    if val > 0:
                        err = abs(val - target)/target
                        if best is None or err < best[0]:
                            best = (err, n, mden, p, q_)
    num['closest_small_form'] = dict(relative_error=best[0], n=best[1], m=best[2], power_of_pi=best[3], half_power_of_2=best[4])
    out['A4_no_small_form_within_measurement'] = best[0] > 1e-6      # measured to parts in 10^10
    known = 4*float(sp.pi)**3 + float(sp.pi)**2 + float(sp.pi)
    num['known_near_miss_4pi3_pi2_pi'] = dict(value=known, relative_difference=(known - target)/target)
    out['A4_known_near_miss_is_refused'] = abs(known - target)/target > 1e-7
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'FS1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o); print(n)
