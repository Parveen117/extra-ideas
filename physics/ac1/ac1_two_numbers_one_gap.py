"""AC1, second part: what stands between the charge number and the mass number.

AC1 said "no relation follows".  This part says exactly what is missing: the mass number is not one number,
it is the product of two mass counts; the charge number does not know the mass at all.  The only thing not
fixed is the mass count of each kind of matter.
"""
import json, math, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def run():
    out, num = {}, {}
    z = lambda e: sp.simplify(e) == 0
    G, c, m, J, e2, kap, m1, m2, r, al, x = sp.symbols('G c m J e2 kappa m1 m2 r alpha x', positive=True)

    # the three lengths of a turning charged centre (SW2, GM1): r_s, a, sqrt(q2)
    rs = 2*G*m/c**2; a = J/(m*c); q2 = G*e2/c**4

    # G1: the turn term r_s a does not know the mass; so the charge number does not know it (nor G)
    out['G1_turn_term_is_mass_free'] = sp.diff(sp.simplify(rs*a), m) == 0
    a_c = sp.simplify(q2/(rs*a))
    out['G1_charge_number_is_mass_free'] = sp.diff(a_c, m) == 0 and sp.diff(a_c, G) == 0 and z(a_c - e2/(2*J*c))

    # G2: the mass number of a pair is the product of two mass counts, n = r_s / (2 l), l^2 = G kappa / c^3
    l = sp.sqrt(G*kap/c**3)
    n = lambda mm: (2*G*mm/c**2)/(2*l)
    out['G2_mass_number_is_a_product_of_two_counts'] = z(n(m1)*n(m2) - G*m1*m2/(kap*c))

    # G3: for one centre with the turn of one reading (J = kappa/2): mass number = r_s/(4a);
    #     charge number / mass number = 4 q2 / r_s^2 = (charge length / half mass length)^2
    aJ = a.subs(J, kap/2)
    a_g = G*m**2/(kap*c)
    out['G3_mass_number_is_rs_over_4a'] = z(a_g - rs/(4*aJ))
    out['G3_ratio_is_two_lengths_squared'] = z(a_c.subs(J, kap/2)/a_g - 4*q2/rs**2)

    # G4: the fall reaches 1 where r^2 + a^2 = r_s r - q2 (GM1).  With x = r_s/a = 4 alpha_g and q2 = alpha r_s a:
    #     a root exists iff x^2 >= 4 + 4 alpha x, i.e. alpha_g >= (alpha + sqrt(1 + alpha^2))/2
    disc = sp.discriminant(r**2 - x*r + 1 + al*x, r)            # lengths in units of a
    out['G4_edge_condition'] = z(disc - (x**2 - 4 - 4*al*x))
    xe = 2*al + 2*sp.sqrt(1 + al**2)
    out['G4_edge_value'] = z((x**2 - 4 - 4*al*x).subs(x, xe)) and z(xe/4 - (al + sp.sqrt(1 + al**2))/2)

    # G5: numbers (constants recalled)
    hbar, cc, GG, me, mp_, e_, k_e = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.67262192369e-27, 1.602176634e-19, 8.9875517923e9
    mP = math.sqrt(hbar*cc/GG)
    ne, np_ = me/mP, mp_/mP
    num['mass_count_electron'] = ne; num['mass_count_proton'] = np_
    num['product'] = ne*np_; num['alpha_g_electron_proton'] = GG*me*mp_/(hbar*cc)
    num['alpha'] = k_e*e_**2/(hbar*cc)
    num['edge_mass_count_at_alpha'] = math.sqrt((num['alpha'] + math.sqrt(1 + num['alpha']**2))/2)
    out['G5_product_is_the_mass_number'] = abs(ne*np_/num['alpha_g_electron_proton'] - 1) < 1e-12
    out['G5_counts_are_tiny'] = ne < 1e-22 and np_ < 1e-19

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'AC1_TWO_NUMBERS_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k, v in o.items(): print(k, v)
    for k, v in n.items(): print(k, v)
