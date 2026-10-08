"""AC1: the largest attraction that still has a first circle; and where a free circle of a mass is on the diagonal.

BR1: a charge on a circle of whole area l kappa has speed alpha/l.  MO1 / OR1: circles of a mass: L^2 = r_s r/(2 - 3x) per g^2.
OA1-T3: count factor^2 = 1 - 3x/2.  DO1: the diagonal is seen = lost.
Symbolic (sympy).  Python 3.12."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
z = lambda e: sp.simplify(e) == 0


def run():
    out, num = {}, {}
    al, l, x, d = sp.symbols('alpha l x d', positive=True)
    # A1: charge: speed alpha/l < 1: a circle of area l kappa exists iff alpha < l ; at alpha -> l its radius -> 0 and count factor -> 0
    radius = l**2*sp.sqrt(1 - al**2/l**2)/al                          # in units kappa/g
    out['A1_charge_limit'] = z(radius.subs(al, l)) and z((1 - al**2/l**2).subs(al, l))
    # A2: mass: (L / g r_s)^2 = 1/(x (2 - 3x)) : least over all circles at x = 1/3: L_min = sqrt 3 g r_s
    L2 = 1/(x*(2 - 3*x))
    xs = sp.solve(sp.diff(L2, x), x)
    out['A2_least_area_at_last_stable_circle'] = xs == [sp.Rational(1, 3)] and z(L2.subs(x, sp.Rational(1, 3)) - 3)
    # with alpha_g := g r_s / (2 kappa): a circle of area l kappa exists iff alpha_g <= l / (2 sqrt 3)
    g_, rs_, kap_ = sp.symbols('g r_s kappa', positive=True)
    grs = sp.solve(sp.Eq(sp.sqrt(3)*g_*rs_, l*kap_), rs_)[0]*g_          # g r_s at which the least area equals l kappa
    out['A2_mass_limit'] = z(grs/(2*kap_) - l/(2*sp.sqrt(3)))
    num['largest_alpha_charge_first_circle'] = 1.0
    num['largest_alpha_mass_first_circle'] = float(1/(2*sp.sqrt(3)))
    # A3: on a free circle of a mass count factor^2 = 1 - 3x/2 : at the last stable circle it is exactly 1/2: seen = lost
    cf2 = 1 - sp.Rational(3, 2)*x
    out['A3_last_stable_circle_is_on_the_diagonal'] = z(cf2.subs(x, sp.Rational(1, 3)) - sp.Rational(1, 2))
    # in d cuts: last stable circle at m = (4-d)/d (CD1), count factor^2 = 1 - (d/2) m = (d-2)/2 : one half only for d = 3
    m_last = (4 - d)/d
    out['A3_only_in_three_cuts'] = z((1 - d*m_last/2) - (d - 2)/2) and sp.solve(sp.Eq((d - 2)/2, sp.Rational(1, 2)), d) == [3]
    # three places of the eight-mark value 1/sqrt 2 in the field of a centre at rest
    out['A3_three_diagonals'] = (z(sp.sqrt(1 - sp.Rational(1, 2)) - 1/sp.sqrt(2))                 # held reading at r = 2 r_s (DG1): N
                                 and z(sp.sqrt(cf2.subs(x, sp.Rational(1, 3))) - 1/sp.sqrt(2)))       # free circle at r = 3 r_s
    # stable free circles see more than they lose: count factor^2 > 1/2  <=>  x < 1/3
    out['A3_stable_means_seen_exceeds_lost'] = bool(cf2.subs(x, sp.Rational(1, 4)) > sp.Rational(1, 2)) and bool(cf2.subs(x, sp.Rational(2, 5)) < sp.Rational(1, 2)) and sp.diff(cf2, x) < 0
    # A4: numbers.  The two pure numbers of the synthesis for an electron and a proton
    hbar, c, G, me, mp_, e_, k_e = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.67262192369e-27, 1.602176634e-19, 8.9875517923e9
    a_c = k_e*e_**2/(hbar*c); a_g = G*me*mp_/(hbar*c)
    num['alpha_charge'] = a_c; num['alpha_mass_electron_proton'] = a_g; num['ratio'] = a_c/a_g
    num['mass_of_a_centre_with_no_first_circle_for_an_electron_kg'] = float(1/(2*sp.sqrt(3)))*hbar*c/(G*me)
    out['A4_both_far_below_their_limits'] = a_c < 1 and a_g < 0.2887
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'AC1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    print(run())
