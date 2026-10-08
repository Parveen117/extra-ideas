"""RD1: the pure numbers of a gas of light from whole counts.

WQ1: a recordable mode has E = n kappa omega.  QC2: the ensemble weighs a reading by exp(-energy / unit); theta = the thermal unit.
ZP1: modes between walls.  EM1: light has two readings per direction (d - 1 = 2 in three cuts).
Exact rational arithmetic, sympy and mpmath.  Explicit tail bounds in the sense of theorum/28 section 4.  Python 3.12."""
import json, os
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 30


def run():
    out, num = {}, {}
    # L1: whole counts with weights x^n, x = exp(-kappa omega / theta): mean count x/(1-x), exactly, with an explicit tail
    ok = True
    for x in (Fr(1, 2), Fr(1, 3), Fr(3, 5)):
        Kmax = 60
        Zs = sum(x**n for n in range(Kmax)); Ns = sum(n*x**n for n in range(Kmax))
        tailZ = x**Kmax/(1 - x)
        ok &= (Zs + tailZ == 1/(1 - x))
        tailN = x**Kmax*(Kmax/(1 - x) + x/(1 - x)**2)
        ok &= (Ns + tailN == x/(1 - x)**2)
        ok &= ((Ns + tailN)/(Zs + tailZ) == x/(1 - x))
    out['L1_mean_count'] = ok
    # L2: integral of y^3/(e^y - 1) = sum_k 6/k^4, with tail below 2/K^3 for every K : value pi^4/15
    K = 200
    part = sum(Fr(6, k**4) for k in range(1, K + 1))
    target = mp.pi**4/15
    out['L2_energy_integral'] = (mp.mpf(part.numerator)/part.denominator < target < mp.mpf(part.numerator)/part.denominator + mp.mpf(2)/K**3)
    part3 = sum(Fr(2, k**3) for k in range(1, K + 1))                # integral of y^2/(e^y - 1) = 2 zeta(3), tail below 1/K^2
    z3 = mp.zeta(3)
    out['L2_count_integral'] = (mp.mpf(part3.numerator)/part3.denominator < 2*z3 < mp.mpf(part3.numerator)/part3.denominator + mp.mpf(1)/K**2)
    # L3: the pure numbers, three cuts, two readings per direction
    num['energy_density_number_pi2_over_15'] = float(mp.pi**2/15)     # u (kappa c)^3 / theta^4
    num['flux_number_pi2_over_60'] = float(mp.pi**2/60)
    num['count_density_number_2zeta3_over_pi2'] = float(2*z3/mp.pi**2)
    num['energy_per_count_in_theta'] = float(mp.pi**4/(30*z3))
    num['entropy_per_count'] = float(2*mp.pi**4/(45*z3))
    out['L3_numbers'] = abs(num['energy_per_count_in_theta'] - 2.70118) < 1e-5 and abs(num['entropy_per_count'] - 3.60157) < 1e-5
    # L4: where the spectrum peaks: y = 3(1 - e^-y) per unit rate ; y = 5(1 - e^-y) per unit length
    y3 = mp.findroot(lambda y: y - 3*(1 - mp.e**(-y)), 2.8)
    y5 = mp.findroot(lambda y: y - 5*(1 - mp.e**(-y)), 4.9)
    num['peak_per_rate'] = float(y3); num['peak_per_length'] = float(y5)
    out['L4_peaks'] = abs(y3 - mp.mpf('2.821439372')) < 1e-8 and abs(y5 - mp.mpf('4.965114232')) < 1e-8
    # L5: from TD1: pressure/energy density = 1/d , entropy x theta / energy = (d+1)/d , memory m = 1
    d = sp.Symbol('d', positive=True)
    S, V = sp.symbols('S V', positive=True)
    U = S**((d + 1)/d)*V**(-1/d)
    out['L5_ratios'] = sp.simplify(-sp.diff(U, V)*V/U - 1/d) == 0 and sp.simplify(sp.diff(U, S)*S/U - (d + 1)/d) == 0
    # d cuts: the energy exponent is d + 1 (four in three cuts): theta ~ U_S ~ (S/V)^(1/d)  =>  U/V ~ theta^(d+1)
    th = sp.simplify(sp.diff(U, S))
    out['L5_fourth_power_in_three_cuts'] = all(sp.simplify(((U/V)/th**(d + 1)).subs(d, n) - sp.Rational(n, n + 1)**(n + 1)) == 0 for n in (1, 2, 3, 4))
    # measured (arXiv:astro-ph/9605054, abstract, read at source): the sky's light gas
    num['measured'] = dict(deviation_from_the_shape_ppm_of_peak='< 50', count_potential_over_theta='< 9e-5 (95%)', temperature_K='2.728 +- 0.004')
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'RD1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    o, n = run(); print(o); print(n)
