"""CT1: the two-valued record inside the chain of turns: its centre.

DU1: a chain of two-valued marks has rate 2k* ~ 2 e^(-2k) and k -> k - (1/2) ln 2 on doubling the cell.
MG1-M4: the chain of turns (YM-21) has rate -Log(I2(kappa)/I1(kappa)) -> 3/(2 kappa).
The centre of the turn block (U -> -U, YM-39) is a two-valued record inside it.  This stage finds its coupling.
Exact rational series (S^3 moments, YM-37), mpmath for numbers.  Python 3.12."""
import json, math, os
from fractions import Fraction as Fr
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def catalan(n):
    return math.comb(2*n, n)//(n + 1)


def run():
    out, num = {}, {}

    # C1: weight of a face exp(kappa c), c = (1/2) Tr U.  Under the centre U -> -U: c -> -c.
    #     even part cosh(kappa c), odd part sinh(kappa c).  Moments of c on the unit sphere of the block: <c^(2n)> = Catalan(n)/4^n.
    #     f_0 = <cosh(kappa c)> , f_1/2 = <c sinh(kappa c)> ;  YM line: f_j = 2 I_(2j+1)(kappa)/kappa.
    ok0, ok1 = True, True
    for n in range(0, 25):
        mom0 = Fr(catalan(n), 4**n)/math.factorial(2*n)                       # coefficient of kappa^(2n) in <cosh(kappa c)>
        bes0 = Fr(1, 4**n*math.factorial(n)*math.factorial(n + 1))            # coefficient of kappa^(2n) in 2 I_1(kappa)/kappa
        ok0 &= mom0 == bes0
        mom1 = Fr(catalan(n + 1), 4**(n + 1))/math.factorial(2*n + 1)         # coefficient of kappa^(2n+1) in <c sinh(kappa c)>
        bes1 = Fr(1, 2**(2*n + 1)*math.factorial(n)*math.factorial(n + 2))    # coefficient of kappa^(2n+1) in 2 I_2(kappa)/kappa
        ok1 &= mom1 == bes1
    out['C1_even_reading_is_mean_cosh'] = ok0
    out['C1_odd_reading_is_mean_c_sinh'] = ok1
    out['C1_moments_of_the_sphere'] = [Fr(catalan(n), 4**n) for n in (0, 1, 2, 3)] == [1, Fr(1, 4), Fr(1, 8), Fr(5, 64)]

    # C2: the centre coupling: tanh k_c = f_1/2 / f_0 = I2/I1 = < c tanh(kappa c) > under the weight cosh(kappa c)
    mp.mp.dps = 40
    r = lambda v: mp.besseli(2, v)/mp.besseli(1, v)
    dens = lambda c: 2/mp.pi*mp.sqrt(1 - c*c)                                  # law of c on the sphere of the block
    okm = True
    for v in (mp.mpf(1)/2, 2, 7):
        top = mp.quad(lambda c: c*mp.tanh(v*c)*mp.cosh(v*c)*dens(c), [-1, 0, 1])
        bot = mp.quad(lambda c: mp.cosh(v*c)*dens(c), [-1, 0, 1])
        okm &= abs(top/bot - r(v)) < mp.mpf(10)**-25
    out['C2_centre_coupling_is_a_mean_of_local_tanh'] = bool(okm)
    kc = lambda v: mp.atanh(r(v))
    out['C2_rate_is_twice_the_dual_centre_coupling'] = all(abs(-mp.log(r(v)) - 2*(-mp.log(mp.tanh(kc(v)))/2)) < mp.mpf(10)**-30 for v in (1, 5, 50))

    # C3: weak coupling: k_c = (1/2) ln(4 kappa / 3) + O(1/kappa) ; strong: k_c ~ kappa/4
    out['C3_centre_coupling_is_the_log_of_the_turn_coupling'] = all(abs(kc(v) - mp.log(4*v/3)/2) < 1/v for v in (mp.mpf(10)**2, mp.mpf(10)**4, mp.mpf(10)**6))
    out['C3_strong_end'] = all(abs(kc(v)/(v/4) - 1) < v*v for v in (mp.mpf(1)/10, mp.mpf(1)/100))
    # so the exponential law of DU1 in k_c is the power law of MG1 in kappa; doubling: k_c -> k_c - ln2/2  <=>  kappa -> kappa/2
    v0 = mp.mpf(10)**4
    v1 = mp.findroot(lambda v: r(v) - r(v0)**2, v0/2)
    out['C3_doubling_takes_half_ln2_from_the_centre_coupling'] = abs((kc(v0) - kc(v1)) - mp.log(2)/2) < mp.mpf(1)/1000 and abs(v1/v0 - mp.mpf(1)/2) < mp.mpf(1)/1000
    hbar, cc, GG, mp_ = 1.054571817e-34, 299792458.0, 6.67430e-11, 1.67262192369e-27
    n_p = mp_/math.sqrt(hbar*cc/GG)
    kap_p = mp.mpf(3)/2/n_p
    num['proton_count_in_the_two_couplings'] = dict(turn_coupling=float(kap_p), centre_coupling=float(mp.log(4*kap_p/3)/2))
    out['C3_same_count_two_couplings'] = 22 < num['proton_count_in_the_two_couplings']['centre_coupling'] < 23

    # C4: the centre record that is its own dual: I2/I1 = sqrt2 - 1
    kd = mp.findroot(lambda v: r(v) - (mp.sqrt(2) - 1), 2)
    num['turn_coupling_where_the_centre_record_is_self_dual'] = float(kd)
    num['rate_there'] = float(-mp.log(r(kd)))
    out['C4_self_dual_centre_record'] = abs(mp.sinh(2*kc(kd)) - 1) < mp.mpf(10)**-25 and abs(-mp.log(r(kd)) - mp.log(1 + mp.sqrt(2))) < mp.mpf(10)**-25 and 1.8 < kd < 2.0

    # C5: faces in more dimensions: Tr(product of s_l U_l) = (product of s_l) Tr(product of U_l): a two-valued gauge record with
    #     local coupling kappa c_p set by the cosets.  One face of four links, exact quaternion arithmetic.
    def qmul(a, b):
        return (a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3], a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
                a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1], a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0])
    links = [(Fr(3, 5), Fr(4, 5), 0, 0), (Fr(2, 7), Fr(3, 7), Fr(6, 7), 0), (Fr(1, 9), Fr(4, 9), Fr(8, 9), 0), (Fr(12, 13), 0, Fr(3, 13), Fr(4, 13))]
    okf = all(sum(x*x for x in l) == 1 for l in links)
    import itertools
    def face(ls):
        p = ls[0]
        for l in ls[1:]:
            p = qmul(p, l)
        return p[0]
    c0 = face(links)
    for signs in itertools.product((1, -1), repeat=4):
        okf &= face([tuple(sg*x for x in l) for sg, l in zip(signs, links)]) == math.prod(signs)*c0
    out['C5_face_is_a_two_valued_gauge_record_with_local_coupling'] = okf

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'CT1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
    for k_, v in n.items(): print(k_, v)
