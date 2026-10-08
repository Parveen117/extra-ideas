"""PM1: the prime-turn series.

Source objects: the catalog draft's "prime lattice" {(log p, 2 pi/p)}, "prime-modular form"
sum_p e^(2 pi i tau/p) and "seam zeta function" sum_p p^(-s) e^(-2 pi i/p).  None of the three is used:
the first series has terms that do not tend to zero, and the twist e^(-2 pi i/p) is neither an exact
turn of the native algebra nor multiplicative (P6).  This stage builds the series that the native
algebra does carry: the same idea -- primes entering through turns -- with the turns of PT1.

Sources read before building: physics PT1 (exact turns; one prime turn u_p = pi/pi-bar for each prime
p = 1 mod 4; 2 is the quarter turn; primes 3 mod 4 give no turn; no exact turn repeats), RKF F00-E
(cut-complex field, iota^2 = -1), Publications LAM-1 T3 (theta seam at a point by directed rational
enclosures), LAM-2 (seam phase derived from finite arithmetic; the character mod 4), LAM-3 (the crossing).

Carrier: whole cut-complex numbers alpha = a + b iota, norm N = a^2 + b^2, four units iota^j.

P1  The power alpha^m is the same for the four unit multiples of alpha exactly when 4 | m.
    For other m every norm shell sums to zero.  So among the powers the unit-free ones are alpha^(4k), and
        A_k(n) = (1/4) sum_(N(alpha) = n) alpha^(4k)        is a whole number (no turn-part).
P2  A_k is multiplicative.  On primes:  A_k(2) = (-4)^k ;  A_k(p) = 0 for p = 3 mod 4 ;
        A_k(p) = 2 p^(2k) rad(u_p^(2k))  for p = 1 mod 4   (u_p the prime turn of PT1) ;
        A_k(p^(e+1)) = A_k(p) A_k(p^e) - chi(p) p^(4k) A_k(p^(e-1)) ,  chi the character mod 4.
P3  |A_k(p)| < 2 p^(2k): the size bound is "the rad-part of a turn is below 1"; never reached (PT1-P4).
P4  Product form: the numbers built from the prime rule P2 alone equal the shell sums for every n
    (two routes).  Formally  sum A_k(n) n^(-s) = prod_p 1/(1 - A_k(p) p^(-s) + chi(p) p^(4k-2s)).
P5  Seam law at a point.  Theta_k(t) = sum_alpha alpha^(4k) Exp(-pi N(alpha) t):
        Theta_k(1/t) = t^(4k+1) Theta_k(t) ,   phase +1 = iota^(-4k) ,
    the two sides agree to more than forty digits in two-sided rational enclosures at t = 2 and t = 3/2 for
    k = 1, 2 (and k = 0, the square of LAM-1's theta); equality itself is the pinned classical anchor.
    Wrong weight and wrong sign are separated, and so is a^4 + b^4, which is unit-free but not a power.
P6  Why the draft's two series are not used (exact reasons).

Exact integer arithmetic; directed rational enclosures for P5.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
from math import isqrt
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = 1500
BITS = 480
SCALE = 1 << BITS


# ---------------------------------------------------------------- whole cut-complex numbers
def gmul(z, w):
    return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])


def gpow(z, n):
    out = (1, 0)
    for _ in range(n):
        out = gmul(out, z)
    return out


def shells(nmax):
    out = {n: [] for n in range(nmax + 1)}
    r = isqrt(nmax)
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = a*a + b*b
            if n <= nmax:
                out[n].append((a, b))
    return out


def shell_sum(shell, m):
    r = t = 0
    for z in shell:
        p = gpow(z, m)
        r += p[0]
        t += p[1]
    return r, t


def is_prime(n):
    if n < 2:
        return False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def two_squares(p):
    for x in range(1, isqrt(p) + 1):
        y = isqrt(p - x*x)
        if y > 0 and x*x + y*y == p:
            return (max(x, y), min(x, y))
    return None


def chi4(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


def factor(n):
    out, d = {}, 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def prime_value(p, k):
    """A_k(p) from the prime rule alone (no shell)"""
    if p == 2:
        return (-4)**k
    if p % 4 == 3:
        return 0
    return 2*gpow(two_squares(p), 4*k)[0]


def from_primes(n, k):
    out = 1
    for p, e in factor(n).items():
        a1 = prime_value(p, k)
        prev, cur = 1, a1
        for _ in range(e - 1):
            prev, cur = cur, a1*cur - chi4(p)*p**(4*k)*prev
        out *= cur if e >= 1 else 1
    return out


def prime_turn_power(p, n):
    """u_p^n as an exact rational point of the unit quadric (PT1: u_p = pi/pi-bar = pi^2/p)"""
    z = gpow(two_squares(p), 2*n)
    return F(z[0], p**n), F(z[1], p**n)


# ---------------------------------------------------------------- directed rational enclosures
def down(x):
    return F((x.numerator*SCALE)//x.denominator, SCALE)


def up(x):
    return F(-((-x.numerator*SCALE)//x.denominator), SCALE)


def arctan_inv(q, terms):
    s, parts = F(0), []
    for n in range(terms):
        s += F((-1)**n, (2*n + 1)*q**(2*n + 1))
        parts.append(s)
    return min(parts[-2:]), max(parts[-2:])


def pi_interval():
    a_lo, a_hi = arctan_inv(5, 110)
    b_lo, b_hi = arctan_inv(239, 34)
    return down(16*a_lo - 4*b_hi), up(16*a_hi - 4*b_lo)


def _series_bounds(ri, terms=90):
    """integers lo, hi with  lo <= SCALE*Exp(-ri/SCALE) <= hi ,  0 <= ri <= SCALE/2  (alternating, decreasing terms:
    a partial sum ending on an odd index is below the value, one ending on an even index is above)"""
    lo = hi = 0
    tlo = thi = SCALE
    lo_odd = hi_even = None
    for n in range(terms):
        if n % 2 == 0:
            lo, hi = lo + tlo, hi + thi
            hi_even = hi
        else:
            lo, hi = lo - thi, hi - tlo
            lo_odd = lo
        tlo = (tlo*ri)//((n + 1)*SCALE)
        thi = -((-thi*ri)//((n + 1)*SCALE))
    return lo_odd, hi_even


def exp_neg(q):
    """two-sided enclosure of Exp(-q), q >= 0 rational: factorial series on q/2^s <= 1/2 in directed
    fixed-point arithmetic, then s squarings"""
    s = 0
    while q > F(1, 2)*(1 << s):
        s += 1
    r = q/(1 << s)
    r_dn = (r.numerator*SCALE)//r.denominator
    r_up = -((-r.numerator*SCALE)//r.denominator)
    lo = _series_bounds(r_up)[0]
    hi = _series_bounds(r_dn)[1]
    for _ in range(s):
        lo, hi = (lo*lo)//SCALE, -((-hi*hi)//SCALE)
    return F(lo, SCALE), F(hi, SCALE)


def exp_neg_pi(c, pi_iv):
    """enclosure of Exp(-pi c), c > 0 rational"""
    return exp_neg(pi_iv[1]*c)[0], exp_neg(pi_iv[0]*c)[1]


def weighted_sum(coefs, e_iv, t_power, tail_degree):
    """enclosure of  t_power * ( coefs[0] + sum_(n>=1) coefs[n] E^n ),  E in e_iv, with the tail bounded by
    6 n^tail_degree E^n beyond the last coefficient.  (A shell of norm n has at most 4 sqrt(n) + 2 <= 6n
    numbers, each reading of degree d has size at most n^(d/2); run() checks the count.)"""
    lo = hi = F(coefs[0])
    plo, phi = F(1), F(1)
    for n in range(1, len(coefs)):
        plo, phi = down(plo*e_iv[0]), up(phi*e_iv[1])
        c = coefs[n]
        if c >= 0:
            lo, hi = lo + c*plo, hi + c*phi
        else:
            lo, hi = lo + c*phi, hi + c*plo
    m = len(coefs)                                         # first omitted index
    rho = F(m + 1, m)**tail_degree*e_iv[1]
    assert rho < 1
    tail = 6*m**tail_degree*up(phi*e_iv[1])/(1 - rho)
    lo, hi = lo - tail, hi + tail
    if t_power < 0:                                         # a negative factor exchanges the two ends
        lo, hi = hi, lo
    return down(t_power*lo), up(t_power*hi)


def overlap(x, y):
    assert x[0] <= x[1] and y[0] <= y[1], 'not an enclosure'
    return max(x[0], y[0]) <= min(x[1], y[1])


def gap(x, y):
    return max(x[0], y[0]) - min(x[1], y[1])


def run():
    out, num = {}, {}
    sh = shells(NMAX)
    ks = (1, 2)

    # P1 -------------------------------------------------------------------------------------------------
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    z = (3, 2)
    out['P1_power_is_unit_free_exactly_for_multiples_of_four'] = all(
        (len({gpow(gmul(u, z), m) for u in units}) == 1) == (m % 4 == 0) for m in range(1, 17))
    out['P1_other_readings_sum_to_zero_on_every_shell'] = all(
        shell_sum(sh[n], m) == (0, 0) for n in range(1, 301) for m in (1, 2, 3, 5, 6, 7))
    out['P1_shell_sums_are_whole_numbers_without_turn_part'] = all(
        shell_sum(sh[n], 4*k)[1] == 0 and shell_sum(sh[n], 4*k)[0] % 4 == 0 for n in range(1, NMAX + 1) for k in ks)
    A = {k: [shell_sum(sh[n], 4*k)[0]//4 for n in range(NMAX + 1)] for k in ks}
    for k in ks:
        A[k][0] = 0
    num['A_1_first_values'] = A[1][1:14]

    # P2 -------------------------------------------------------------------------------------------------
    from math import gcd
    out['P2_multiplicative'] = all(
        A[k][m*n] == A[k][m]*A[k][n]
        for k in ks for m in range(1, 39) for n in range(1, NMAX//m + 1) if gcd(m, n) == 1)
    primes = [p for p in range(2, NMAX + 1) if is_prime(p)]
    out['P2_two_is_the_quarter_turn'] = all(A[k][2] == (-4)**k for k in ks) and gmul((1, 1), (1, 1)) == (0, 2)
    out['P2_primes_three_mod_four_are_silent'] = all(A[k][p] == 0 for k in ks for p in primes if p % 4 == 3)
    ok = True
    for k in ks:
        for p in primes:
            if p % 4 == 1:
                rad, turn = prime_turn_power(p, 2*k)
                ok &= rad*rad + turn*turn == 1 and A[k][p] == 2*p**(2*k)*rad
    out['P2_split_primes_read_the_prime_turn'] = ok
    ok = True
    for k in ks:
        for p in primes:
            e = 1
            while p**(e + 1) <= NMAX:
                ok &= A[k][p**(e + 1)] == A[k][p]*A[k][p**e] - chi4(p)*p**(4*k)*A[k][p**(e - 1) if e > 1 else 1]
                e += 1
    out['P2_prime_power_rule'] = ok
    num['prime_turns'] = {str(p): [str(x) for x in prime_turn_power(p, 1)] for p in (5, 13, 17, 29)}

    # P3 -------------------------------------------------------------------------------------------------
    out['P3_size_bound_strict'] = all(
        abs(A[k][p]) < 2*p**(2*k) for k in ks for p in primes if p % 4 == 1)

    # P4 -------------------------------------------------------------------------------------------------
    out['P4_product_form_equals_shell_sums'] = all(
        from_primes(n, k) == A[k][n] for k in ks for n in range(1, NMAX + 1))

    # P5 -------------------------------------------------------------------------------------------------
    pi_iv = pi_interval()
    out['P5_pi_enclosure'] = pi_iv[0] < F(314159265358979323847, 10**20) and \
        F(314159265358979323846, 10**20) < pi_iv[1] and pi_iv[1] - pi_iv[0] < F(1, 10**70)
    cut = 110
    theta_coefs = {0: [1] + [len(sh[n]) for n in range(1, cut)]}
    for k in ks:
        theta_coefs[k] = [0] + [4*A[k][n] for n in range(1, cut)]
    quartic = [0] + [sum(a**4 + b**4 for a, b in sh[n]) for n in range(1, cut)]     # not a unit-free reading

    def side(coefs, c, weight, deg):
        return weighted_sum(coefs, exp_neg_pi(c, pi_iv), weight, deg)

    ok_seam = ok_weight = ok_sign = True
    widths, gaps = [], []
    for k in (0, 1, 2):
        deg = 2*k + 1
        for t in (F(2), F(3, 2)):
            left = side(theta_coefs[k], 1/t, F(1), deg)
            right = side(theta_coefs[k], t, t**(4*k + 1), deg)
            ok_seam &= overlap(left, right)
            widths.append(max(left[1] - left[0], right[1] - right[0]))
            for wrong in (4*k, 4*k + 2):
                bad = side(theta_coefs[k], t, t**wrong, deg)
                ok_weight &= not overlap(left, bad)
                gaps.append(gap(left, bad))
            neg = side(theta_coefs[k], t, -t**(4*k + 1), deg)
            ok_sign &= not overlap(left, neg)
            if k == 1:
                num[f'Theta_1(1/t), t={t}'] = str(round(float(left[0]), 15))
    out['P5_seam_law_at_the_points'] = ok_seam
    out['P5_enclosures_agree_to_forty_digits'] = max(widths) < F(1, 10**40)
    out['P5_wrong_weight_is_separated'] = ok_weight and min(gaps) > F(1, 10**3)
    out['P5_wrong_sign_is_separated'] = ok_sign
    out['P5_shell_count_bound_used_in_the_tails'] = all(len(sh[n]) <= 6*n for n in range(1, NMAX + 1))
    # a^4 + b^4 = (3/4) N^2 + (1/4) rad(alpha^4): unit-free, but not a power
    probes = [(3, 2), (1, 5), (7, -4)]
    out['P5_a4_plus_b4_is_unit_free_but_not_a_power'] = all(
        len({gmul(u, zz)[0]**4 + gmul(u, zz)[1]**4 for u in units}) == 1
        and 4*(zz[0]**4 + zz[1]**4) == 3*(zz[0]**2 + zz[1]**2)**2 + gpow(zz, 4)[0] for zz in probes)
    ql = side(quartic, F(1, 2), F(1), 3)
    qr = side(quartic, F(2), F(2)**5, 3)
    out['P5_unit_free_reading_that_is_not_a_power_does_not_cross'] = not overlap(ql, qr) and gap(ql, qr) > F(1, 100)
    num['P5_gap_for_a4_plus_b4'] = str(round(float(gap(ql, qr)), 6))

    # P6 -------------------------------------------------------------------------------------------------
    # (a) the terms of sum_p e^(2 pi i tau/p) have modulus Exp(-2 pi y/p), which stays above 1/2 for p > 14 y
    ok = True
    for y in (F(1, 10), F(1), F(5)):
        for p in [q for q in primes if q > 14*y][:40]:
            ok &= exp_neg(pi_iv[1]*2*y/p)[0] > F(1, 2)
    out['P6_draft_prime_series_terms_do_not_tend_to_zero'] = ok
    # (b) the twist e^(-2 pi i/p) has finite order p, while no exact turn other than a quarter turn repeats
    #     (PT1-P4): checked here for the prime turns
    def turn_order(u, limit=24):
        cur = u
        for n in range(1, limit + 1):
            if cur == (F(1), F(0)):
                return n
            cur = (cur[0]*u[0] - cur[1]*u[1], cur[0]*u[1] + cur[1]*u[0])
        return None
    out['P6_prime_turns_never_repeat'] = all(turn_order(prime_turn_power(p, 1)) is None for p in primes[:60] if p % 4 == 1)

    # controls ----------------------------------------------------------------------------------------------
    bad = list(A[1])
    bad[65] += 1
    out['control_planted_coefficient_breaks_the_product_form'] = any(from_primes(n, 1) != bad[n] for n in range(1, 100))
    out['control_rule_with_the_wrong_character_fails'] = any(
        A[1][p*p] != A[1][p]*A[1][p] + chi4(p)*p**4 for p in (3, 5, 7))
    out['control_third_power_is_not_unit_free'] = len({gpow(gmul(u, z), 3) for u in units}) == 4

    out['pass'] = all(v for v in out.values() if isinstance(v, bool))
    return out, num


if __name__ == '__main__':
    out, num = run()
    with open(os.path.join(HERE, 'PM1_RESULT.json'), 'w') as fh:
        json.dump({'checks': out, 'numbers': num}, fh, indent=1, sort_keys=True)
        fh.write('\n')
    for k, v in out.items():
        print(f'{k}: {v}')
    print(num)
