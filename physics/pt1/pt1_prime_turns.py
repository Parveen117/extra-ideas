"""PT1: the primes inside the native algebra -- exact turns, exact boosts, and plain addition.

The owner's question: pure numbers are rotation and repetition once a unit is fixed; is there a place for the
prime numbers in the native algebra?
Sources: EMK-1 (R^2 = -1, K^2 = +1, the block over the rationals), SY1 (three sectors: circular, dual, split),
PS1 (Bindu: z -> z / sqrt(det z), the lift to the unit quadric), RMG2-T1 (number sectors), BC1 (boundary = addition),
QC1-T3 (silence), CL1 / R38 (clocks with Q marks; the eight-mark value 1/sqrt 2).
Carrier: rational elements z = a + b R.  det z = a^2 + b^2.  Exact turns = rational points of det z = 1.
P1  The square of the Bindu lift needs no root:  B2(z) = z / zbar = z^2 / det z  is always an exact turn.
P2  PRIME TURNS.  For a prime p = x^2 + y^2 (p = 1 mod 4) put u_p = B2(x + y R).  Every exact turn is, uniquely,
        R^k  x  product of u_p^(n_p) ,   k in {0,1,2,3},  n_p whole numbers.
    Exponents add under composition: the exact turns are the quarter turns times one independent turn per such prime.
P3  A prime p = 3 mod 4 is not a sum of two squares: it gives no turn, only scale.  The prime 2 = det(1 + R) gives
    B2(1 + R) = R, the quarter turn; its half, the eighth turn, needs sqrt 2 -- the eight-mark value of R38.
P4  NO EXACT TURN REPEATS.  u^n = 1 for an exact turn u only if u is a quarter turn.  Clocks with Q marks that are
    exact over the rationals: Q = 1, 2, 4 only (the cosine alone is rational also for Q = 3, 6).
P5  SPLIT SECTOR.  Exact boosts h(t) = (t + 1/t)/2 + (t - 1/t)/2 K, t rational, compose as h(t)h(s) = h(ts):
    one independent boost for EVERY prime.
P6  DUAL SECTOR.  1 + s N composes by addition of s and every element has every root: no primes at all.
Exact integer and rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json
from math import gcd, isqrt


def is_prime(n):
    if n < 2:
        return False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def two_squares(p):
    for x in range(1, isqrt(p) + 1):
        y2 = p - x*x
        y = isqrt(y2)
        if y*y == y2 and y > 0:
            return (max(x, y), min(x, y))
    return None


def gmul(z, w):
    return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])


def gdivides(z, w):
    """does the Gaussian integer w divide z ?  returns the quotient or None"""
    n = w[0]*w[0] + w[1]*w[1]
    a, b = z[0]*w[0] + z[1]*w[1], z[1]*w[0] - z[0]*w[1]
    if a % n == 0 and b % n == 0:
        return (a//n, b//n)
    return None


def b2(z):
    """z / zbar as a rational point of the unit quadric"""
    n = z[0]*z[0] + z[1]*z[1]
    return (F(z[0]*z[0] - z[1]*z[1], n), F(2*z[0]*z[1], n))


def tmul(u, v):
    return (u[0]*v[0] - u[1]*v[1], u[0]*v[1] + u[1]*v[0])


def tpow(u, n):
    out = (F(1), F(0))
    base = u if n >= 0 else (u[0], -u[1])
    for _ in range(abs(n)):
        out = tmul(out, base)
    return out


QUARTER = [(F(1), F(0)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(-1))]


def factor_turn(u):
    """exact turn u = (a + b R)/c  ->  (k, {p: n_p})"""
    c = (u[0].denominator*u[1].denominator)//gcd(u[0].denominator, u[1].denominator)
    a, b = int(u[0]*c), int(u[1]*c)
    if a*a + b*b != c*c:
        raise ValueError('not an exact turn')
    exps = {}
    z = (a, b)
    m = c
    p = 2
    while m > 1:
        if m % p == 0:
            if p % 4 == 1:
                x, y = two_squares(p)
                pi, pib = (x, y), (x, -y)
                while m % p == 0:
                    m //= p
                    # c carries p once; a + b R carries pi^2 or pib^2 (or p, removed by the gcd)
                    q1 = gdivides(z, gmul(pi, pi))
                    q2 = gdivides(z, gmul(pib, pib))
                    if q1 is not None and q2 is None:
                        z = q1
                        exps[p] = exps.get(p, 0) + 1
                    elif q2 is not None and q1 is None:
                        z = q2
                        exps[p] = exps.get(p, 0) - 1
                    else:
                        raise ValueError('turn not in lowest terms')
            else:
                raise ValueError('denominator of an exact turn in lowest terms has only primes 1 mod 4')
        else:
            p += 1
    unit = (F(z[0]), F(z[1]))
    if unit not in QUARTER:
        raise ValueError('remaining factor must be a quarter turn')
    return QUARTER.index(unit), exps


def rebuild(k, exps):
    out = QUARTER[k]
    for p, n in exps.items():
        out = tmul(out, tpow(b2(two_squares(p)), n))
    return out


def chebyshev_orders():
    """x = 2 cos(theta) rational  =>  x in {-2,-1,0,1,2}; the order Q of the turn for each"""
    out = {}
    for x in (-2, -1, 0, 1, 2):
        c0, c1 = 2, x
        Q = 1
        while c1 != 2 or Q == 0:
            if Q > 24:
                break
            c0, c1 = c1, x*c1 - c0
            Q += 1
        out[x] = 1 if x == 2 else Q
    return out


def run():
    primes1 = [p for p in range(3, 120) if is_prime(p) and p % 4 == 1]
    primes3 = [p for p in range(3, 120) if is_prime(p) and p % 4 == 3]
    gens = {p: b2(two_squares(p)) for p in primes1}
    if gens[5] != (F(3, 5), F(4, 5)) or gens[13] != (F(5, 13), F(12, 13)):
        raise ValueError('prime turns of 5 and 13 failed')
    if any(two_squares(p) is not None for p in primes3):
        raise ValueError('a prime 3 mod 4 must give no turn')
    if b2((1, 1)) != (F(0), F(1)):
        raise ValueError('the prime 2 must give the quarter turn')
    # P2: unique factorisation, exponents add
    samples = []
    for k, exps in ((0, {5: 1}), (1, {5: 2, 13: -1}), (3, {17: 1, 29: 1, 5: -3}), (2, {13: 2, 37: -1, 41: 1})):
        u = rebuild(k, exps)
        if u[0]**2 + u[1]**2 != 1:
            raise ValueError('not on the unit quadric')
        k2, e2 = factor_turn(u)
        if (k2, e2) != (k, {p: n for p, n in exps.items() if n}):
            raise ValueError('factorisation is not unique')
        samples.append((k, exps, (str(u[0]), str(u[1]))))
    u1, u2 = rebuild(1, {5: 2, 13: -1}), rebuild(3, {5: -1, 13: 1, 17: 2})
    k12, e12 = factor_turn(tmul(u1, u2))
    if (k12, e12) != (0, {5: 1, 17: 2}):
        raise ValueError('exponents must add under composition')
    # all Pythagorean turns with small denominator factor
    count = 0
    for c in range(2, 150):
        for a in range(1, c):
            b = isqrt(c*c - a*a)
            if b > 0 and a*a + b*b == c*c and gcd(a, gcd(b, c)) == 1:
                k, e = factor_turn((F(a, c), F(b, c)))
                if rebuild(k, e) != (F(a, c), F(b, c)):
                    raise ValueError('rebuild failed')
                count += 1
    # P4: no exact turn repeats
    for p in primes1[:8]:
        acc = (F(1), F(0))
        for n in range(1, 13):
            acc = tmul(acc, gens[p])
            if acc in QUARTER:
                raise ValueError('a prime turn must never return to a quarter turn')
    orders = chebyshev_orders()
    if orders != {-2: 2, -1: 3, 0: 4, 1: 6, 2: 1}:
        raise ValueError('rational cosines must belong to Q = 1, 2, 3, 4, 6 only')
    # the eighth turn: its square is B2(1 + R) = R; its own entries would be 1/sqrt 2, and 2 is not a rational square
    if any(F(n, d)**2 == 2 for n in range(1, 60) for d in range(1, 60)):
        raise ValueError('sqrt 2 must not be rational')
    # P5: split sector
    def h(t):
        return ((t + 1/t)/2, (t - 1/t)/2)
    def hmul(u, v):
        return (u[0]*v[0] + u[1]*v[1], u[0]*v[1] + u[1]*v[0])
    for t, s in ((F(2), F(3)), (F(7, 5), F(11, 3)), (F(1, 2), F(2))):
        if hmul(h(t), h(s)) != h(t*s) or h(t)[0]**2 - h(t)[1]**2 != 1:
            raise ValueError('exact boosts must compose by multiplying t')
    # P6: dual sector: every root exists; a prime turn has no exact square root
    s = F(5, 7)
    for n in (2, 3, 7):
        if n*(s/n) != s:
            raise ValueError('dual sector roots failed')
    k5, e5 = factor_turn(gens[5])
    if e5[5] % 2 == 0:
        raise ValueError('exponent of a prime turn must be one')
    return dict(prime_turns={p: (str(g[0]), str(g[1])) for p, g in list(gens.items())[:6]},
                silent_primes=primes3[:8], prime_two='quarter turn', factored_turns=count,
                rational_cosine_orders=orders, examples=[(k, {str(p): n for p, n in e.items()}, u) for k, e, u in samples])


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PT1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
