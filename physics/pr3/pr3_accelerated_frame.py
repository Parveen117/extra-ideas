"""PR3: the sector law in a uniformly accelerated frame.

Inertial law (c = 1):  (d_t - A d_x - g R) psi = 0,  A = H = diag(1, -1), R = K H.
Null coordinates u = x + t, v = x - t;  frame coordinates rho^2 = u v, e^{2 eta} = u / v.
Functions are finite sums of monomials u^a v^b with rational exponents; every identity
below is checked on monomials, hence on all such sums.  Exact, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json

HALF = F(1, 2)


# ------------------------------------------------------- monomial calculus
def clean(f):
    return {k: c for k, c in f.items() if c}


def fadd(f, g, s=1):
    out = dict(f)
    for k, c in g.items():
        out[k] = out.get(k, F(0))+s*c
    return clean(out)


def scale(c, f):
    return clean({k: c*v for k, v in f.items()})


def mul(f, a, b, c=F(1)):
    """multiply by c u^a v^b"""
    return clean({(k[0]+a, k[1]+b): c*v for k, v in f.items()})


def du(f):
    return clean({(k[0]-1, k[1]): v*k[0] for k, v in f.items()})


def dv(f):
    return clean({(k[0], k[1]-1): v*k[1] for k, v in f.items()})


def d_eta(f):                       # u d_u - v d_v
    return clean({k: v*(k[0]-k[1]) for k, v in f.items()})


def d_s(f):                         # rho d_rho = u d_u + v d_v
    return clean({k: v*(k[0]+k[1]) for k, v in f.items()})


def rho(f, power=1, c=F(1)):
    return mul(f, F(power, 2), F(power, 2), c)


# ------------------------------------------------------------- the two laws
def inertial(psi, g):
    """(d_t - A d_x - g R) psi with d_t - d_x = -2 d_v, d_t + d_x = 2 d_u, R psi = (-psi2, psi1)."""
    p1, p2 = psi
    return (fadd(scale(F(-2), dv(p1)), scale(g, p2)), fadd(scale(F(2), du(p2)), scale(g, p1), -1))


def frame_law(chi, g, split=F(0)):
    """(1/rho) d_eta chi - A (d_rho + split/rho) chi - g R chi, d_rho = (1/rho) d_s."""
    c1, c2 = chi
    a1 = fadd(d_s(c1), scale(split, c1))
    a2 = fadd(d_s(c2), scale(split, c2))
    r1 = fadd(rho(fadd(d_eta(c1), a1, -1), -1), scale(g, c2))
    r2 = fadd(rho(fadd(d_eta(c2), a2), -1), scale(g, c1), -1)
    return r1, r2


TESTS = [({(F(1), F(2)): F(1)}, {(F(0), F(3)): F(2)}),
         ({(F(3, 2), F(-1)): F(-1), (F(0), F(0)): F(5)}, {(F(2), F(1, 2)): F(1, 3)}),
         ({(F(-1, 2), F(5, 2)): F(7)}, {(F(1), F(1)): F(-2), (F(4), F(-3, 2)): F(1)})]


def boost_only_control(g):
    """psi = Exp(-eta A/2) phi: the frame law carries the split term A/(2 rho)."""
    for phi in TESTS:
        psi = (mul(phi[0], F(-1, 4), F(1, 4)), mul(phi[1], F(1, 4), F(-1, 4)))     # E^-1/2, E^1/2
        left = inertial(psi, g)
        law = frame_law(phi, g, split=HALF)
        right = (mul(law[0], F(1, 4), F(-1, 4)), mul(law[1], F(-1, 4), F(1, 4)))   # Exp(eta A/2)
        if left != right:
            raise ValueError('boosted law with split term A/(2 rho) failed')
        bare = frame_law(phi, g, split=F(0))
        if (mul(bare[0], F(1, 4), F(-1, 4)), mul(bare[1], F(-1, 4), F(1, 4))) == left:
            raise ValueError('split term should be needed without the half-density')
    return dict(split_term='A / (2 rho)', needed=True)


def half_density_control(g):
    """psi_1 = u^-1/2 chi_1, psi_2 = v^-1/2 chi_2: the frame law is the inertial one with d_t -> (1/rho) d_eta."""
    for chi in TESTS:
        psi = (mul(chi[0], -HALF, F(0)), mul(chi[1], F(0), -HALF))
        left = inertial(psi, g)
        law = frame_law(chi, g, split=F(0))
        right = (mul(law[0], F(0), -HALF), mul(law[1], -HALF, F(0)))
        if left != right:
            raise ValueError('half-density law failed')
        withsplit = frame_law(chi, g, split=HALF)
        if (mul(withsplit[0], F(0), -HALF), mul(withsplit[1], -HALF, F(0))) == left:
            raise ValueError('split term should be absent in half-density variables')
    return dict(split_term='0', law='(1/rho) d_eta chi = (A d_rho + g R) chi')


# ---------------------------------------- generator G = A d_s + g rho R and G^2
def G(chi, g):
    c1, c2 = chi
    return (fadd(d_s(c1), rho(c2, 1, g), -1), fadd(scale(F(-1), d_s(c2)), rho(c1, 1, g)))


def square_control(g):
    for chi in TESTS:
        c1, c2 = chi
        left = G(G(chi, g), g)
        # d_s^2 - g^2 rho^2 - g rho K,  K chi = (chi2, chi1)
        right = (fadd(fadd(d_s(d_s(c1)), rho(c1, 2, g*g), -1), rho(c2, 1, g), -1),
                 fadd(fadd(d_s(d_s(c2)), rho(c2, 2, g*g), -1), rho(c1, 1, g), -1))
        if left != right:
            raise ValueError('G^2 = d_s^2 - g^2 rho^2 - g rho K failed')
        # in the K eigen-combinations c = chi1 +- chi2 the operator -G^2 factorizes
        for sign in (1, -1):
            c = fadd(c1, c2, sign)
            minus_G2 = fadd(fadd(scale(F(-1), d_s(d_s(c))), rho(c, 2, g*g)), rho(c, 1, sign*g))
            if sign == 1:      # (d_s + W)(-d_s + W) c,  W = g rho
                inner = fadd(scale(F(-1), d_s(c)), rho(c, 1, g))
                fact = fadd(d_s(inner), rho(inner, 1, g))
            else:              # (-d_s + W)(d_s + W) c
                inner = fadd(d_s(c), rho(c, 1, g))
                fact = fadd(scale(F(-1), d_s(inner)), rho(inner, 1, g))
            if minus_G2 != fact:
                raise ValueError('factorization of -G^2 failed')
    # the well of the minus branch: V = g rho (g rho - 1), minimum -1/4 at g rho = 1/2
    well = [(x, x*(x-1)) for x in (F(1, 4), F(1, 2), F(3, 4), F(1), F(2))]
    if min(v for _, v in well) != F(-1, 4) or dict(well)[HALF] != F(-1, 4):
        raise ValueError('well depth is not 1/4')
    return dict(square='d_s^2 - g^2 rho^2 - g rho K', well_depth='1/4 at g rho = 1/2',
                factorized='(-d_s + W)(d_s + W) and (d_s + W)(-d_s + W), W = g rho')


def antimass_control(g):
    """A G_g A = G_{-g}: the reversed turn obeys the same frame law with the same clock factor."""
    for chi in TESTS:
        flipped = (chi[0], scale(F(-1), chi[1]))                 # A chi
        out = G(flipped, g)
        left = (out[0], scale(F(-1), out[1]))                    # A G_g A chi
        if left != G(chi, -g):
            raise ValueError('mass and antimass do not share the frame law')
    return dict(relation='A G_g A = G_{-g}', clock_factor='rho for both sheets')


def lapse_control():
    """rho^{n/2} d_rho rho^{n/2} = rho^n d_rho + (n/2) rho^{n-1} on powers of rho."""
    for n in (F(1), F(2), F(3, 2), F(-1)):
        for k in (F(0), F(1), F(5, 2), F(-2)):
            f = {(k/2, k/2): F(1)}                               # rho^k
            d_rho = lambda h: rho(d_s(h), -1)
            left = mul(d_rho(mul(f, n/4, n/4)), n/4, n/4)
            right = fadd(mul(d_rho(f), n/2, n/2), mul(f, (n-1)/2, (n-1)/2, n/2))
            if left != right:
                raise ValueError('lapse lemma failed')
    return dict(lemma='N^{1/2} d N^{1/2} = N d + N\'/2', checked='N = rho^n')


def run(g=F(3, 2)):
    return dict(boost_only=boost_only_control(g), half_density=half_density_control(g),
                square=square_control(g), antimass=antimass_control(g), lapse=lapse_control(),
                redshift=dict(rest_turn_per_frame_time='g rho', coordinate_speed='rho',
                              ratio_between_two_heights='rho_1 / rho_2'))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('PR3_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
