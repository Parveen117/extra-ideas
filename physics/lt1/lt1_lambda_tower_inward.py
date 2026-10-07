"""LT1: the lambda-tower on the gravity field -- moving inward, and what sits at the end.

The owner's statement: raising lambda in the tower moves inward in the diagram; any lambda can be taken as the
reference; at infinite lambda the secret of counting may hide.
Sources: RMG2 (lambda-tower T_lambda(H) = H + lambda H^2; acts radially; l' = l + log((1+lambda l+)/(1+lambda l-)),
between l and 2l; fixed sets: the centre and the nilpotent boundary; T1: number sectors), NC1 (response element
G = Exp(psi n.C), psi = 2 eta, tanh^2 eta = m = r_s/r), GR1 (N = sech eta), IN1 (invariant n^2 - r.r), SY1
(dual sector: powers count exactly).
Write rho = tanh(psi) = 2 beta / (1 + beta^2) for the odd ratio of G (RMG's rho), q = e^psi.
T1  One generation at infinite lambda is squaring: rapidity doubles, m' = 4m/(1+m)^2, and the place moves
        r' = (r + r_s)^2 / (4 r) :   inward for r > r_s, fixed at r = r_s and at infinity.
T2  Finite lambda interpolates: eigenvalue ratio q^2 -> q^3 (1 + lambda q)/(q + lambda), increasing in lambda from
    q^2 (lambda = 0) to q^4 (lambda = infinity).  Changing lambda moves continuously inward.
T3  Two references.  With r = 4 s the map of T1 is  r' = s (1 + r_s/(4 s))^2 : the relation between the radius in
    which space is conformally flat (s) and the radius of MO1 in which space is flat (r').  One generation of the
    tower carries the law "memory = r_s / radius" from one reference to the other; N = (1 - x)/(1 + x), x = r_s/4s.
T4  The invariant of the normalised element is 1 - rho^2; a generation sends it to ((1 - rho^2)/(1 + rho^2))^2.
    Iterating, it goes to zero: the element becomes rank one -- a single null reading (IN1) -- at the horizon.
T5  At that boundary the quarter-turned element squares to zero (RMG2-T1c), and there powers count exactly:
    (1 + t N)^k = 1 + k t N.
Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json


def square_m(m):
    return 4*m/(1+m)**2


def radius_map(r, rs):
    return (r+rs)**2/(4*r)


def ratio_after(q, lam):
    """eigenvalue ratio of T_lambda(G), G with eigenvalues q, 1/q."""
    return (q*(1+lam*q))/((1/q)*(1+lam/q))


def run():
    rs = F(1)
    # T1
    rows = []
    for r in (F(2), F(5), F(40), F(3, 2)):
        m = rs/r
        m2 = square_m(m)
        r2 = radius_map(r, rs)
        if rs/r2 != m2 or not (rs < r2 < r):
            raise ValueError('generation must move inward and stay outside the horizon')
        rows.append((str(r), str(r2)))
    if radius_map(rs, rs) != rs:
        raise ValueError('horizon must be fixed')
    # doubling law: tanh(2 eta) = 2 t/(1+t^2) with t^2 = m
    for t in (F(1, 3), F(2, 5), F(7, 9)):
        if (2*t/(1+t*t))**2 != square_m(t*t):
            raise ValueError('rapidity doubling failed')
    orbit = [F(40)]
    for _ in range(6):
        orbit.append(radius_map(orbit[-1], rs))
    if any(not (rs < b < a) for a, b in zip(orbit, orbit[1:])):
        raise ValueError('iteration must decrease to the horizon')
    # T2
    q = F(3, 2)
    prev = None
    lam_rows = []
    for lam in (F(0), F(1, 10), F(1), F(10), F(1000)):
        R = ratio_after(q, lam)
        if not (q**2 <= R < q**4) or (prev is not None and R <= prev):
            raise ValueError('finite lambda must interpolate monotonically')
        if R != q**3*(1+lam*q)/(q+lam):
            raise ValueError('ratio formula failed')
        prev = R
        lam_rows.append((str(lam), str(R)))
    if ratio_after(q, F(0)) != q**2 or q**4 - ratio_after(q, F(10**9)) > F(1, 10**6):
        raise ValueError('end points failed')
    # T3
    iso = []
    for s in (F(1), F(3, 4), F(5), F(1, 3)):
        x = rs/(4*s)
        areal = s*(1+x)**2
        if radius_map(4*s, rs) != areal:
            raise ValueError('the generation is not the change of reference radius')
        N = (1-x)/(1+x)
        if x < 1 and N*N != 1-rs/areal:
            raise ValueError('clock factor in the conformally flat reference failed')
        # the lower level: tanh^2(eta/2) = x = r_s / (4 s); one generation doubles it to tanh^2(eta) = r_s / areal
        if square_m(x) != rs/areal:
            raise ValueError('doubling from the lower level failed')
        iso.append((str(s), str(areal), str(N)))
    # T4
    inv = []
    rho = F(3, 5)
    for _ in range(5):
        new_rho = 2*rho/(1+rho*rho)
        if 1-new_rho**2 != ((1-rho**2)/(1+rho**2))**2:
            raise ValueError('invariant law of a generation failed')
        inv.append(str(1-rho**2))
        rho = new_rho
    if not all(F(a) > F(b) for a, b in zip(inv, inv[1:])):
        raise ValueError('invariant must decrease')
    # rank one at rho = 1: (1 + W)/2 with W^2 = 1 is idempotent
    W = [[F(3, 5), F(4, 5)], [F(4, 5), F(-3, 5)]]
    P = [[(int(i == j)+W[i][j])/2 for j in range(2)] for i in range(2)]
    PP = [[sum(P[i][k]*P[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    det = P[0][0]*P[1][1]-P[0][1]*P[1][0]
    if PP != P or det != 0:
        raise ValueError('boundary element must be rank one')
    # T5: dual numbers count exactly
    Nn = [[F(0), F(1)], [F(0), F(0)]]
    def mul(a, b):
        return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    t = F(2, 7)
    U = [[F(1), t], [F(0), F(1)]]
    acc = [[F(1), F(0)], [F(0), F(1)]]
    for k in range(1, 8):
        acc = mul(acc, U)
        if acc != [[F(1), k*t], [F(0), F(1)]]:
            raise ValueError('powers must count exactly in the dual sector')
    if mul(Nn, Nn) != [[F(0), F(0)], [F(0), F(0)]]:
        raise ValueError('nilpotent failed')
    return dict(inward=rows, orbit=[str(x) for x in orbit[:5]], lambda_ratios=lam_rows, two_references=iso,
                invariant_per_generation=inv, far_field_factor=str(F(1, 4)))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('LT1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
