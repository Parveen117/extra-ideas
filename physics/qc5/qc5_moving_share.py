"""QC5: a share that moves in scale - the winding reading and the logistic law.

Two flat readings on the quaternion carrier: the trivial one and the winding
one L = g^-1 dg, g = x/|x|.  A share p(r) between them gives A = p L and
  F = dp ^ L - p(1-p) L ^ L        (gradient term + variance term, QC4-T1 with a moving share).
Reduced action in s = log r:  3 Int [ p_s^2 + 4 p^2 (1-p)^2 ] ds
                            = 3 Int ( p_s - 2p(1-p) )^2 ds + 2 [3p^2 - 2p^3].
Exact rational quaternion arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
from itertools import permutations
import json

ZERO = (F(0), F(0), F(0), F(0))
UNIT = [tuple(F(int(i == k)) for i in range(4)) for k in range(4)]     # e_0 = 1, e_1, e_2, e_3


def qmul(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


def conj(a):
    return (a[0], -a[1], -a[2], -a[3])


def im(a):
    return (F(0), a[1], a[2], a[3])


def lin(*terms):
    out = [F(0)]*4
    for c, q in terms:
        for i in range(4):
            out[i] += c*q[i]
    return tuple(out)


def norm2(a):
    return sum(v*v for v in a)


def dot(a, b):
    return sum(u*v for u, v in zip(a, b))


def sign(perm):
    s, p = 1, list(perm)
    for i in range(4):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


EPS = {perm: sign(perm) for perm in permutations(range(4))}


def field(x, q, dq):
    """F_mu_nu for A_mu = q(rho) Im(xbar e_mu), rho = |x|^2, dq = dq/drho."""
    a = [im(qmul(conj(x), e)) for e in UNIT]
    Fmn = {}
    for m in range(4):
        for n in range(4):
            dA = lin((2*dq*x[m], a[n]), (-2*dq*x[n], a[m]),
                     (q, im(qmul(conj(UNIT[m]), UNIT[n]))), (-q, im(qmul(conj(UNIT[n]), UNIT[m]))))
            c = lin((q*q, qmul(a[m], a[n])), (-q*q, qmul(a[n], a[m])))
            Fmn[(m, n)] = lin((1, dA), (1, c))
    return Fmn


def dual(Fmn):
    out = {}
    for m in range(4):
        for n in range(4):
            acc = ZERO
            for r in range(4):
                for s in range(4):
                    if len({m, n, r, s}) == 4:
                        acc = lin((1, acc), (F(EPS[(m, n, r, s)], 2), Fmn[(r, s)]))
            out[(m, n)] = acc
    return out


def densities(Fmn):
    D = dual(Fmn)
    pairs = [(m, n) for m in range(4) for n in range(m+1, 4)]
    return (sum(norm2(Fmn[k]) for k in pairs), sum(dot(Fmn[k], D[k]) for k in pairs), D)


SHARES = {
    # name: (q, dq/drho) as functions of rho and the scale lam2 = lambda^2
    'logistic': lambda rho, l2: (1/(rho+l2), -1/(rho+l2)**2),
    'steeper': lambda rho, l2: (rho/(rho*rho+l2*l2), (l2*l2-rho*rho)/(rho*rho+l2*l2)**2),
    'equal': lambda rho, l2: (1/(2*rho), -1/(2*rho*rho)),
    'winding only': lambda rho, l2: (1/rho, -1/(rho*rho)),
}


def point_control(x, l2):
    rho = norm2(x)
    rows = {}
    for name, law in SHARES.items():
        q, dq = law(rho, l2)
        p = q*rho
        velocity = 2*rho*(dq*rho+q)                         # p_s = r dp/dr
        variance = 2*p*(1-p)
        Fmn = field(x, q, dq)
        action, topo, D = densities(Fmn)
        if action*rho*rho != 3*(velocity**2+variance**2):
            raise ValueError('action density is not 3 (velocity^2 + (2 p (1-p))^2) / r^4')
        if abs(topo)*rho*rho != 6*abs(velocity*variance):
            raise ValueError('topological density is not 6 velocity * variance / r^4')
        if action*rho*rho != 3*(velocity-variance)**2+6*velocity*variance:
            raise ValueError('square completion failed')
        selfdual = all(Fmn[k] == D[k] for k in Fmn) or all(Fmn[k] == lin((-1, D[k])) for k in Fmn)
        if selfdual != (velocity == variance or velocity == -variance):
            raise ValueError('self-duality is not velocity = variance')
        rows[name] = dict(share=str(p), velocity=str(velocity), twice_variance=str(variance),
                          action_density_r4=str(action*rho*rho), self_dual=selfdual)
    if not rows['logistic']['self_dual'] or rows['steeper']['self_dual']:
        raise ValueError('logistic share must be self-dual and the steeper one must not')
    if rows['winding only']['action_density_r4'] != '0':
        raise ValueError('a pure reading must have no field')
    if rows['equal']['velocity'] != '0' or rows['equal']['action_density_r4'] != '3/4':
        raise ValueError('equal share must be pure variance with density 3/(4 r^4)')
    return rows


def lift_control(x, l2):
    """Logistic share as a record lift u = (x, lambda): F = du^dagger Q du / N, Q = 1 - u u^dagger / N."""
    rho = norm2(x)
    N = rho+l2
    q, dq = SHARES['logistic'](rho, l2)
    Fmn = field(x, q, dq)
    for m in range(4):
        for n in range(4):
            em, en = UNIT[m], UNIT[n]
            cut = lin((1, qmul(conj(em), en)), (-1, qmul(conj(en), em)))        # du^dag du, antisymmetrized
            through = lin((rho/N, qmul(conj(em), en)), (-rho/N, qmul(conj(en), em)))  # du^dag P du
            memory = lin((1/N, cut), (-1/N, through))
            if memory != Fmn[(m, n)]:
                raise ValueError('field is not the memory form of the complement cut')
            closed = lin((l2/(N*N), cut))
            if memory != closed:
                raise ValueError('closed form lambda^2 (..)/N^2 failed')
    action, _, _ = densities(Fmn)
    if action != 24*l2*l2/N**4:
        raise ValueError('action density is not 24 lambda^4/(r^2 + lambda^2)^4 in the native coefficient norm')
    return dict(point=[str(v) for v in x], scale_squared=str(l2), action_density=str(action))


def reduced_action():
    """On the logistic law p_s = 2p(1-p): 3 Int (p_s^2 + 4p^2(1-p)^2) ds = Int_0^1 12 p (1-p) dp."""
    # exact antiderivative 6p^2 - 4p^3 = 2 (3p^2 - 2p^3)
    charge = lambda p: 3*p*p-2*p*p*p
    total = 2*(charge(F(1))-charge(F(0)))
    half = charge(F(1, 2))
    if total != 2 or half != F(1, 2):
        raise ValueError('reduced action or half share value failed')
    # scale independence: the same total for every lambda (the share runs 0 -> 1 whatever the scale)
    samples = []
    for l2 in (F(1, 9), F(1), F(25, 4)):
        for rho in (F(1, 100), F(1), F(100)):
            p = rho/(rho+l2)
            samples.append(charge(p))
        if not samples[-3] < samples[-2] < samples[-1]:
            raise ValueError('charge is not increasing along the share')
    return dict(minimum_reduced_action='2', share_change='0 -> 1', half_share_charge='1/2',
                classical_normalization='reduced 2 = 4 pi^2 in the native coefficient norm = 8 pi^2 in the trace norm (named classical value)')


def run():
    pts = [(F(1, 2), F(-1, 3), F(2, 5), F(1)), (F(0), F(3, 4), F(0), F(-2)), (F(2), F(1), F(-1), F(1, 7))]
    return dict(points=[dict(x=[str(v) for v in x], laws=point_control(x, F(4, 9))) for x in pts],
                lift=[lift_control(x, l2) for x in pts for l2 in (F(4, 9), F(3))],
                reduced=reduced_action())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('QC5_RESULT.json', 'w'), indent=1)
    for name, row in res['points'][0]['laws'].items():
        print(name, row)
    print(res['lift'][0])
    print(res['reduced'])
