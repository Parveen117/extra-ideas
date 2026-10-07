"""QC2: the fluctuation scale as the unit of return.

Exact part (stdlib): in the fluctuation ensemble of scale kappa the intensive
offset acts as kappa d/dw; Wick-ordered readings carry a number N and a
content q; the return generator has spectrum kappa*q.  Fluid part (numpy +
CoolProp): the return in the entropy representation, where the scale is the
constant k_B.  Python 3.12.
"""
from fractions import Fraction as F
import json
import math
from math import comb, factorial
import sys


# ------------------------------------------------------ Gaussian moments
def moment(a, b, C):
    """E[w1^a w2^b] for a centred Gaussian with covariance C (sum over pairings)."""
    if (a+b) % 2:
        return F(0)
    total = F(0)
    for k in range(min(a, b)+1):
        if (a-k) % 2 or (b-k) % 2:
            continue
        i, j = (a-k)//2, (b-k)//2
        total += (F(factorial(a)*factorial(b), factorial(k)*factorial(i)*factorial(j)*2**(i+j))
                  * C[0][0]**i*C[1][1]**j*C[0][1]**k)
    return total


def expect(poly, C):
    return sum(c*moment(a, b, C) for (a, b), c in poly.items())


def pmul(p, q):
    out = {}
    for (a, b), c in p.items():
        for (d, e), f in q.items():
            out[(a+d, b+e)] = out.get((a+d, b+e), 0)+c*f
    return {k: v for k, v in out.items() if v}


def pdiff(p, i):
    out = {}
    for (a, b), c in p.items():
        n = (a, b)[i]
        if n:
            key = (a-1, b) if i == 0 else (a, b-1)
            out[key] = out.get(key, 0)+c*n
    return out


def derivative_control(H, kappa, degree):
    """E[p_i f] = kappa E[d_i f] with p = H w, for every monomial f up to the degree."""
    det = H[0][0]*H[1][1]-H[0][1]**2
    C = [[kappa*H[1][1]/det, -kappa*H[0][1]/det], [-kappa*H[0][1]/det, kappa*H[0][0]/det]]
    checked = 0
    for a in range(degree+1):
        for b in range(degree+1-a):
            f = {(a, b): F(1)}
            for i in range(2):
                p_i = {(1, 0): H[i][0], (0, 1): H[i][1]}
                if expect(pmul(p_i, f), C) != kappa*expect(pdiff(f, i), C):
                    raise ValueError('intensive offset is not kappa d/dw in the ensemble')
                checked += 1
    pairing = expect({(2, 0): H[0][0], (1, 1): 2*H[0][1], (0, 2): H[1][1]}, C)
    if pairing != 2*kappa:
        raise ValueError('pairing p.w does not carry two units')
    return dict(identities=checked, mean_pairing=str(pairing))


# --------------------------------- Wick-ordered readings in the probe plane
def cx(re, im=0):
    return (F(re), F(im))


def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def zmoment(a, b, kappa):
    """E[z^a zbar^b], z = y1 + i y2, y isotropic of variance kappa: via the y expansion."""
    total = (F(0), F(0))
    iso = [[kappa, F(0)], [F(0), kappa]]
    for r in range(a+1):
        for s in range(b+1):
            coeff = cmul((F(comb(a, r)*comb(b, s)), F(0)), cmul(_ipow(r), _ipow(-s)))
            m = moment(a-r+b-s, r+s, iso)
            total = (total[0]+coeff[0]*m, total[1]+coeff[1]*m)
    return total


def _ipow(n):
    return [cx(1), cx(0, 1), cx(-1), cx(0, -1)][n % 4]


def wick(a, b, kappa):
    """Wick-ordered z^a zbar^b: polynomial {(a', b'): coefficient}."""
    return {(a-k, b-k): F((-2)**k*factorial(k)*comb(a, k)*comb(b, k))*kappa**k for k in range(min(a, b)+1)}


def zexpect(poly_a, poly_b_conj, kappa):
    """E[P conj(Q)] for real-coefficient polynomials in (z, zbar)."""
    total = (F(0), F(0))
    for (a, b), c in poly_a.items():
        for (d, e), f in poly_b_conj.items():
            m = zmoment(a+e, b+d, kappa)                   # conj swaps the two degrees
            total = (total[0]+c*f*m[0], total[1]+c*f*m[1])
    return total


def number_operator(poly, kappa):
    """(z d_z + zbar d_zbar - 4 kappa d_z d_zbar) P: minus the relaxation generator kappa*Laplacian - y.grad."""
    out = {}
    for (a, b), c in poly.items():
        out[(a, b)] = out.get((a, b), 0)+c*(a+b)
        if a and b:
            out[(a-1, b-1)] = out.get((a-1, b-1), 0)-4*kappa*c*a*b
    return {k: v for k, v in out.items() if v}


def content_operator(poly):
    return {k: v*(k[0]-k[1]) for k, v in poly.items() if v*(k[0]-k[1])}


def wick_control(kappa, degree):
    index = [(a, b) for a in range(degree+1) for b in range(degree+1-a)]
    for (a, b) in index:
        P = wick(a, b, kappa)
        if number_operator(P, kappa) != {k: v*(a+b) for k, v in P.items() if a+b}:
            raise ValueError('Wick reading is not a number state')
        if content_operator(P) != {k: v*(a-b) for k, v in P.items() if a-b}:
            raise ValueError('Wick reading has no definite content')
        if any(k[0]-k[1] != a-b for k in P):
            raise ValueError('Wick ordering mixed contents')
        for (c, d) in index:
            val = zexpect(P, wick(c, d, kappa), kappa)
            want = F(factorial(a)*factorial(b))*(2*kappa)**(a+b) if (a, b) == (c, d) else F(0)
            if val != (want, F(0)):
                raise ValueError('Wick readings are not orthogonal with norm a! b! (2 kappa)^(a+b)')
    table = {}
    for (a, b) in index:
        table.setdefault(a+b, []).append(a-b)
    return dict(readings=len(index), unit_per_quantum=str(2*kappa),
                contents_by_number={str(n): sorted(q) for n, q in sorted(table.items())})


def exact():
    H = [[F(5, 2), F(-3, 4)], [F(-3, 4), F(7, 3)]]
    out = dict(derivative=[derivative_control(H, k, 5) for k in (F(1), F(2, 3))],
               wick=[wick_control(k, 4) for k in (F(1), F(2, 3))])
    # the scale is visible: Wick ordering at one unit is not centred in an ensemble of another
    if zexpect(wick(1, 1, F(1)), {(0, 0): F(1)}, F(2, 3)) == (F(0), F(0)):
        raise ValueError('unit of return is invisible')
    return out


# -------------------------------------------------------------- fluid part
def fluid():
    import numpy as np
    sys.path.insert(0, '../ph1')
    import ph1_real_fluid_holonomy as p
    from CoolProp import CoolProp as CP

    def entropy_form(name, kind=None):
        st = CP.AbstractState('HEOS', name)
        R = st.gas_constant()/st.molar_mass()

        def G(T, rho):
            v = 1/rho
            if kind == 'ideal':
                st.update(CP.DmassT_INPUTS, 1e-6, T)
                cv = st.cp0mass()-R
                P, pT, pv = R*T/v, R/v, -R*T/v**2
            else:
                st.update(CP.DmassT_INPUTS, rho, T)
                cv, P = st.cvmass(), st.p()
                pT = st.first_partial_deriv(CP.iP, CP.iT, CP.iDmass)
                pv = -rho*rho*st.first_partial_deriv(CP.iP, CP.iDmass, CP.iT)
            Tv_u = -(T*pT-P)/cv
            s_uu = -1/(T*T*cv)
            s_uv = -Tv_u/(T*T)
            s_vv = (pv+pT*Tv_u)/T-P*Tv_u/(T*T)
            return -np.array([[s_uu, s_uv], [s_uv, s_vv]])
        return G

    rows = []
    for name, cyc in (('CarbonDioxide', (310.0, 350.0, 200.0, 700.0)), ('Argon', (155.0, 200.0, 250.0, 800.0))):
        T, rho = p.loop(8000, *cyc)
        H, _, _ = p.real_fluid(name)
        energy = p.holonomy([H(float(t), float(r)) for t, r in zip(T, rho)])
        G = entropy_form(name)
        Gs = [G(float(t), float(r)) for t, r in zip(T, rho)]
        if min(np.linalg.eigvalsh(g)[0] for g in Gs) <= 0:
            raise ValueError('entropy form is not positive on the cycle')
        balance = (1/math.sqrt(Gs[0][0, 0]), 1/math.sqrt(Gs[0][1, 1]))
        entropy = p.holonomy(Gs, scale=balance)
        Gi = entropy_form(name, 'ideal')
        ideal = p.holonomy([Gi(float(t), float(r)) for t, r in zip(T, rho)], scale=balance)
        rows.append(dict(fluid=name, cycle=cyc, theta_energy=energy['theta_loop'],
                         theta_entropy=entropy['theta_loop'], theta_entropy_transport=entropy['theta_transport'],
                         theta_entropy_ideal_cv0=ideal['theta_loop'],
                         length_energy=energy['shape_length'], length_entropy=entropy['shape_length']))
    return rows


if __name__ == '__main__':
    out = dict(exact=exact(), fluid=fluid())
    json.dump(out, open('QC2_RESULT.json', 'w'), indent=1)
    print(out['exact']['derivative'])
    print(out['exact']['wick'][1])
    for r in out['fluid']:
        print(r)
