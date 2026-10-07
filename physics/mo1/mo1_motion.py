"""MO1: motion in the memory field - the count of a moving reading and its largest-count histories.

Setting (the thesis plus one assumption, stated in the document):
  at radius r the locally pure (flat) frame of GR2-G1 moves inward with beta^2 = m = r_s / r, and
  these frames share one flat space and one time t.  The count of a leg (CL2-T1) is then g d(tau),
        d(tau)^2 = dt^2 - (dr + beta dt)^2 - r^2 d(phi)^2 .
M1  Same form in static time dT = dt - beta dr / N^2:  N^2 dT^2 - dr^2 / N^2 - r^2 d(phi)^2,  N^2 = 1 - m.
M2  A reading at rest counts N dt (GR1's clock factor); one moving with the frame counts dt.
M3  Conserved E, L and the radial law  rdot^2 = E^2 - N^2 (1 + L^2 / r^2).
M4  Circular histories: (d phi/dt)^2 = r_s / (2 r^3);  count factor^2 = 1 - 3m/2;  light circles at m = 2/3;
    the last stable circle at r = 3 r_s with L^2 = 3 r_s^2.
M5  Orbit equation u'' + u = r_s/(2 L^2) + (3/2) r_s u^2 (u = 1/r); near-circular: (radial/angular)^2 = 1 - 3m.
M6  Light: first-order solution and deflection 2 r_s / b.
Exact rational arithmetic for M1-M6; floats only in the illustration.  Python 3.12.
"""
from fractions import Fraction as F
import json
import math

PYTH = [(F(3, 5), F(4, 5)), (F(4, 5), F(3, 5)), (F(5, 13), F(12, 13)), (F(8, 17), F(15, 17))]   # (beta, N)


def form_control():
    rows = []
    for beta, N in PYTH:
        if N*N != 1-beta*beta:
            raise ValueError('N^2 + memory = 1 failed')
        for dt, dr, dphi, r in ((F(1), F(1, 3), F(1, 7), F(5)), (F(2), F(-1, 2), F(0), F(9)), (F(1), F(0), F(1, 4), F(3))):
            river = dt*dt-(dr+beta*dt)**2-r*r*dphi*dphi
            dT = dt-beta*dr/(N*N)
            static = N*N*dT*dT-dr*dr/(N*N)-r*r*dphi*dphi
            if river != static:
                raise ValueError('river and static forms of the count differ')
        at_rest = F(1)-(beta*1)**2                                 # dt = 1, dr = dphi = 0
        if at_rest != N*N:
            raise ValueError('a reading at rest must count N dt')
        with_frame = F(1)-(-beta+beta)**2                          # dr = -beta dt
        if with_frame != 1:
            raise ValueError('a reading moving with the frame must count dt')
        rows.append(dict(memory=str(beta*beta), N=str(N), count_at_rest=str(N), count_with_the_frame='1'))
    return rows


def radial_law_control():
    checked = 0
    for beta, N in PYTH:
        for r, L, rdot, E in ((F(5), F(2), F(1, 3), F(3, 2)), (F(9), F(0), F(-1, 2), F(1)), (F(4), F(7, 2), F(0), F(2))):
            # u = rdot + beta tdot;  E = tdot - beta u;  solve the two linear relations for tdot, u
            # tdot = E + beta u,  rdot = u - beta tdot  ->  u = (rdot + beta E) / N^2
            u = (rdot+beta*E)/(N*N)
            tdot = E+beta*u
            if rdot != u-beta*tdot or E != tdot-beta*u:
                raise ValueError('conserved quantity relations failed')
            norm = tdot*tdot-u*u-L*L/(r*r)                          # d(tau)^2 per unit parameter
            # the radial law in the form valid for any normalisation:
            if N*N*norm != E*E-rdot*rdot-N*N*L*L/(r*r):
                raise ValueError('radial law failed')
            checked += 1
    return checked


def circular_control():
    rows = []
    for m, r in ((F(1, 9), F(18)), (F(3, 25), F(5)), (F(2, 9), F(9)), (F(1, 3), F(6)), (F(2, 3), F(3))):
        r_s = m*r
        # radial stationarity for a circle (rdot = 0, tdot = 1): -beta beta' - r omega^2 = 0, beta beta' = -m/(2r)
        def residual(w2):
            return m/(2*r)-r*w2
        omega2 = r_s/(2*r**3)
        if residual(omega2) != 0 or residual(omega2*2) == 0:
            raise ValueError('circular law failed')
        factor2 = 1-m-r*r*omega2
        if factor2 != 1-F(3, 2)*m:
            raise ValueError('count factor on a circle failed')
        rows.append(dict(memory=str(m), omega_squared=str(omega2), count_factor_squared=str(factor2),
                         light_like=(factor2 == 0)))
    if not rows[-1]['light_like'] or any(r['light_like'] for r in rows[:-1]):
        raise ValueError('light must circle exactly at memory 2/3')
    # last stable circle: V = (1 - r_s/r)(1 + L^2/r^2), V' = V'' = 0 at r = 3 r_s, L^2 = 3 r_s^2 (units r_s = 1)
    r, L2 = F(3), F(3)
    V1 = 1/r**2-2*L2/r**3+3*L2/r**4
    V2 = -2/r**3+6*L2/r**4-12*L2/r**5
    if V1 != 0 or V2 != 0:
        raise ValueError('last stable circle is not at r = 3 r_s')
    r = F(4)
    L2 = r*r/(2*r-3)                                               # circular condition V' = 0
    if 1/r**2-2*L2/r**3+3*L2/r**4 != 0 or -2/r**3+6*L2/r**4-12*L2/r**5 <= 0:
        raise ValueError('a circle outside 3 r_s must be stable')
    return rows


def orbit_control():
    rows = []
    for r_s, L, E, u in ((F(1), F(3), F(1), F(1, 7)), (F(1, 2), F(2), F(9, 10), F(1, 5)), (F(2), F(5), F(1), F(1, 20))):
        L2 = L*L
        # F(u) = (E^2 - (1 - r_s u)(1 + L^2 u^2)) / L^2,  u'' = F'(u) / 2
        dF = (r_s*(1+L2*u*u)-2*L2*u*(1-r_s*u))/L2
        if dF/2 != r_s/(2*L2)-u+F(3, 2)*r_s*u*u:
            raise ValueError('orbit equation failed')
        rows.append(dict(r_s=str(r_s), L=str(L)))
    # near-circular: epsilon'' + (1 - 3 r_s u0) epsilon = 0
    near = []
    for m, root in ((F(3, 25), F(4, 5)), (F(1, 4), F(1, 2)), (F(5, 27), F(2, 3))):
        if root*root != 1-3*m:
            raise ValueError('radial/angular ratio failed')
        near.append(dict(memory=str(m), ratio=str(root), advance_per_turn_over_2pi=str(1/root-1), first_order=str(F(3, 2)*m)))
    return dict(orbit_equation_points=len(rows), near_circular=near)


def light_control():
    # first-order solution of u'' + u = (3/2) r_s u^2 about u0 = sin(phi)/b:
    # u1 = (r_s / (2 b^2)) (1 + cos^2 phi);  check with rational (cos, sin)
    for c, s in ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(1), F(0)), (F(0), F(1))):
        cos2 = 2*c*c-1
        left = (1+c*c)-2*cos2                                      # (u1'' + u1) * 2 b^2 / r_s
        right = 3*s*s
        if left != right:
            raise ValueError('first-order light solution failed')
    # asymptote: u0 + u1 = 0 at phi = -delta: -delta/b + (r_s/(2 b^2)) * 2 = 0
    out = []
    for r_s, b in ((F(1), F(100)), (F(3), F(10**6))):
        delta = r_s/b
        if -delta/b+(r_s/(2*b*b))*2 != 0:
            raise ValueError('asymptote failed')
        out.append(dict(r_s=str(r_s), b=str(b), deflection=str(2*delta), clock_only_value=str(delta)))
    return out


def illustration():
    G_N, c, Msun = 6.67430e-11, 299792458.0, 1.98847e30
    r_s = 2*G_N*Msun/c**2
    a, e = 5.7909e10, 0.20563
    per_orbit = 3*math.pi*r_s/(a*(1-e*e))
    arcsec = 180*3600/math.pi
    return dict(sun_r_s_m=r_s, mercury_advance_arcsec_per_century=per_orbit*arcsec*(36525/87.969),
                light_at_solar_limb_arcsec=2*r_s/6.957e8*arcsec,
                earth_orbit_count_factor_minus_one=-0.75*r_s/1.496e11)


def run():
    return dict(form=form_control(), radial_law_points=radial_law_control(), circular=circular_control(),
                orbit=orbit_control(), light=light_control(), illustration=illustration())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('MO1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
