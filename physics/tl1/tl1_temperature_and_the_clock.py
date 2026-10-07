"""TL1: heat and the clock factor -- the unit of fluctuation is carried by the clock.

Sources: R43.3 (exchange conserves the summed readout iff the unit factors are equal), DC1-K1 (re-certified),
GR1-T2 / CL1-T2 (a reading at clock factor N counts g N per unit of static time), PR3 (accelerated frame: the
law is the inertial one with d_t -> (1/rho) d_eta), QC2 (fluctuation unit kappa: weights exp(-energy/kappa)).
L1  A turn kept in static time arrives with local rate  w_B = w_A N_A / N_B.
L2  Two ledgers at places A, B that exchange through static time conserve  N_A E_A + N_B E_B.  With local unit
    factors u_A, u_B the summed readout is kept for every preparation  iff  u_A N_A = u_B N_B   (R43.3 with
    b = u N).
L3  Two fluctuation ensembles exchanging quanta of static size e have no net flow  iff  kappa_A N_A = kappa_B N_B.
    Otherwise the flow runs toward the smaller kappa N.
L4  Hence in equilibrium  kappa(r) N(r) = constant.   Weak field: d ln kappa / dh = g / c^2.
    Accelerated frame (PR3): kappa * rho = constant.
L5  A ledger of local energy U at clock factor N has static energy N U; its weight is U g / c^2.
Exact rational arithmetic; floats only for the illustration.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../dc1')
import dc1_memory_between_masses as dc

c, ident, kron, madd, scal, mmul, msub, frob2, swap, trace = dc.c, dc.ident, dc.kron, dc.madd, dc.scal, dc.mmul, dc.msub, dc.frob2, dc.swap, dc.trace


def clock(num, den):
    """a rational clock factor N = num/den, with memory m = 1 - N^2 (GR1)."""
    n = F(num, den)
    if not 0 < n <= 1:
        raise ValueError('clock factor must lie in (0, 1]')
    return n


def redshift_control():
    rows = []
    for na, nb, w in ((clock(3, 5), clock(4, 5), F(7)), (clock(5, 13), clock(1, 1), F(2, 3)), (clock(12, 13), clock(3, 5), F(11, 4))):
        static = w*na                       # rate per unit of static time (GR1/CL1)
        wb = static/nb                      # read locally at B
        if wb*nb != w*na or wb != w*na/nb:
            raise ValueError('transport of a turn failed')
        rows.append((str(na), str(nb), str(w), str(wb)))
    return rows


def unit_control():
    E = [[c(2), c(1, 1)], [c(1, -1), c(-1)]]
    d = 2
    S = swap(d)
    bracket = trace(mmul(E, E))[0]-trace(E)[0]**2/d
    rows = []
    for na, nb in ((clock(3, 5), clock(4, 5)), (clock(5, 13), clock(12, 13))):
        for ua, ub in ((F(1), F(1)), (nb, na), (F(4), F(3)), (2*nb, 2*na)):
            ba, bb = ua*na, ub*nb
            C = madd(scal(ba, kron(E, ident(d))), scal(bb, kron(ident(d), E)))
            comm = msub(mmul(S, C), mmul(C, S))
            if frob2(comm) != 2*d*(ba-bb)**2*bracket:
                raise ValueError('calibration identity with clock factors failed')
            kept = frob2(comm) == 0
            if kept != (ua*na == ub*nb):
                raise ValueError('units must satisfy u_A N_A = u_B N_B exactly when the sum is kept')
            rows.append((str(na), str(nb), str(ua), str(ub), kept))
    if not any(r[4] and r[2] != r[3] for r in rows):
        raise ValueError('a kept case with unequal local units must be present')
    return rows


def net_flow(xa, xb, levels=6):
    """two ledgers with geometric weights xa^i, xb^j (x = exp(-e / (kappa N)) per static quantum);
    one quantum moves A->B at rate proportional to the occupation of the giving ledger's upper level.
    Returns the net A->B current of the product weight (zero iff detailed balance)."""
    za = sum(xa**i for i in range(levels))
    zb = sum(xb**j for j in range(levels))
    cur = F(0)
    for i in range(levels):
        for j in range(levels):
            w = xa**i*xb**j/(za*zb)
            if i >= 1 and j+1 < levels:
                cur += w               # A gives one quantum
            if j >= 1 and i+1 < levels:
                cur -= w               # B gives one quantum
    return cur


def equilibrium_control():
    rows = []
    for xa, xb in ((F(1, 2), F(1, 2)), (F(2, 3), F(2, 3)), (F(1, 2), F(1, 3)), (F(1, 4), F(3, 5))):
        cur = net_flow(xa, xb)
        if (cur == 0) != (xa == xb):
            raise ValueError('no net flow must hold exactly for equal weights per static quantum')
        if xa != xb and (cur > 0) != (xa > xb):
            raise ValueError('flow must run from the larger kappa N to the smaller')
        rows.append((str(xa), str(xb), str(cur)))
    return rows


def weight_control():
    """static energy N U with N^2 = 1 - r_s/r.  Exact difference quotient of (N U)^2 in r:
    [ (N U)^2(r+h) - (N U)^2(r-h) ] / 2h = U^2 r_s / (r^2 - h^2)  ->  U^2 r_s / r^2 ,
    so d(N U)/dr = U r_s / (2 r^2 N): the pull on a local energy U is that on a mass U / c^2."""
    rows = []
    for rs, r, U, Nn in ((F(9, 25), F(1), F(5), F(4, 5)), (F(32, 25), F(2), F(3), F(3, 5))):
        if Nn*Nn != 1-rs/r:
            raise ValueError('clock factor fixture failed')
        sq = lambda rr: (1-rs/rr)*U*U
        for h in (F(1, 10), F(1, 1000)):
            if (sq(r+h)-sq(r-h))/(2*h) != U*U*rs/(r*r-h*h):
                raise ValueError('difference quotient failed')
        slope = U*rs/(2*r*r*Nn)
        newton = (U/1)*(rs/(2*r*r))/Nn           # mass U (c = 1) times r_s/(2 r^2), over N
        if slope != newton:
            raise ValueError('weight of a local energy failed')
        rows.append((str(rs), str(r), str(U), str(slope)))
    return rows


def illustration():
    g, cl, G_N, M_sun, R_sun = 9.80665, 299792458.0, 6.67430e-11, 1.98841e30, 6.957e8
    per_metre = g/cl**2
    sun = G_N*M_sun/(R_sun*cl**2)
    return dict(earth_fraction_per_metre=per_metre, earth_kelvin_per_km_at_300K=300*per_metre*1000,
                sun_surface_to_far_fraction=sun, weight_of_one_joule_newton=g/cl**2)


def run():
    return dict(redshift=redshift_control(), unit=unit_control(), equilibrium=equilibrium_control(),
                weight=weight_control(), illustration=illustration())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('TL1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v)
