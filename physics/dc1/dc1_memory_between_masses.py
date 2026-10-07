"""DC1: the common unit, and gravitational memory between two masses.

Reads the owner's R43 (exchange, relative calibration, pair memory) and the information-invariance
ledger's D4 (decoherence) with the gravity thesis.
K1  (R43.3 re-certified)  For C = bA E(x)1 + bB 1(x)E and the swap S:
        |[S, C]|^2 = 2 d (bA - bB)^2 [ Tr E^2 - (Tr E)^2 / d ] ,
    so exchange conserves the summed readout for every preparation iff bA = bB (E not scalar):
    two frames that exchange share one unit.  The common value stays free (R43.7, SC1).
K2  Two two-branch ledgers A, B with shares p, q and a branch-pair phase turn Delta
    (Delta = phi11 - phi12 - phi21 + phi22):
        memory of A  M = 1 - Tr rho_A^2 = (1/2) [4p(1-p)] [4q(1-q)] sin^2(Delta/2) ,
        IN1's invariant of A's reading  n^2 - r.r = 2 M .
K3  With the thesis clock factor (GR1, weak field) a reading of rest turn g1 at distance d from mass 2
    counts g1 (1 - r_s2 / (2 d)) per unit time, so phi = G m1 m2 tau / (hbar d).
K4  One ledger alone (q = 0 or 1, or Delta = 0): M = 0 at all times.
Exact arithmetic over the Gaussian rationals for K1, K2, K4; floats only for the illustration.  Python 3.12.
"""
from fractions import Fraction as F
import json
import math
import sys

sys.path.insert(0, '../in1')
import in1_the_invariant as in1

c, cmul, cadd, conj = in1.c, in1.cmul, in1.cadd, in1.conj
ZC = (F(0), F(0))


def kron(a, b):
    n, m = len(a), len(b)
    return [[cmul(a[i//m][j//m], b[i % m][j % m]) for j in range(n*m)] for i in range(n*m)]


def mmul(a, b):
    n = len(a)
    out = [[ZC]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = ZC
            for k in range(n):
                s = cadd(s, cmul(a[i][k], b[k][j]))
            out[i][j] = s
    return out


def msub(a, b):
    return [[(x[0]-y[0], x[1]-y[1]) for x, y in zip(r1, r2)] for r1, r2 in zip(a, b)]


def scal(z, a):
    return [[cmul(c(z), x) for x in row] for row in a]


def madd(a, b):
    return [[cadd(x, y) for x, y in zip(r1, r2)] for r1, r2 in zip(a, b)]


def frob2(a):
    return sum(x[0]**2+x[1]**2 for row in a for x in row)


def trace(a):
    s = ZC
    for i in range(len(a)):
        s = cadd(s, a[i][i])
    return s


def ident(d):
    return [[c(int(i == j)) for j in range(d)] for i in range(d)]


def swap(d):
    S = [[ZC]*(d*d) for _ in range(d*d)]
    for i in range(d):
        for j in range(d):
            S[j*d+i][i*d+j] = c(1)
    return S


def unit_control():
    rows = []
    for E in ([[c(2), c(1, 1)], [c(1, -1), c(-1)]], [[c(1), c(0), c(0, 2)], [c(0), c(3), c(1)], [c(0, -2), c(1), c(F(1, 2))]]):
        d = len(E)
        S = swap(d)
        bracket = (trace(mmul(E, E))[0]-trace(E)[0]**2/d)
        if bracket <= 0:
            raise ValueError('readout must not be scalar')
        for bA, bB in ((F(1), F(1)), (F(3), F(1)), (F(2, 3), F(5, 4))):
            C = madd(scal(bA, kron(E, ident(d))), scal(bB, kron(ident(d), E)))
            comm = msub(mmul(S, C), mmul(C, S))
            if frob2(comm) != 2*d*(bA-bB)**2*bracket:
                raise ValueError('R43.3 calibration identity failed')
            if (frob2(comm) == 0) != (bA == bB):
                raise ValueError('conservation must hold exactly for equal units')
        rows.append(dict(roles=d, bracket=str(bracket)))
    return rows


def pair_state(a, b, zeta):
    """amplitudes (a0, a1) for A, (b0, b1) for B; phase zeta on the branch pair (1, 1) only."""
    amp = [[cmul(c(a[i]), c(b[j])) for j in range(2)] for i in range(2)]
    amp[1][1] = cmul(amp[1][1], zeta)
    return amp


def reduced(amp):
    rho = [[ZC, ZC], [ZC, ZC]]
    for i in range(2):
        for k in range(2):
            s = ZC
            for j in range(2):
                s = cadd(s, cmul(amp[i][j], conj(amp[k][j])))
            rho[i][k] = s
    return rho


def memory_control():
    rows = []
    for a in ((F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(1), F(0))):
        for b in ((F(3, 5), F(4, 5)), (F(8, 17), F(15, 17)), (F(0), F(1))):
            for cz, sz in ((F(3, 5), F(4, 5)), (F(-7, 25), F(24, 25)), (F(1), F(0))):
                p, q = a[0]**2, b[0]**2
                rho = reduced(pair_state(a, b, (cz, sz)))
                tr = trace(rho)[0]
                purity = trace(mmul(rho, rho))[0]
                M = 1-purity
                sin2_half = (1-cz)/2
                want = F(1, 2)*(4*p*(1-p))*(4*q*(1-q))*sin2_half
                if tr != 1 or M != want:
                    raise ValueError('pair memory law failed')
                n, r = in1.readings(rho)
                if in1.form(n, r) != 2*M:
                    raise ValueError('reading invariant is not twice the memory')
                rows.append((str(p), str(q), str(cz), str(M)))
    alone = [row for row in rows if row[1] in ('0', '1') or row[2] == '1']
    if any(row[3] != '0' for row in alone):
        raise ValueError('one ledger alone, or no relative phase, must leave no memory')
    best = max(F(row[3]) for row in rows)
    return dict(cases=len(rows), largest_memory=str(best), isolated_or_no_phase_memory='0',
                example=dict(p='9/25', q='9/25', cos_delta='3/5',
                             memory=str(F(1, 2)*(4*F(9, 25)*F(16, 25))**2*F(1, 5))))


def phase_illustration():
    G_N, hbar, cl = 6.67430e-11, 1.054571817e-34, 299792458.0

    def delta(m1, m2, tau, d11, d12, d21, d22):
        return G_N*m1*m2*tau/hbar*(1/d11-1/d12-1/d21+1/d22)
    # two 1e-14 kg masses, each split by 250 micrometres along the line joining them, centres 450 micrometres apart
    L, D, m, tau = 250e-6, 450e-6, 1e-14, 2.5
    d = delta(m, m, tau, D, D+L, D-L, D)
    M = 0.5*math.sin(d/2)**2
    # the owner's example: 1e-17 kg, 1e-7 m
    m2, R = 1e-17, 1e-7
    scale = G_N*m2*m2/(hbar*R)
    rows = []
    for t in (10.0, 100.0, 1000.0):
        x = scale*t
        rows.append(dict(seconds=t, phase=x, pair_memory=0.5*math.sin(x/2)**2,
                         exponential_law_deficit=0.5*(1-math.exp(-2*x))))
    # clock factor check: g1 (1 - r_s2/(2d)) -> phase rate G m1 m2 / (hbar d)
    g1 = m*cl*cl/hbar
    r_s2 = 2*G_N*m/cl**2
    return dict(two_masses=dict(mass_kg=m, split_m=L, centres_m=D, seconds=tau, delta_rad=d, memory=M),
                clock_factor_phase_rate=g1*r_s2/(2*D), direct_phase_rate=G_N*m*m/(hbar*D),
                owner_example=dict(mass_kg=m2, length_m=R, rate_scale_per_s=scale, rows=rows))


def run():
    return dict(unit=unit_control(), memory=memory_control(), illustration=phase_illustration())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('DC1_RESULT.json', 'w'), indent=1)
    print(res['unit']); print(res['memory'])
    ill = res['illustration']
    print(ill['two_masses']); print(ill['clock_factor_phase_rate'], ill['direct_phase_rate'])
    print(ill['owner_example']['rate_scale_per_s'])
    for r in ill['owner_example']['rows']:
        print(r)
