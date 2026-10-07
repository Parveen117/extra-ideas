"""MS1: more speed, less mass - the state law, the flow law and the combination law.

State   psi = (psi1, psi2): n = psi1^2 + psi2^2, j = psi1^2 - psi2^2, sigma = 2 psi1 psi2 (PR4).
  T1  v = j / n = 2p - 1,  memory = 1 - v^2 = (sigma / n)^2,  n = sigma / sqrt(memory).
Flow    static law d_x psi = g K psi (from (A d_x + g R) psi = 0, A R = -K):
  T2  psi(x) = Exp(G K) psi0, G the accumulated flip; j is constant; a light-like inflow (1, 0)
      becomes n = cosh 2G, sigma = sinh 2G, v = sech 2G, memory = tanh^2 2G.
Coins   rest turns add under composition of coin turns:
  T3  speed c12 = c1 c2 - s1 s2;  M = 2 sin(theta/2):  M12 = M1 E2 + M2 E1 <= M1 + M2,
      t = tan(theta/2) combine by (t1 + t2)/(1 - t1 t2);  N equal coins: M_N = M_1 chi_N(theta/2) <= N M_1.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json
import sys

sys.path.insert(0, '../pr1')
sys.path.insert(0, '../qc3')
import pr1_sector_speed_and_mass as pr1
import qc3_even_odd_arrow as qc3

ONE, H, K, R, mm, add, sc = pr1.ONE, pr1.H, pr1.K, pr1.R, pr1.mm, pr1.add, pr1.sc


def vec(m, v):
    return (m[0][0]*v[0]+m[0][1]*v[1], m[1][0]*v[0]+m[1][1]*v[1])


def state(psi):
    n = psi[0]**2+psi[1]**2
    j = psi[0]**2-psi[1]**2
    sigma = 2*psi[0]*psi[1]
    return n, j, sigma


def state_control():
    rows = []
    for psi in ((F(3), F(1, 2)), (F(1), F(1)), (F(2), F(0)), (F(1, 3), F(5))):
        n, j, sigma = state(psi)
        v = j/n
        p = psi[0]**2/n
        memory = 4*p*(1-p)
        if v != 2*p-1 or 1-v*v != memory or (sigma/n)**2 != memory:
            raise ValueError('state law failed')
        rows.append(dict(share=str(p), speed=str(v), memory=str(memory)))
    if rows[2]['memory'] != '0' or rows[2]['speed'] != '1' or rows[1]['speed'] != '0':
        raise ValueError('pure reading must be light-like; balanced state must be at rest')
    return rows


def boostK(b):
    """Exp(G K) with e^G = b: cosh G + sinh G K."""
    ch, sh = (b+1/b)/2, (b-1/b)/2
    return add(sc(ch, ONE), sc(sh, K)), ch, sh


def flow_control():
    if mm(H, R) != sc(F(-1), K):
        raise ValueError('A R = -K failed')
    rows = []
    psi0 = (F(1), F(0))                                       # light-like inflow, pure reading
    prev = None
    for b in (F(1), F(3, 2), F(2), F(3), F(5)):
        E, ch, sh = boostK(b)
        psi = vec(E, psi0)
        n, j, sigma = state(psi)
        ch2, sh2 = ch*ch+sh*sh, 2*ch*sh                       # cosh 2G, sinh 2G
        if (n, j, sigma) != (ch2, F(1), sh2):
            raise ValueError('static flow is not (cosh 2G, 1, sinh 2G)')
        v = j/n
        if v != 1/ch2 or 1-v*v != (sh2/ch2)**2:
            raise ValueError('speed is not sech 2G')
        # generator: (Exp(GK) - Exp(-GK))/2 = sinh G K, and the semigroup law
        Em = boostK(1/b)[0]
        if sc(F(1, 2), add(E, Em, -1)) != sc(sh, K):
            raise ValueError('K is not the generator of the static flow')
        if prev is not None and mm(boostK(b/prev)[0], boostK(prev)[0]) != E:
            raise ValueError('accumulated flips do not compose')
        if prev is not None and not v < rows[-1]['_v']:
            raise ValueError('speed must drop as the flip accumulates')
        rows.append(dict(exp_G=str(b), density=str(n), current=str(j), proper_density=str(sigma),
                         speed=str(v), memory=str(1-v*v), _v=v))
        prev = b
    # a general inflow keeps its current and obeys sigma = j sqrt(1 - v^2) / v in squares
    psi0 = (F(3), F(1))
    for b in (F(2), F(1, 3)):
        psi = vec(boostK(b)[0], psi0)
        n, j, sigma = state(psi)
        if j != state(psi0)[1] or sigma*sigma*j*j != j*j*(n*n-j*j):
            raise ValueError('current is not constant along the static flow')
        v = j/n
        if (sigma/j)**2 != (1-v*v)/(v*v):
            raise ValueError('proper density is not j sqrt(1 - v^2)/v')
    for r in rows:
        del r['_v']
    return rows


def combination_control():
    halves = [(F(4, 5), F(3, 5)), (F(12, 13), F(5, 13)), (F(40, 41), F(9, 41)), (F(15, 17), F(8, 17))]
    rows = []
    for i, (E1, O1) in enumerate(halves):
        for (E2, O2) in halves[i:]:
            c1, s1 = E1*E1-O1*O1, 2*E1*O1
            c2, s2 = E2*E2-O2*O2, 2*E2*O2
            c12, s12 = c1*c2-s1*s2, s1*c2+c1*s2
            E12, O12 = E1*E2-O1*O2, O1*E2+E1*O2              # half turn of the composite
            if (E12*E12-O12*O12, 2*E12*O12) != (c12, s12):
                raise ValueError('composite half turn failed')
            M1, M2, M12 = 2*O1, 2*O2, 2*O12
            if M12 != M1*E2+M2*E1:
                raise ValueError('mass combination law failed')
            defect = M1+M2-M12
            if defect != M1*(1-E2)+M2*(1-E1) or defect <= 0:
                raise ValueError('mass defect is not M1 (1 - E2) + M2 (1 - E1)')
            t1, t2 = O1/E1, O2/E2
            if E12 != 0 and O12/E12 != (t1+t2)/(1-t1*t2):
                raise ValueError('circular addition law failed')
            if c12 >= min(c1, c2):
                raise ValueError('combining rest turns must reduce the speed')
            walk = pr1.walk_control(c12, s12)                 # the composite coin is a coin
            if walk['cone_speed'] != str(c12):
                raise ValueError('composite walk speed failed')
            rows.append(dict(speeds=[str(c1), str(c2)], composite_speed=str(c12),
                             masses=[str(M1), str(M2)], composite_mass=str(M12), mass_defect=str(defect)))
    return rows


def tower_control(E1, O1, count=6):
    """N equal coins: half turn (E_N, O_N); M_N = M_1 chi_N, chi_N = U_{N-1}(E_1) <= N."""
    rows = []
    E, O = F(1), F(0)
    stopped = None
    for N in range(1, count+1):
        E, O = E*E1-O*O1, O*E1+E*O1
        c = E*E-O*O
        chi = qc3.cheb_U(N-1, E1)
        if 2*O != 2*O1*chi:
            raise ValueError('M_N = M_1 chi_N failed')
        if chi > N or (chi == N) != (N == 1):
            raise ValueError('circular character must stay below N')
        if stopped is None and c <= 0:
            stopped = N
        rows.append(dict(coins=N, speed=str(c), mass=str(2*O), binding=str(2*O1*(N-chi))))
    return dict(one_coin_speed=str(E1*E1-O1*O1), rows=rows, speed_reaches_zero_by=stopped)


def run():
    return dict(state=state_control(), flow=flow_control(), combination=combination_control(),
                tower=tower_control(F(40, 41), F(9, 41)))


if __name__ == '__main__':
    res = run()
    json.dump(res, open('MS1_RESULT.json', 'w'), indent=1)
    print(res['state'])
    for r in res['flow']:
        print(r)
    print(res['combination'][0]); print(res['combination'][2])
    for r in res['tower']['rows']:
        print({k: (v if len(str(v)) < 22 else round(float(F(v)), 5)) for k, v in r.items()})
    print('stops by', res['tower']['speed_reaches_zero_by'])
