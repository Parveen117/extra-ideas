"""CL1: the native clock answers the clock question - mass^2 is clock curvature.

Docked to R41 (relational clock): tick T (mark successor), phase readout D e_r = zeta^r e_r,
count N e_r = r e_r;  R41.3: D T D^-1 T^-1 = zeta,  [D,T]^dagger [D,T] = 2 - zeta - zeta^-1.
Docked to PR1 (coin-turn walk U = Shift . Exp(theta R)).
T1  clock curvature of a tick phase zeta = c + iota s is 2 - 2c.
T2  the coin turn with the same zeta has rest defect 2 - (U + U^-1) = 2 - 2c at the uniform
    sector: mass^2 = clock curvature;  cone speed = 1 - (clock curvature)/2.
T3  [N, O_T] = E_T and [N, E_T] = O_T for the tick's even and odd parts (away from a wrap).
T4  position-dependent coin: (U + U^-1) = cos-weighted shifts + (s_x - s_{x-1}) exchange.
T5  history residue of Q ticks of a turn: |W^Q - I|^2 / Q = 2 (1 - E(Q phi)) / Q.
Exact rational arithmetic, stdlib only.  Python 3.12.
"""
from fractions import Fraction as F
import json


def cmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def cpow(z, n):
    out = (F(1), F(0))
    base = z if n >= 0 else (z[0], -z[1])
    for _ in range(abs(n)):
        out = cmul(out, base)
    return out


# ------------------------------------------------------------ T1: the Weyl pair
def weyl_control(c, s, marks=6):
    if c*c+s*s != 1:
        raise ValueError('not a turn')
    zeta = (c, s)
    for r in range(marks):
        # D T e_r = zeta^(r+1) e_(r+1);  T D e_r = zeta^r e_(r+1)
        dt, td = cpow(zeta, r+1), cpow(zeta, r)
        if dt != cmul(zeta, td):
            raise ValueError('D T D^-1 T^-1 = zeta failed')
        diff = (dt[0]-td[0], dt[1]-td[1])                   # [D, T] e_r coefficient
        if diff[0]**2+diff[1]**2 != 2-2*c:
            raise ValueError('clock curvature is not 2 - zeta - zeta^-1')
    return 2-2*c


# ---------------------------------------- T2 / T4: the walk with a coin field
def walk_sum(coins, psi):
    """(U + U^-1) psi on interior sites; coins[x] = (c_x, s_x); psi[x] = (a_x, b_x).

    U psi (x) = ( [c a - s b](x-1), [s a + c b](x+1) );  U^-1 phi (x) = C_x^-1 (phi_1(x+1), phi_2(x-1)).
    """
    n = len(coins)
    U = {}
    for x in range(1, n-1):
        cl, sl = coins[x-1]
        cr, sr = coins[x+1]
        U[x] = (cl*psi[x-1][0]-sl*psi[x-1][1], sr*psi[x+1][0]+cr*psi[x+1][1])
    Ui = {}
    for x in range(1, n-1):
        c, s = coins[x]
        p1, p2 = psi[x+1][0], psi[x-1][1]
        Ui[x] = (c*p1+s*p2, -s*p1+c*p2)
    return {x: (U[x][0]+Ui[x][0], U[x][1]+Ui[x][1]) for x in U}


def uniform_control(c, s):
    coins = [(c, s)]*7
    out = []
    for comp in (0, 1):
        psi = [(F(int(comp == 0)), F(int(comp == 1)))]*7     # the uniform sector
        got = walk_sum(coins, psi)[3]
        want = (2*c*F(int(comp == 0)), 2*c*F(int(comp == 1)))
        if got != want:
            raise ValueError('uniform sector: U + U^-1 = 2 cos(theta) failed')
        out.append(2-got[comp])
    if out[0] != out[1]:
        raise ValueError('rest defect differs between components')
    return out[0]


def tower_control():
    """Tick phases zeta_Q for Q = 2, 4 (rational) and the squares for Q = 8: curvature, speed."""
    rows = []
    for Q, (c, s) in ((2, (F(-1), F(0))), (4, (F(0), F(1)))):
        curv = weyl_control(c, s)
        if cpow((c, s), Q) != (F(1), F(0)):
            raise ValueError('tick phase does not close after Q marks')
        rows.append(dict(marks=Q, clock_curvature=str(curv), cone_speed=str(1-curv/2)))
    # Q = 8: c^2 = 1/2 (positive half-root of the quarter turn): (2 - curvature)^2 = 4 c^2 = 2
    rows.append(dict(marks=8, clock_curvature='2 - sqrt(2)', cone_speed_squared='1/2',
                     relation='(2 - curvature)^2 = 4 speed^2 = 2'))
    return rows


def varying_coin_control():
    triples = [(F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(4, 5), F(3, 5)), (F(1), F(0)),
               (F(8, 17), F(15, 17)), (F(12, 13), F(5, 13)), (F(0), F(1))]
    psi = [(F(k+1), F(2*k-3, 2)) for k in range(7)]
    got = walk_sum(triples, psi)
    for x in range(1, 6):
        c0, s0 = triples[x-1]
        c1, s1 = triples[x]
        c2, s2 = triples[x+1]
        want1 = c0*psi[x-1][0]+c1*psi[x+1][0]+(s1-s0)*psi[x-1][1]
        want2 = c2*psi[x+1][1]+c1*psi[x-1][1]+(s2-s1)*psi[x+1][0]
        if got[x] != (want1, want2):
            raise ValueError('varying-coin identity failed')
    flat = walk_sum([triples[0]]*7, psi)
    c, s = triples[0]
    for x in range(1, 6):
        if flat[x] != (c*(psi[x-1][0]+psi[x+1][0]), c*(psi[x-1][1]+psi[x+1][1])):
            raise ValueError('constant coin should have no exchange term')
    return dict(exchange_term='(s_x - s_{x-1}) psi_2(x-1) in component 1, (s_{x+1} - s_x) psi_1(x+1) in component 2',
                constant_coin='no exchange term')


# ------------------------------------------------- T3: count against the tick
def count_control(marks=7):
    """N, T as matrices on a window of marks; identities checked on interior basis vectors."""
    T = [[F(int(i == j+1)) for j in range(marks)] for i in range(marks)]       # T e_j = e_(j+1)
    Ti = [[F(int(i == j-1)) for j in range(marks)] for i in range(marks)]      # T^-1 e_j = e_(j-1)
    N = [[F(i) if i == j else F(0) for j in range(marks)] for i in range(marks)]

    def mul(a, b):
        return [[sum(a[i][k]*b[k][j] for k in range(marks)) for j in range(marks)] for i in range(marks)]

    def lin(ca, a, cb, b):
        return [[ca*a[i][j]+cb*b[i][j] for j in range(marks)] for i in range(marks)]
    E = lin(F(1, 2), T, F(1, 2), Ti)
    O = lin(F(1, 2), T, F(-1, 2), Ti)
    NO = lin(F(1), mul(N, O), F(-1), mul(O, N))
    NE = lin(F(1), mul(N, E), F(-1), mul(E, N))
    for j in range(1, marks-1):
        if [NO[i][j] for i in range(marks)] != [E[i][j] for i in range(marks)]:
            raise ValueError('[N, O] = E failed')
        if [NE[i][j] for i in range(marks)] != [O[i][j] for i in range(marks)]:
            raise ValueError('[N, E] = O failed')
    return dict(count_with_odd='even part', count_with_even='odd part')


# ------------------------------------------------------ T5: history residue
def residue_control():
    rows = []
    for name, turn in (('quarter turn', (F(0), F(1))), ('rational turn (3 + 4R)/5', (F(3, 5), F(4, 5)))):
        vals = []
        for Q in range(1, 9):
            w = cpow(turn, Q)
            res = ((w[0]-1)**2+w[1]**2)/Q
            if res != 2*(1-w[0])/Q:
                raise ValueError('history residue is not 2 (1 - E) / Q')
            vals.append(res)
        rows.append(dict(turn=name, residues=[str(v) for v in vals], closes_at=[q+1 for q, v in enumerate(vals) if v == 0]))
    if rows[0]['closes_at'] != [4, 8] or rows[1]['closes_at']:
        raise ValueError('closure pattern failed')
    return rows


def run():
    pairs = []
    for c, s in ((F(3, 5), F(4, 5)), (F(4, 5), F(3, 5)), (F(5, 13), F(12, 13)), (F(1), F(0)), (F(0), F(1))):
        curv, rest = weyl_control(c, s), uniform_control(c, s)
        if curv != rest or 1-curv/2 != c:
            raise ValueError('mass^2 is not the clock curvature')
        pairs.append(dict(tick_phase=[str(c), str(s)], clock_curvature=str(curv), rest_defect=str(rest),
                          cone_speed=str(c)))
    return dict(mass_is_clock_curvature=pairs, tower=tower_control(), varying_coin=varying_coin_control(),
                count=count_control(), history=residue_control())


if __name__ == '__main__':
    res = run()
    json.dump(res, open('CL1_RESULT.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, v if k != 'history' else [(r['turn'], r['closes_at'], r['residues'][:4]) for r in v])
