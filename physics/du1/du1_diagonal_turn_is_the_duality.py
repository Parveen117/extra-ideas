"""DU1: the diagonal turn exchanges the two readings of a two-valued record; the rate is the dual coupling.

EMK block: cuts C1 = K (which mark), C2 = RK (the flip).  The turn through their diagonal, H = (C1 + C2)/sqrt2,
exchanges them (DO1-P1).  A chain of two-valued marks with coupling k has block a + b C2, a = e^k, b = e^-k.
MG1 asked for a native exponential law; YM-21 gave the geometric law r^(cells).  sympy, exact enumeration, numpy for
finite strips (evidence only).  Python 3.12."""
import itertools, json, math, os
from collections import Counter
import numpy as np
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def even_subgraphs(m, n):
    """m x n grid of marks, free edges.  Returns (count of even subgraphs by size, count of cuts of the dual by size)."""
    V = [(i, j) for i in range(m) for j in range(n)]
    idx = {v: a for a, v in enumerate(V)}
    OUT = (m - 1)*(n - 1)
    face = lambda i, j: i*(n - 1) + j if 0 <= i < m - 1 and 0 <= j < n - 1 else OUT
    edges, dual = [], []
    for i in range(m):
        for j in range(n - 1):                                   # horizontal edge (i,j)-(i,j+1): faces above and below
            edges.append((idx[i, j], idx[i, j + 1])); dual.append((face(i - 1, j), face(i, j)))
    for i in range(m - 1):
        for j in range(n):                                       # vertical edge (i,j)-(i+1,j): faces left and right
            edges.append((idx[i, j], idx[i + 1, j])); dual.append((face(i, j - 1), face(i, j)))
    E = len(edges)
    masks = [(1 << a) | (1 << b) for a, b in edges]
    high = Counter()
    for sub in range(1 << E):
        par, cnt, s, e = 0, 0, sub, 0
        while s:
            if s & 1:
                par ^= masks[e]; cnt += 1
            s >>= 1; e += 1
        if par == 0:
            high[cnt] += 1
    low = Counter()
    for sig in range(1 << (OUT + 1)):
        low[sum(((sig >> f) ^ (sig >> g)) & 1 for f, g in dual)] += 1
    for kk in low:
        low[kk] //= 2
    return high, low, E, len(V)


def strip_gap(k, L):
    """rate (log of the ratio of the two largest readings) of the transfer block of a ring of L marks, same coupling both ways."""
    conf = np.array(list(itertools.product((1, -1), repeat=L)))
    d = np.exp(0.5*k*np.sum(conf*np.roll(conf, 1, axis=1), axis=1))
    one = np.array([[math.exp(k), math.exp(-k)], [math.exp(-k), math.exp(k)]])
    V2 = one
    for _ in range(L - 1):
        V2 = np.kron(V2, one)
    ev = np.linalg.eigvalsh(d[:, None]*V2*d[None, :])
    return math.log(ev[-1]/ev[-2])


def run():
    out, num = {}, {}
    z = lambda e: sp.simplify(e) == 0
    a, b, k = sp.symbols('a b k', positive=True)

    # U1: the diagonal turn
    C1 = sp.diag(1, -1); Rm = sp.Matrix([[0, -1], [1, 0]]); C2 = Rm*C1
    H = (C1 + C2)/sp.sqrt(2)
    out['U1_cuts'] = Rm**2 == -sp.eye(2) and C1**2 == sp.eye(2) and C2**2 == sp.eye(2) and C1*C2 + C2*C1 == sp.zeros(2)
    out['U1_diagonal_turn_exchanges_the_cuts'] = sp.simplify(H*H) == sp.eye(2) and sp.simplify(H*C1*H) == C2 and sp.simplify(H*C2*H) == C1
    B = a*sp.eye(2) + b*C2
    out['U1_block_read_on_the_diagonal'] = sp.simplify(H*B*H) == sp.diag(a + b, a - b)

    # U2: the dual record: weights (a + b, a - b):  e^(-2 k*) = tanh k
    ks = -sp.log(sp.tanh(k))/2
    out['U2_dual_is_an_involution'] = z(sp.simplify((-sp.log(sp.tanh(ks))/2 - k).rewrite(sp.exp)))
    out['U2_product_of_sinh'] = z(sp.simplify((sp.sinh(2*k)*sp.sinh(2*ks)).rewrite(sp.exp)) - 1)
    kc = sp.log(1 + sp.sqrt(2))/2
    t8 = sp.sqrt(2) - 1                                             # e^(-2 kc)
    out['U2_self_dual_point'] = z((1/t8 - t8)/2 - 1) and z((1 - t8)/(1 + t8) - t8)      # sinh 2kc = 1 ; tanh kc = e^(-2kc)
    out['U2_block_points_along_the_sixteenth'] = z(sp.exp(-2*kc) - t8) and z(sp.tan(sp.pi/8) - t8) and sp.simplify(H*sp.Matrix([sp.cos(sp.pi/8), sp.sin(sp.pi/8)]) - sp.Matrix([sp.cos(sp.pi/8), sp.sin(sp.pi/8)])) == sp.zeros(2, 1)
    th2, ch2 = sp.simplify(sp.tanh(2*kc).rewrite(sp.exp)), sp.simplify(sp.cosh(2*kc).rewrite(sp.exp))
    out['U2_self_dual_is_seen_equals_lost'] = z(th2 - 1/sp.sqrt(2)) and z(th2**2 - sp.Rational(1, 2)) and z(1/ch2**2 - sp.Rational(1, 2))
    num['self_dual_coupling'] = float(kc); num['one_over_the_self_dual_coupling'] = float(1/kc)

    # U3: the chain.  ring of N marks, exact: <s0 sn> = (t^n + t^(N-n)) / (1 + t^N), t = tanh k = e^(-2 k*)
    t = sp.symbols('t', positive=True)
    N = 6
    Z, corr = 0, {n_: 0 for n_ in range(1, N)}
    for s in itertools.product((1, -1), repeat=N):
        w = sp.Integer(1)
        for i in range(N):
            w *= (1 + t*s[i]*s[(i + 1) % N])
        Z += w
        for n_ in corr:
            corr[n_] += s[0]*s[n_]*w
    Z = sp.expand(Z)
    out['U3_chain_decay_is_geometric'] = all(z(sp.expand(corr[n_])/Z - (t**n_ + t**(N - n_))/(1 + t**N)) for n_ in corr)
    out['U3_rate_is_twice_the_dual_coupling'] = z(sp.simplify((-sp.log(sp.tanh(k)) - 2*ks)))
    # doubling the cell: sum over the middle mark
    B2 = sp.expand((a*sp.eye(2) + b*C2)**2)
    a2, b2 = a**2 + b**2, 2*a*b
    out['U3_doubling_squares_the_ratio'] = B2 == a2*sp.eye(2) + b2*C2 and z((a2 - b2)/(a2 + b2) - ((a - b)/(a + b))**2)
    kp = sp.log(sp.cosh(2*k))/2                                    # e^(2k') = a2/b2 = cosh 2k
    out['U3_doubling_in_the_coupling'] = z(sp.simplify((sp.exp(2*kp) - (a2/b2).subs({a: sp.exp(k), b: sp.exp(-k)})).rewrite(sp.exp)))
    out['U3_doubling_doubles_the_dual_coupling'] = z(sp.simplify((sp.tanh(kp) - sp.tanh(k)**2).rewrite(sp.exp)))
    mp.mp.dps = 80
    out['U3_coupling_runs_as_a_log'] = all(0 < mp.log(mp.cosh(2*mp.mpf(v)))/2 - (v - mp.log(2)/2) < mp.e**(-4*mp.mpf(v)) for v in (1, 2, 5, 10, 20))

    # U4: the exponential law.  rate = 2 artanh(e^-2k):  2 e^-2k < rate < 2 e^-2k / (1 - e^-4k)
    ok = True
    for v in (mp.mpf(1)/2, 1, 3, 10, 22):
        rate = -mp.log(mp.tanh(v)); lo = 2*mp.e**(-2*v)
        ok &= lo < rate < lo/(1 - mp.e**(-4*v))
    out['U4_exponential_law'] = bool(ok)
    hbar, cc, GG, me, mp_ = 1.054571817e-34, 299792458.0, 6.67430e-11, 9.1093837015e-31, 1.67262192369e-27
    mP = math.sqrt(hbar*cc/GG)
    num['coupling_for_the_proton_count'] = 0.5*math.log(2/(mp_/mP)); num['coupling_for_the_electron_count'] = 0.5*math.log(2/(me/mP))
    out['U4_a_number_near_22_not_a_large_number'] = 22 < num['coupling_for_the_proton_count'] < 23
    # contrast: the chain of turns (YM-21): ratio I2/I1; doubling halves the coupling (a power, no log)
    r = lambda v: mp.besseli(2, v)/mp.besseli(1, v)
    v0 = mp.mpf(1000); v1 = mp.findroot(lambda v: r(v) - r(v0)**2, v0/2)
    num['turn_chain_doubling_factor'] = float(v1/v0)
    out['U4_turn_chain_runs_as_a_power'] = abs(v1/v0 - mp.mpf(1)/2) < mp.mpf(1)/100

    # U5: the plane.  even subgraphs of a grid = cuts of its dual, size by size (exact enumeration)
    okp, sizes = True, {}
    for (m_, n_) in [(2, 2), (2, 3), (3, 3), (3, 4)]:
        high, low, E, Nv = even_subgraphs(m_, n_)
        okp &= high == low
        sizes[f'{m_}x{n_}'] = dict(edges=E, even_subgraphs=sum(high.values()))
    out['U5_plane_duality_exact'] = okp
    num['grids'] = sizes

    # U6 (evidence, not proof): rings of width L with the same coupling both ways.  Known result for the infinite plane: rate = 2 |k* - k|
    kv = 0.3; ksv = -0.5*math.log(math.tanh(kv))
    g_off = {L: strip_gap(kv, L) for L in (4, 6, 8, 10)}
    kcv = float(kc)
    g_on = {L: L*strip_gap(kcv, L) for L in (4, 6, 8, 10)}
    num['strip_rate_off_the_diagonal'] = g_off; num['known_plane_rate'] = 2*(ksv - kv)
    num['width_times_rate_on_the_diagonal'] = g_on; num['quarter_of_pi'] = math.pi/4
    out['U6_strips_approach_twice_the_distance_to_the_dual'] = abs(g_off[10] - 2*(ksv - kv)) < 1e-2 and abs(g_off[10] - 2*(ksv - kv)) < abs(g_off[4] - 2*(ksv - kv))
    out['U6_on_the_diagonal_rate_falls_as_one_over_width'] = abs(g_on[10] - math.pi/4) < 0.01 and abs(g_on[10] - math.pi/4) < abs(g_on[4] - math.pi/4)

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'DU1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
    for k_, v in n.items(): print(k_, v)
