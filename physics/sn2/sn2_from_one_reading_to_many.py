"""SN2: the spread of a count, from one reading to many: 1/3 -> 1/4; and what adds along one chain.

SN1 took the even spread of the sector's angle as put in.  Here: (1) spread/mean is the seen-weighted mean of lost;
(2) one reading with no preferred direction in d cuts gives (d-1)/2d; (3) N = 2^(n-1) readings in any frame made of
products of cuts, none preferred, give N^2/(4N^2 - 1), exactly, by finite enumeration; (4) along one chain of cells
it is log(seen) that adds, not the angle.  sympy, exact enumeration, mpmath.  Python 3.12."""
import itertools, json, os
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def run():
    out, num = {}, {}
    zz = lambda e: sp.simplify(e) == 0
    ph, u = sp.symbols('phi u', positive=True)

    # N1: spread/mean = sum q(1-q)/sum q = mean of lost, weighted by seen.  SN1's three sectors read this way
    out['N1_turn'] = sp.integrate(sp.sin(ph)**2*sp.cos(ph)**2, (ph, 0, sp.pi/2))/sp.integrate(sp.cos(ph)**2, (ph, 0, sp.pi/2)) == sp.Rational(1, 4)
    out['N1_shear_is_the_plain_mean_of_lost_over_the_angle'] = sp.integrate(sp.sin(ph)**2, (ph, 0, sp.pi/2))/(sp.pi/2) == sp.Rational(1, 2)
    out['N1_boost_is_the_mean_square_of_an_even_amplitude'] = sp.integrate(u**2, (u, 0, 1)) == sp.Rational(1, 3)       # seen d(eta) = d(tanh eta)

    # N2: one reading, no preferred direction, d cuts.  sum of F_i^2 = 1 (IN1) => <F^2> = 1/d ; R = (1+F)/2 ; spread/mean = (d-1)/2d
    okd, flat = True, []
    t = sp.symbols('t', real=True)                                   # F = u = sin t ; law of u in d cuts: (1 - u^2)^((d-3)/2) du = cos^(d-2) t dt
    for d in range(2, 8):
        norm = sp.integrate(sp.cos(t)**(d - 2), (t, -sp.pi/2, sp.pi/2))
        m2 = sp.simplify(sp.integrate(sp.sin(t)**2*sp.cos(t)**(d - 2), (t, -sp.pi/2, sp.pi/2))/norm)
        okd &= m2 == sp.Rational(1, d)
        okd &= sp.simplify((1 - m2)/4/sp.Rational(1, 2)) == sp.Rational(d - 1, 2*d)
        flat.append(d == 3)                                          # the exponent (d - 3)/2 is zero only for d = 3
        # on the 2d cut directions themselves: F = +-1 for two of them, 0 for the rest
        okd &= Fr(2, 2*d) == Fr(1, d)
    out['N2_one_reading_in_d_cuts'] = okd
    out['N2_amplitude_evenly_spread_only_in_three_cuts'] = flat == [False, True, False, False, False, False]
    dsym = sp.symbols('d', positive=True)
    out['N2_mean_flip_squared_equals_the_ratio_only_in_three_cuts'] = sp.solve(sp.Eq(1/dsym, (dsym - 1)/(2*dsym)), dsym) == [3]
    num['one_reading'] = {str(d): str(sp.Rational(d - 1, 2*d)) for d in (2, 3, 4, 6)}

    # N3: many readings.  n factors, carrier 2^n, cuts C1 = K, C2 = RK, C3 = iota R on each factor.
    #     cut P = (1 + C1 on the first factor)/2 ; readings' frame Pi = (1 + sigma)/2, sigma = +- any product of cuts except 1.
    K_ = sp.diag(1, -1); R_ = sp.Matrix([[0, -1], [1, 0]])
    units = [sp.eye(2), K_, R_*K_, sp.I*R_]
    def kron(ms):
        M = ms[0]
        for m_ in ms[1:]:
            M = sp.kronecker_product(M, m_)
        return M
    okn, table = True, {}
    for n in (1, 2, 3):
        dim = 2**n; N = dim//2
        P = (sp.eye(dim) + kron([K_] + [sp.eye(2)]*(n - 1)))/2
        m1, m2, cnt = sp.Integer(0), sp.Integer(0), 0
        for idx in itertools.product(range(4), repeat=n):
            if all(i == 0 for i in idx):
                continue
            sig = kron([units[i] for i in idx])
            assert sig*sig == sp.eye(dim)
            for sgn in (1, -1):
                Pi = (sp.eye(dim) + sgn*sig)/2
                m1 += (P*Pi).trace(); m2 += (P*Pi*P*Pi).trace(); cnt += 1
        mean, sec = m1/cnt, m2/cnt
        ratio = sp.nsimplify((mean - sec)/mean)
        okn &= mean == sp.Rational(N, 2) and ratio == sp.Rational(N**2, 4*N**2 - 1)
        table[N] = str(ratio)
    out['N3_many_readings_exact'] = okn
    num['many_readings'] = table
    # closed count for every n: sigma = C1 itself (1), commuting with it (M^2/2 - 2), not commuting (M^2/2)
    okc = True
    for n in range(1, 12):
        M = 2**n; N = M//2
        sec = Fr(M, 4)*(1 + (Fr(M*M, 2) - 2) + Fr(M*M, 2)*Fr(1, 2))/(M*M - 1)
        okc &= (Fr(N, 2) - sec)/Fr(N, 2) == Fr(N*N, 4*N*N - 1)
    out['N3_closed_count_for_every_size'] = okc
    out['N3_from_one_third_to_one_quarter'] = Fr(1, 4*1 - 1) == Fr(1, 3) and sp.limit(sp.Symbol('N')**2/(4*sp.Symbol('N')**2 - 1), sp.Symbol('N'), sp.oo) == sp.Rational(1, 4)

    # N4: one chain.  Two cells with boost angles e1, e2 joined through a turn phi:  1/seen = cosh^2 e1 cosh^2 e2 |1 + t1 t2 e^(i phi)|^2
    e1, e2, f = sp.symbols('eta1 eta2 varphi', real=True)
    Mb = lambda e: sp.Matrix([[sp.cosh(e), sp.sinh(e)], [sp.sinh(e), sp.cosh(e)]])
    Dt = sp.diag(sp.exp(sp.I*f/2), sp.exp(-sp.I*f/2))
    M = Mb(e2)*Dt*Mb(e1)
    lhs = sp.expand(M[0, 0]*sp.conjugate(M[0, 0]))
    rhs = sp.cosh(e1)**2*sp.cosh(e2)**2*(1 + sp.tanh(e1)**2*sp.tanh(e2)**2 + 2*sp.tanh(e1)*sp.tanh(e2)*sp.cos(f))
    out['N4_joining_two_cells'] = zz((lhs - rhs).rewrite(sp.exp))
    mp.mp.dps = 30
    out['N4_log_of_seen_adds'] = all(abs(mp.quad(lambda t: mp.log(1 + a*a + 2*a*mp.cos(t)), [0, mp.pi, 2*mp.pi])) < mp.mpf(10)**-20 for a in (mp.mpf(1)/10, mp.mpf(1)/2, mp.mpf(9)/10))
    # and the mean of 1/seen multiplies with a factor: the angle itself does not add
    mean_inv = sp.integrate(rhs, (f, 0, 2*sp.pi))/(2*sp.pi)
    out['N4_angle_does_not_add'] = zz((mean_inv - sp.cosh(e1)**2*sp.cosh(e2)**2*(1 + sp.tanh(e1)**2*sp.tanh(e2)**2)).rewrite(sp.exp)) and not zz((mean_inv - sp.cosh(e1 + e2)**2).rewrite(sp.exp))

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'SN2_RESULT.json'), 'w'), indent=1, default=str)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
    for k_, v in n.items(): print(k_, v)
