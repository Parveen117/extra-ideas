"""FD1: one function for every memory: det(1 + zX) = sum z^n e_n(X).

MM1: every memory of the line is e_2.  SH1: the N-fold memory of a record is e_N.  The function that carries all
e_n at once is the determinant det(1 + zX) on the cut carrier.  This stage builds it exactly and reads from it:
the split S = R + D at second order, closed circuits, the split across a cut, counts through a cut, rates, and an
infinite ladder with a declared tail.  sympy and exact rationals.  Python 3.12."""
import itertools, json, math, os, random
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def series_coeffs(expr, z, order):
    return [sp.simplify(sp.diff(expr, z, n).subs(z, 0)/sp.factorial(n)) for n in range(order + 1)]


def closed_walks(K, n):
    """closed n-step walks on the graph with adjacency K, as tuples of vertices."""
    N = K.shape[0]
    out = []
    for w in itertools.product(range(N), repeat=n):
        if all(K[w[i], w[(i + 1) % n]] for i in range(n)):
            out.append(w)
    return out


def primitive_circuits(K, n):
    """closed n-step walks up to rotation that are not repeats of a shorter one."""
    seen, count = set(), 0
    for w in closed_walks(K, n):
        if w in seen:
            continue
        rots = {w[i:] + w[:i] for i in range(n)}
        seen |= rots
        if len(rots) == n:
            count += 1
    return count


def run():
    out, num = {}, {}
    z = sp.symbols('z')
    zz = lambda e: sp.simplify(e) == 0

    # F1: the function and its first coefficients
    X = sp.Matrix(3, 3, sp.symbols('x0:9'))
    co = sp.Poly(sp.expand((sp.eye(3) + z*X).det()), z).all_coeffs()[::-1]
    e1, e2, e3 = X.trace(), (X.trace()**2 - (X*X).trace())/2, X.det()
    out['F1_coefficients_are_e_n'] = zz(co[1] - e1) and zz(co[2] - e2) and zz(co[3] - e3)
    n_, r1, r2, r3 = sp.symbols('n r1 r2 r3', real=True)
    K_ = sp.diag(1, -1); R_ = sp.Matrix([[0, -1], [1, 0]]); C = [K_, R_*K_, sp.I*R_]
    rho = (n_*sp.eye(2) + r1*C[0] + r2*C[1] + r3*C[2])/2
    out['F1_reading_tensor'] = zz((sp.eye(2) + z*rho).det() - (1 + z*n_ + z**2*(n_**2 - r1**2 - r2**2 - r3**2)/4))     # IN1
    # the reciprocal function carries h_n; second order: h_2 = R + D, e_2 = R - D, R = (tr)^2/2, D = tr(X^2)/2
    X2 = sp.Matrix(2, 2, sp.symbols('y0:4'))
    hser = series_coeffs(1/(sp.eye(2) - z*X2).det(), z, 3)
    eser = sp.Poly(sp.expand((sp.eye(2) + z*X2).det()), z).all_coeffs()[::-1] + [0]
    Rr, Dd = X2.trace()**2/2, (X2*X2).trace()/2
    out['F1_second_order_is_S_R_F'] = zz(hser[2] - (Rr + Dd)) and zz(eser[2] - (Rr - Dd)) and zz(series_coeffs(sp.exp(z*X2.trace()), z, 2)[2] - Rr)
    out['F1_the_two_functions_are_reciprocal'] = zz(sp.series(sp.expand((sp.eye(2) + z*X2).det())*(1/(sp.eye(2) + z*X2).det()), z, 0, 4).removeO() - 1) and all(zz(sum(eser[j]*(-1)**(m - j)*hser[m - j] for j in range(m + 1))) for m in (1, 2, 3))
    # gravity's vacuum law (SE1, DO1) is the vanishing of the z^2 coefficient: stretch (rho, 1, ..., 1) in d cuts
    rho_s = sp.symbols('rho')
    ok = True
    for d in (3, 4, 5, 6):
        St = sp.diag(rho_s, *([1]*(d - 1)))
        c2 = sp.Poly(sp.expand((sp.eye(d) + z*St).det()), z).all_coeffs()[::-1][2]
        ok &= sp.solve(c2, rho_s) == [-sp.Rational(d - 2, 2)]
    out['F1_vacuum_law_is_second_coefficient_zero'] = ok

    # F2: closed circuits.  -log det(1 - zK) = sum z^n tr K^n / n ; for a graph: an Euler product over primitive circuits
    Kg = sp.Matrix(2, 2, sp.symbols('k0:4'))
    lhs = sp.series(-sp.log((sp.eye(2) - z*Kg).det()), z, 0, 6).removeO()
    out['F2_log_det_is_the_sum_over_circuits'] = zz(sp.expand(lhs - sum(z**m*(Kg**m).trace()/m for m in range(1, 6))))
    okc, oke, prim = True, True, {}
    for name, A in (('two marks', sp.Matrix([[1, 1], [1, 0]])), ('triangle', sp.Matrix([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))):
        M = 9 if A.shape[0] == 2 else 7
        P = {m: primitive_circuits(A, m) for m in range(1, M + 1)}
        okc &= all((A**m).trace() == len(closed_walks(A, m)) for m in range(1, M + 1))
        prod = sp.Integer(1)
        for m, p in P.items():
            prod *= (1 - z**m)**p
        diff = sp.series(sp.expand(prod) - (sp.eye(A.shape[0]) - z*A).det(), z, 0, M + 1).removeO()
        oke &= sp.expand(diff) == 0
        prim[name] = list(P.values())
    out['F2_trace_counts_closed_walks'] = okc
    out['F2_euler_product_over_primitive_circuits'] = oke
    num['primitive_circuits'] = prim

    # F3: across a cut P + Q = 1:  det = det(seen block) x det(memory block dressed by the exchange)
    random.seed(7)
    A6 = sp.Matrix(6, 6, lambda i, j: sp.Rational(random.randint(-5, 5), random.randint(1, 4)))
    a, b, c_, d_ = A6[:3, :3], A6[:3, 3:], A6[3:, :3], A6[3:, 3:]
    lhs = (sp.eye(6) - z*A6).det()
    rhs = (sp.eye(3) - z*a).det()*(sp.eye(3) - z*d_ - z**2*c_*(sp.eye(3) - z*a).inv()*b).det()
    out['F3_split_across_a_cut'] = zz(sp.cancel(lhs - rhs))
    # a Gram form S = Z*Z (T24): e_k(Z*Z) = e_k(ZZ*); second order: e2 = e2(seen) + e2(lost) + seen x lost - exchange^2
    Z = sp.Matrix(5, 3, lambda i, j: sp.Rational(random.randint(-4, 4), random.randint(1, 3)))
    G1, G2 = Z.T*Z, Z*Z.T
    e2f = lambda M: (M.trace()**2 - (M*M).trace())/2
    As, Bs, Cs = G2[:2, :2], G2[:2, 2:], G2[2:, 2:]
    out['F3_gram_second_order'] = e2f(G1) == e2f(G2) == e2f(As) + e2f(Cs) + As.trace()*Cs.trace() - (Bs*Bs.T).trace()
    v = sp.Matrix(4, 1, sp.symbols('v0:4', real=True))
    Gv = v*v.T
    out['F3_single_reading_seen_times_lost_is_exchange_squared'] = zz(Gv[:2, :2].trace()*Gv[2:, 2:].trace() - (Gv[:2, 2:]*Gv[:2, 2:].T).trace()) and zz(e2f(Gv))

    # F4: counts through a cut.  N independent readings (columns of A) as one exclusive record (SH1-S4: squared wedge).
    def count_law(A, p):
        m, N = A.shape
        G = (A.T*A).det()
        prob = [sp.Integer(0)]*(N + 1)
        for S in itertools.combinations(range(m), N):
            prob[sum(1 for i in S if i < p)] += A[list(S), :].det()**2
        return [q/G for q in prob]
    s = sp.symbols('s')
    okd = True
    for (m, N, p) in ((4, 2, 2), (5, 3, 2), (6, 3, 3)):
        A = sp.Matrix(m, N, lambda i, j: sp.Rational(random.randint(-3, 3), random.randint(1, 2)))
        while (A.T*A).det() == 0:
            A = sp.Matrix(m, N, lambda i, j: sp.Rational(random.randint(-3, 3), random.randint(1, 2)))
        prob = count_law(A, p)
        Pm = sp.diag(*([1]*p + [0]*(m - p)))
        KP = (A.T*A).inv()*A.T*Pm*A
        gen = sp.expand((sp.eye(N) + (s - 1)*KP).det())
        okd &= sum(prob) == 1 and sp.expand(sum(q*s**k for k, q in enumerate(prob)) - gen) == 0
        okd &= sum(k*q for k, q in enumerate(prob)) == KP.trace()
        okd &= sum(k*k*q for k, q in enumerate(prob)) - KP.trace()**2 == (KP*(sp.eye(N) - KP)).trace()
    out['F4_count_through_a_cut_is_a_determinant'] = okd
    # on the diagonal (every reading seen = lost) the count is that of fair coins
    N = 3
    Ad = sp.Matrix(2*N, N, lambda i, j: 1 if (i == j or i == j + N) else 0)
    out['F4_diagonal_cut_counts_fair_coins'] = count_law(Ad, N) == [sp.Rational(math.comb(N, k), 2**N) for k in range(N + 1)]
    # the three kinds of count in one form: G = det(1 - eta w N)^(-eta), w = s - 1; second cumulant n (1 + eta n)
    w, nn, lam = sp.symbols('w n lambda', positive=True)
    okk = True
    for eta in (1, -1):
        Gf = (1 - eta*(sp.exp(lam) - 1)*nn)**(-eta)
        okk &= zz(sp.diff(sp.log(Gf), lam, 2).subs(lam, 0) - nn*(1 + eta*nn)) and zz(sp.diff(sp.log(Gf), lam).subs(lam, 0) - nn)
    out['F4_three_counts_one_form'] = okk and zz(sp.limit((1 - sp.Symbol('e')*w*nn)**(-1/sp.Symbol('e')), sp.Symbol('e'), 0) - sp.exp(w*nn))

    # F5: rates are ratios of zeros.  The two-valued chain of DU1
    a_, b_, k = sp.symbols('a b k', positive=True)
    C2 = R_*K_; H = (K_ + C2)/sp.sqrt(2)
    B = a_*sp.eye(2) + b_*C2
    fd = sp.factor((sp.eye(2) - z*B).det())
    out['F5_zeros_are_inverse_readings'] = zz(fd - (1 - z*(a_ + b_))*(1 - z*(a_ - b_)))
    out['F5_rate_is_log_ratio_of_zeros'] = zz(sp.simplify((sp.log((a_ + b_)/(a_ - b_)).subs({a_: sp.exp(k), b_: sp.exp(-k)}) + sp.log(sp.tanh(k))).rewrite(sp.exp)))
    out['F5_diagonal_turn_keeps_the_function'] = zz((sp.eye(2) - z*H*B*H).det() - (sp.eye(2) - z*B).det())
    out['F5_doubling_the_cell'] = zz((sp.eye(2) - z**2*B*B).det() - (sp.eye(2) - z*B).det()*(sp.eye(2) + z*B).det())

    # F6: an infinite ladder with a declared tail (shape of theorum/28 sections 4 and 9): readings q^i, i = 0, 1, 2, ...
    q = sp.symbols('q')
    M = 6
    polyM = sp.Poly(sp.expand(sp.prod([1 + z*q**i for i in range(M)])), z).all_coeffs()[::-1]
    def qbin(M_, n):
        return sp.cancel(sp.prod([(1 - q**(M_ - n + j)) for j in range(1, n + 1)])/sp.prod([(1 - q**j) for j in range(1, n + 1)]))
    out['F6_finite_ladder_exact'] = all(zz(polyM[n] - q**(n*(n - 1)//2)*qbin(M, n)) for n in range(M + 1))
    okt = True
    for qv in (Fr(1, 2), Fr(1, 3), Fr(9, 10)):
        for Mv in (8, 14, 40):
            for n in range(1, 6):
                einf = qv**(n*(n - 1)//2)/math.prod(1 - qv**j for j in range(1, n + 1))
                eM = einf*math.prod(1 - qv**(Mv - n + j) for j in range(1, n + 1))
                okt &= 0 <= einf - eM <= einf*qv**(Mv - n + 1)/(1 - qv)
    out['F6_declared_tail_bounds_every_coefficient'] = okt
    num['ladder_tail_example'] = dict(q=0.5, M=40, relative_tail_bound_for_e5=float(Fr(1, 2)**36/(1 - Fr(1, 2))))

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'FD1_RESULT.json'), 'w'), indent=1, default=float)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
    for k_, v in n.items(): print(k_, v)
