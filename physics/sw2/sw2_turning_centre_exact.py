"""SW2: the turning centre exactly.  The centre at rest, displaced by iota*a along its axis.

Frame (coordinates t, r, x = cos theta, phi; flat space in the oblate layout  X + iY = sqrt(r^2+a^2) sin(theta) e^{i phi}, Z = r x):
    e^0 = dt ,   e^1 = (rho/s) dr + b nu ,   e^2 = rho dtheta ,   e^3 = s sin(theta) dphi ,
    nu = dt - a sin^2(theta) dphi ,   rho^2 = r^2 + a^2 x^2 ,   s^2 = r^2 + a^2 ,   b = K(r)/rho ,   K^2 = P(r) free.
With b = 0 this is flat space; with a = 0 it is the frame of MO1 with memory m = P/r^2.
The law is CV1's: zero contracted curvature of the frame (the connection from the order defects alone).
The square roots rho, s, sin(theta), K, n are carried as symbols and reduced by their squares; everything is exact.
Symbolic (sympy).  Python 3.12.  Takes several minutes."""
import json, os, sys, time
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))

r, x, a, rs = sp.symbols('r x a r_s', positive=True)
rho, s, sig, K, n = sp.symbols('rho s sigma K n', positive=True)
P, P1, P2, P3 = sp.symbols('P P1 P2 P3')
SQ = [(rho, r**2 + a**2*x**2), (s, r**2 + a**2), (sig, 1 - x**2), (K, P), (n, r**2 + a**2*x**2 - P)]
ETA = sp.diag(1, -1, -1, -1)


def red(e, squares=None):
    squares = SQ if squares is None else squares
    num, den = sp.fraction(sp.cancel(e))
    def lower(p):
        p = sp.expand(p)
        for sym, sq in squares:
            pp = sp.Poly(p, sym)
            p = sp.expand(sum(co*sq**(k//2)*sym**(k % 2) for (k,), co in pp.terms()))
        return p
    return sp.cancel(lower(num)/lower(den))


def Dr(f):
    return (sp.diff(f, r) + sp.diff(f, rho)*r/rho + sp.diff(f, s)*r/s + sp.diff(f, K)*P1/(2*K) + sp.diff(f, n)*(2*r - P1)/(2*n)
            + sp.diff(f, P)*P1 + sp.diff(f, P1)*P2 + sp.diff(f, P2)*P3)


def Dx(f):
    return sp.diff(f, x) + sp.diff(f, rho)*a**2*x/rho - sp.diff(f, sig)*x/sig + sp.diff(f, n)*a**2*x/n


def frame(b):
    """rows of E: the frame directions; rows of TH: the coframe (columns t, r, x, phi)."""
    E = sp.Matrix([[1, -b*s/rho, 0, 0], [0, s/rho, 0, 0], [0, 0, -sig/rho, 0], [0, a*sig*b/rho, 0, 1/(s*sig)]])
    TH = sp.Matrix([[1, 0, 0, 0], [b, rho/s, 0, -a*sig**2*b], [0, 0, -rho/sig, 0], [0, 0, 0, s*sig]])
    return E, TH


def ap(v, f):
    return v[1]*Dr(f) + v[2]*Dx(f)


def order_defect(E, TH):
    c = [[[sp.Integer(0)]*4 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            comm = [ap(E.row(i), E[j, m]) - ap(E.row(j), E[i, m]) for m in range(4)]
            for k in range(4):
                c[i][j][k] = red(sum(comm[m]*TH[k, m] for m in range(4)))
                c[j][i][k] = -c[i][j][k]
    return c


def connection(c):
    """CV1's formula: w[i][j][k] = w^i_{jk}, from the order defects alone; keeps eta, no torsion."""
    low = lambda i, j, k: ETA[k, k]*c[i][j][k]
    return [[[red(ETA[i, i]*sp.Rational(1, 2)*(low(k, j, i) - low(j, i, k) + low(i, k, j))) for k in range(4)] for j in range(4)] for i in range(4)]


def riemann(E, c, w, i, j, k, l):
    """CV1's formula: curvature = order defect of the connected directions."""
    v = ap(E.row(k), w[i][j][l]) - ap(E.row(l), w[i][j][k])
    v += sum(w[i][e][k]*w[e][j][l] - w[i][e][l]*w[e][j][k] for e in range(4))
    v -= sum(c[k][l][e]*w[i][j][e] for e in range(4))
    return v


def contracted(E, c, w):
    Ric = sp.zeros(4)
    for j in range(4):
        for l in range(j, 4):
            Ric[j, l] = Ric[l, j] = red(sum(riemann(E, c, w, i, j, i, l) for i in range(4)))
    return Ric


def flat_space_identities():
    out = {}
    z = lambda e: sp.simplify(e) == 0
    X, Y, Z = sp.symbols('X Y Z', real=True)
    aa = sp.Symbol('a', positive=True)
    R2 = X**2 + Y**2 + (Z - sp.I*aa)**2
    rr, xx, phi = sp.symbols('r x phi', real=True)
    sub = {X: sp.sqrt(rr**2 + aa**2)*sp.sqrt(1 - xx**2)*sp.cos(phi), Y: sp.sqrt(rr**2 + aa**2)*sp.sqrt(1 - xx**2)*sp.sin(phi), Z: rr*xx}
    # F1: the distance from the displaced centre is r - iota a x
    out['F1_complex_distance'] = z(sp.expand(R2.subs(sub)) - sp.expand((rr - sp.I*aa*xx)**2))
    # F2: 1/R is harmonic on flat space (MC1's least-cost law in three dimensions, continued)
    inv = 1/sp.sqrt(R2)
    out['F2_harmonic'] = z(sum(sp.diff(inv, v, 2) for v in (X, Y, Z)))
    # F3: (grad 1/R).(grad 1/R) = 1/R^4
    out['F3_square_of_gradient'] = z(sum(sp.diff(inv, v)**2 for v in (X, Y, Z)) - 1/R2**2)
    # F4: every multipole is (iota a)^l : 1/R = sum (iota a)^l P_l(cos) / d^(l+1)   (checked to l = 5 on the axis and off it)
    d, cth, eps = sp.symbols('d c epsilon', positive=True)
    gen = 1/sp.sqrt(1 - 2*cth*(sp.I*eps) + (sp.I*eps)**2)                     # d/R with eps = a/d
    ser = sp.series(gen, eps, 0, 6).removeO()
    out['F4_multipoles'] = z(sp.expand(ser - sum((sp.I*eps)**l*sp.legendre(l, cth) for l in range(6))))
    # the same generating function is 1/R: R^2 = d^2 - 2 iota a d c - a^2
    out['F4_generating'] = z(sp.expand(R2.subs({X: d*sp.sqrt(1 - cth**2), Y: 0, Z: d*cth}) - d**2*(1 - 2*cth*sp.I*aa/d + (sp.I*aa/d)**2)))
    return out


def first_order_is_sw1():
    """the exact count form to first order in a, after the relabelling phi -> phi + a F(r), F' = beta/r^2, is SW1's with W = r_s a / r^3."""
    aa, rr, th = sp.symbols('a r theta', positive=True)
    beta = sp.sqrt(rs/rr)
    dt, dr, dth, dph = sp.symbols('dt dr dtheta dphi')
    rho_ = sp.sqrt(rr**2 + aa**2*sp.cos(th)**2); s_ = sp.sqrt(rr**2 + aa**2)
    b_ = sp.sqrt(rs*rr)/rho_
    dphi_old = dph + aa*beta/rr**2*dr
    e1 = rho_/s_*dr + b_*(dt - aa*sp.sin(th)**2*dphi_old)
    form = dt**2 - e1**2 - rho_**2*dth**2 - s_**2*sp.sin(th)**2*dphi_old**2
    first = sp.series(sp.expand(form), aa, 0, 2).removeO()
    W = rs*aa/rr**3
    sw1 = dt**2 - (dr + beta*dt)**2 - rr**2*dth**2 - rr**2*sp.sin(th)**2*(dph - W*dt)**2      # SW1: e_0 = d_t - beta d_r + W d_phi
    sw1_first = sp.series(sp.expand(sw1), aa, 0, 2).removeO()
    return sp.simplify(sp.expand(first - sw1_first)) == 0


def wide_family_consequences():
    """reads sw2_wide_family.json (fall function free in r and x, written by sw2_wide_family.py) and the r-only result."""
    out = {}
    wide = {k: sp.sympify(v) for k, v in json.load(open(os.path.join(HERE, 'sw2_wide_family.json'))).items()}
    narrow = {k: sp.sympify(v) for k, v in json.load(open(os.path.join(HERE, 'sw2_contracted_curvature.json'))).items()}
    Pr, Px, Prr, Prx, Pxx = sp.symbols('Pr Px Prr Prx Pxx')
    rho2 = r**2 + a**2*x**2
    sq = [(rho, rho2), (s, r**2 + a**2), (sig, 1 - x**2), (K, P)]
    rd = lambda e: red(e, sq)
    # W1: radial + theta-theta components:  P_rr / (2 rho^2)
    out['W1_P_rr'] = rd(wide['11'] + wide['22'] - Prr/(2*rho2)) == 0
    # W2: a combination of the (time, radial) and (radial, round) components:  4 P P_x x rho^2
    A01 = wide['01']*8*K*P*rho*rho2**2
    A13 = -wide['13']*8*K*P*rho*s*rho2**2/(a*sig)
    out['W2_P_x'] = rd(A13 - A01 - 4*P*Px*x*rho2) == 0
    # W3: with P_x = 0 the wide family is the r-only family
    back = {Px: 0, Pxx: 0, Prx: 0, Pr: P1, Prr: P2}
    out['W3_reduces'] = all(rd(wide[k].subs(back) - narrow[k]) == 0 for k in wide)
    return out


def run():
    out, t0 = {}, time.time()
    out.update(flat_space_identities())
    out['L1_first_order_is_SW1'] = first_order_is_sw1()
    # the flat frame (b = 0): every curvature component vanishes
    E0, TH0 = frame(sp.Integer(0))
    c0 = order_defect(E0, TH0); w0 = connection(c0)
    out['G0_flat_layout'] = all(red(riemann(E0, c0, w0, i, j, k, l)) == 0 for i in range(4) for j in range(4) for k in range(4) for l in range(k + 1, 4))
    # the falling frame with a free fall function
    b = K/rho
    E, TH = frame(b)
    out['G0_frame_and_coframe'] = (E*TH.T).applyfunc(red) == sp.eye(4)
    c = order_defect(E, TH); w = connection(c)
    # G1: one time and free fall for every P
    out['G1_free_fall'] = all(w[i][0][0] == 0 for i in range(4))
    Ric = contracted(E, c, w)
    json.dump({'%d%d' % (j, l): sp.srepr(Ric[j, l]) for j in range(4) for l in range(j, 4)}, open(os.path.join(HERE, 'sw2_contracted_curvature.json'), 'w'), indent=0)
    rho2 = r**2 + a**2*x**2
    # G2: the contracted curvature is linear in P = (memory) x rho^2, with these components
    lin = lambda e: sp.expand(sp.numer(sp.together(e))).as_poly(P, P1, P2).total_degree() <= 1
    out['G2_linear_in_memory'] = all(lin(Ric[j, l]) for j in range(4) for l in range(4) if Ric[j, l] != 0)
    out['G2_theta_theta'] = red(Ric[2, 2] - (r*P1 - P)/rho2**2) == 0
    out['G2_radial'] = red(Ric[1, 1] - (2*P - 2*r*P1 + rho2*P2)/(2*rho2**2)) == 0
    out['G2_zero_components'] = all(Ric[j, l] == 0 for (j, l) in ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3)))
    # G3: the law (all components zero)  <=>  r P' = P , i.e. P = r_s r : memory = Re( r_s / (r - iota a x) )
    vac = {P: rs*r, P1: rs, P2: 0, P3: 0}
    SQV = [(sym, sq.subs(vac)) for sym, sq in SQ]
    redv = lambda e: red(sp.sympify(e).subs(vac), SQV)
    out['G3_law_holds_exactly'] = all(redv(Ric[j, l]) == 0 for j in range(4) for l in range(4))
    out['G3_memory_is_real_part'] = sp.simplify(sp.re(rs/(r - sp.I*a*x)) - (rs*r/rho2)) == 0
    # a = 0: CV1's own components, with memory m = P/r^2:  ((r m)''/2r , (r m)'/r^2), signs as in EG1's pattern
    m = sp.Function('m')(r)
    cv_sub = {P: r**2*m, P1: sp.diff(r**2*m, r), P2: sp.diff(r**2*m, r, 2)}
    at0 = Ric.subs(a, 0).subs(cv_sub).applyfunc(sp.simplify)
    rm1, rm2 = sp.diff(r*m, r), sp.diff(r*m, r, 2)
    out['G4_reduces_to_CV1'] = (at0 - sp.diag(-rm2/(2*r), rm2/(2*r), rm1/r**2, rm1/r**2)).applyfunc(sp.simplify) == sp.zeros(4)
    # G5: the shell r is light-like where r^2 + a^2 = r_s r ; none when a > r_s/2
    grr = red(E[0, 1]**2 - E[1, 1]**2 - E[3, 1]**2)
    out['G5_lightlike_shell'] = red(grr - (P - r**2 - a**2)/rho2) == 0
    # H: a held reading u = (e_0 + b e_1)/N , N^2 = 1 - b^2 = 1 - Re Phi : acceleration and turn of the held frame
    N = n/rho
    u = [1/N, b/N, 0, 0]
    Du = [[red(ap(E.row(cc), u[aa]) + sum(w[aa][bb][cc]*u[bb] for bb in range(4))) for cc in range(4)] for aa in range(4)]
    acc = [red(ETA[aa, aa]*sum(u[cc]*Du[aa][cc] for cc in range(4))) for aa in range(4)]                         # a_a
    Om = [red(sp.Rational(1, 2)*sum(sp.LeviCivita(aa, bb, cc, dd)*u[bb]*ETA[cc, cc]*Du[dd][cc]
                                    for bb in range(4) for cc in range(4) for dd in range(4))) for aa in range(4)]   # Omega_a
    re_phi, im_phi = rs*r/rho**2, rs*a*x/rho**2
    out['H0_real_and_imaginary_parts'] = sp.simplify(rs/(r - sp.I*a*x) - (rs*r + sp.I*rs*a*x)/rho2) == 0
    gre = [red(ap(E.row(i), re_phi)) for i in range(4)]
    gim = [red(ap(E.row(i), im_phi)) for i in range(4)]
    out['H1_acceleration_is_gradient_of_real_part'] = all(redv(acc[i]*2*N**2 - gre[i]) == 0 for i in range(4))
    out['H2_turn_is_gradient_of_imaginary_part'] = all(redv(Om[i]*2*N**2 - gim[i]) == 0 for i in range(4))
    out['H3_clock_factor'] = redv(N**2 - (1 - re_phi)) == 0
    if os.path.exists(os.path.join(HERE, 'sw2_wide_family.json')):
        out.update(wide_family_consequences())
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, seconds=round(time.time() - t0)), open(os.path.join(HERE, 'SW2_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
