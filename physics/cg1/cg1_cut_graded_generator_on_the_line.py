"""CG1: the cut-graded generator (RKF theorum/41) on the physics line.

theorum/41: for a cut J (J^2 = 1) and a generator G:  G_e = (G + JGJ)/2 keeps the two sheets, G_o = (G - JGJ)/2 carries
between them;  for cut-odd G:  J U_t J = U_-t,  E_t^2 - O_t^2 = 1,  (U_t + U_-t)(U_t - U_-t) = U_2t - U_-2t;  for any G:
log(J U_t J U_t) = 2t G_e + t^2 [G_e, G_o] + O(t^3).
Here these are checked on the cut-complex block (only J^2 = 1 is used; no dagger), and the line's stages of today are read
as instances.  sympy.  Python 3.12."""
import json, os, random
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def exp_series(G, t, order):
    n = G.shape[0]
    out, term = sp.zeros(n), sp.eye(n)
    for k in range(order + 1):
        out += term*t**k/sp.factorial(k)
        term = term*G
    return out


def coeffs(M, t, order):
    return [M.applyfunc(lambda e: sp.expand(e).coeff(t, k)) for k in range(order + 1)]


def run():
    out, num = {}, {}
    Z = lambda M: M.applyfunc(sp.simplify) == sp.zeros(*M.shape)
    zz = lambda e: sp.simplify(e) == 0
    t = sp.symbols('t')
    K_ = sp.diag(1, -1); R_ = sp.Matrix([[0, -1], [1, 0]]); I2 = sp.eye(2)
    C1, C2, C3 = K_, R_*K_, sp.I*R_
    out['G0_cuts'] = C1*C2 == sp.I*C3 and C1*C2*C3 == sp.I*I2 and all(c*c == I2 for c in (C1, C2, C3))

    # G1: theorum/41 (i) and (iv) on the block, cut J = C1, general G
    g0, g1, g2, g3 = sp.symbols('g0:4')
    G = g0*I2 + g1*C1 + g2*C2 + g3*C3
    Ge, Go = (G + C1*G*C1)/2, (G - C1*G*C1)/2
    out['G1_grading'] = Z(Ge - (g0*I2 + g1*C1)) and Z(Go - (g2*C2 + g3*C3)) and Z(C1*Ge*C1 - Ge) and Z(C1*Go*C1 + Go)
    U = exp_series(G, t, 3)
    loop = coeffs(C1*U*C1*U, t, 2)
    out['G1_cut_loop_series'] = Z(loop[0] - I2) and Z(loop[1] - 2*Ge) and Z(loop[2] - (2*Ge*Ge + Ge*Go - Go*Ge))
    X = loop[1]*t + loop[2]*t**2
    logl = coeffs(X - X*X/2, t, 2)
    out['G1_log_of_the_cut_loop'] = Z(logl[1] - 2*Ge) and Z(logl[2] - (Ge*Go - Go*Ge))
    # and on a larger carrier with exact rationals
    random.seed(41)
    J4 = sp.diag(1, 1, -1, -1)
    G4 = sp.Matrix(4, 4, lambda i, j: sp.Rational(random.randint(-4, 4), random.randint(1, 3)))
    Ge4, Go4 = (G4 + J4*G4*J4)/2, (G4 - J4*G4*J4)/2
    U4 = exp_series(G4, t, 3)
    l4 = coeffs(J4*U4*J4*U4, t, 2)
    out['G1_larger_carrier'] = Z(l4[1] - 2*Ge4) and Z(l4[2] - (2*Ge4*Ge4 + Ge4*Go4 - Go4*Ge4)) and Ge4[:2, 2:] == sp.zeros(2) and Go4[:2, :2] == sp.zeros(2)

    # G2: instances already on the line
    eta, k = sp.symbols('eta k', positive=True)
    Ub = (sp.cosh(eta)*I2 + sp.sinh(eta)*C2)                      # flow of the cut-odd generator C2 (a boost)
    E, O = (Ub + C1*Ub*C1)/2, (Ub - C1*Ub*C1)/2
    out['G2_reversal_and_split'] = Z(C1*Ub*C1 - Ub.subs(eta, -eta)) and Z((E*E - O*O).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp))) - I2)
    out['G2_seen_and_lost_of_DG1'] = zz((1/sp.cosh(eta)**2 + sp.tanh(eta)**2 - 1).rewrite(sp.exp))      # E^2 - O^2 = 1 divided by E^2
    # the two-valued chain of DU1: its block is the flow at time k* (the dual coupling)
    a, b = sp.exp(k), sp.exp(-k)
    ks = -sp.log(sp.tanh(k))/2
    B = a*I2 + b*C2
    flow = sp.sqrt(a**2 - b**2)*(sp.cosh(ks)*I2 + sp.sinh(ks)*C2)
    out['G2_chain_block_is_the_flow_at_the_dual_coupling'] = all(abs(complex((B - flow)[i, j].subs(k, v))) < 1e-12 for i in range(2) for j in range(2) for v in (0.3, 1.1, 2.7)) and zz((sp.tanh(ks) - sp.exp(-2*k)).rewrite(sp.exp))
    # (3.13): join x cut = cut at double time: the doubling of the cell
    Jt, Ct = 2*sp.cosh(eta), 2*sp.sinh(eta)
    out['G2_join_times_cut_is_doubling'] = zz((Jt*Ct - 2*sp.sinh(2*eta)).rewrite(sp.exp)) and zz((Jt**2 - Ct**2 - 4).rewrite(sp.exp))
    # the centre of the turn block (CT1): J: c -> -c ; weight exp(kappa c) ; even part cosh, odd part sinh
    c, kap = sp.symbols('c kappa', real=True)
    w = sp.exp(kap*c)
    out['G2_centre_split_of_CT1'] = zz((w + w.subs(c, -c))/2 - sp.cosh(kap*c)) and zz((w - w.subs(c, -c))/2 - sp.sinh(kap*c))

    # G3: one generator, three sectors.  G = (turn about the cut) + (boost across it) = i a C1 + b C2
    aa, bb, om = sp.symbols('a b omega', positive=True)
    Gs = sp.I*aa*C1 + bb*C2
    out['G3_square_is_a_number'] = Z(Gs*Gs - (bb**2 - aa**2)*I2)
    s = sp.symbols('s')                                           # s = sinh(omega)/omega, ch = cosh(omega), omega^2 = b^2 - a^2
    ch = sp.symbols('ch')
    Us = ch*I2 + s*Gs                                             # exp(G) for G^2 = omega^2
    rel = {ch**2: 1 + (bb**2 - aa**2)*s**2}                        # cosh^2 - omega^2 s^2 = 1
    det = sp.expand(Us.det()).subs(rel)
    out['G3_flow_has_unit_determinant'] = zz(sp.expand(Us.det() - (ch**2 - (bb**2 - aa**2)*s**2)))
    seen_inv = sp.expand(Us[0, 0]*Us[1, 1])                       # the two sheets' own entries
    lost_over_seen = sp.expand(Us[0, 1]*Us[1, 0])                 # what crosses, there and back
    out['G3_lost_over_seen'] = zz(lost_over_seen - bb**2*s**2) and zz(seen_inv - ch**2 - aa**2*s**2) and zz((seen_inv - 1 - lost_over_seen).subs(ch**2, 1 + (bb**2 - aa**2)*s**2))
    # the three sectors: s = sinh(w)/w (b > a), 1 (b = a), sin(W)/W (a > b)
    x = sp.symbols('x', positive=True)
    out['G3_boost_sector'] = zz((bb**2*(sp.sinh(om)/om)**2).subs(om, bb) - sp.sinh(bb)**2)
    out['G3_shear_sector'] = sp.limit(sp.sinh(om)/om, om, 0) == 1
    W = sp.symbols('Omega', positive=True)
    out['G3_turn_sector_closes_at_whole_half_turns'] = zz((sp.sinh(sp.I*W)/(sp.I*W)) - sp.sin(W)/W) and all((sp.sin(n*sp.pi)/(n*sp.pi) == 0) and (sp.cos(n*sp.pi) == (-1)**n) for n in (1, 2, 3, 4))
    # exact numbers through mpmath-free substitution: compare with the matrix exponential
    okn = True
    for (av, bv) in ((0, sp.Rational(3, 4)), (sp.Rational(1, 2), sp.Rational(1, 2)), (sp.Rational(5, 4), sp.Rational(1, 3)), (sp.Rational(1, 3), sp.Rational(7, 5))):
        Ue = (sp.I*av*C1 + bv*C2).exp()
        w2 = bv**2 - av**2
        sv = 1 if w2 == 0 else (sp.sinh(sp.sqrt(w2))/sp.sqrt(w2))
        okn &= abs(complex(sp.N(Ue[0, 1]*Ue[1, 0] - bv**2*sv**2, 30))) < 1e-20
    out['G3_checked_against_the_exponential'] = okn

    # G4: the cut loop, exactly:  J U J U = alpha + beta G_e + gamma [G_e, G_o],  gamma = s^2 ;  lost/seen = b^2 gamma
    Ge_, Go_ = sp.I*aa*C1, bb*C2
    Um = ch*I2 + s*(Ge_ - Go_)                                    # J U J
    H = (Um*Us).applyfunc(sp.expand)
    comm = Ge_*Go_ - Go_*Ge_
    alpha = ch**2 - (aa**2 + bb**2)*s**2
    out['G4_cut_loop_closed_form'] = Z(H - (alpha*I2 + 2*ch*s*Ge_ + s**2*comm))
    out['G4_curvature_lies_along_the_third_cut'] = Z(comm + 2*aa*bb*C3)
    out['G4_curvature_coefficient_is_lost_over_seen'] = zz(lost_over_seen - bb**2*s**2)
    out['G4_loop_is_identity_iff_odd_or_closed'] = zz(alpha.subs(ch**2, 1 + (bb**2 - aa**2)*s**2) - (1 - 2*aa**2*s**2))

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'CG1_RESULT.json'), 'w'), indent=1, default=str)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
