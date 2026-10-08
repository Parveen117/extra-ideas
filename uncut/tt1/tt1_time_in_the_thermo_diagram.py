"""TT1: where time sits in the thermo diagram - the turn-part of the cut is the timelike direction, and one
cycle of the diagram with turn-part w is a boost.

The owner's statement: relativity has become a special case and was connected even before thermo; now use the
thermo space together with time in place of space, and see what relation between time and thermo comes out.

Sources read (unchanged):
  response-geometry RMG10 T2, T4, T6 (response space = SL(2,R), signature (2,1), t = w/m timelike; causal sectors;
                    the tower keeps t/rho), RMG11 T3 (the three native exponentials are the geodesics), RMG9 T1, T2
  extra-ideas uncut/up4 (w), up6 (cross-corner pair), up7 (cuts: u^2 + v^2 - w^2 = 1; four-step walk)
  extra-ideas physics/rb1 (return = boost), lt1 (tower generation x -> 2x/(1+x^2)), gb1 (unit block Exp(psi n)), tl1
No measurement, no physical constant.  sympy, exact.  Python 3.12.
"""
import json

import sympy as sp

I2 = sp.eye(2)
R = sp.Matrix([[0, -1], [1, 0]])
K = sp.Matrix([[1, 0], [0, -1]])
S = R*K


def parts(M):
    """coefficients (a, u, v, w) of M = a + uK + vS + wR"""
    return ((M[0, 0] + M[1, 1])/2, (M[0, 0] - M[1, 1])/2, (M[0, 1] + M[1, 0])/2, (M[1, 0] - M[0, 1])/2)


def Exp_turn(th):
    return sp.cos(th)*I2 + sp.sin(th)*R


def Exp_cut(eta, n):
    return sp.cosh(eta)*I2 + sp.sinh(eta)*n


def zero_m(M):
    return sp.simplify(M.applyfunc(lambda q: sp.expand(q.rewrite(sp.exp)))) == sp.zeros(*M.shape)


def tower(x):
    return 2*x/(1 + x*x)


def run():
    res = {}
    u, v, w, th, eta = sp.symbols('u v w theta eta', real=True)
    X = u*K + v*S + w*R

    # T1 the space of cuts has two space directions and one time direction; the three generators act as a turn and two boosts
    a, u1, v1, w1 = parts(sp.simplify(Exp_turn(th)*X*Exp_turn(-th)))
    if not (sp.simplify(w1 - w) == 0 and sp.simplify(u1 - (sp.cos(2*th)*u - sp.sin(2*th)*v)) == 0
            and sp.simplify(v1 - (sp.sin(2*th)*u + sp.cos(2*th)*v)) == 0):
        raise ValueError('turn')
    a, u2, v2, w2 = parts(sp.simplify(Exp_cut(eta, K)*X*Exp_cut(-eta, K)))
    if not (sp.simplify(u2 - u) == 0 and sp.simplify((v2 - (sp.cosh(2*eta)*v - sp.sinh(2*eta)*w)).rewrite(sp.exp)) == 0
            and sp.simplify((w2 - (-sp.sinh(2*eta)*v + sp.cosh(2*eta)*w)).rewrite(sp.exp)) == 0):
        raise ValueError('boost along K')
    a, u3, v3, w3 = parts(sp.simplify(Exp_cut(eta, S)*X*Exp_cut(-eta, S)))
    if not (sp.simplify(v3 - v) == 0 and sp.simplify((u3 - (sp.cosh(2*eta)*u + sp.sinh(2*eta)*w)).rewrite(sp.exp)) == 0
            and sp.simplify((w3 - (sp.sinh(2*eta)*u + sp.cosh(2*eta)*w)).rewrite(sp.exp)) == 0):
        raise ValueError('boost along S')
    for q in ((u1, v1, w1), (u2, v2, w2), (u3, v3, w3)):
        if sp.simplify((q[0]**2 + q[1]**2 - q[2]**2 - (u*u + v*v - w*w)).rewrite(sp.exp)) != 0:
            raise ValueError('invariant')
    res['cut_space'] = 'u^2 + v^2 - w^2 kept ; Exp(theta R) turns (u, v) by 2 theta ; Exp(eta K), Exp(eta S) boost (v, w), (u, w) by 2 eta'

    # T2 one cycle of the diagram with turn-part w = sinh(eta) is the boost Exp(-2 eta n)
    c, s = sp.symbols('c s', real=True)
    n0 = c*K + s*S
    kap = sp.cosh(eta)*n0 + sp.sinh(eta)*R
    chi = R*kap
    n = R*n0
    on = {s: sp.sqrt(1 - c*c)}
    if not zero_m((kap*kap - I2).subs(on)) or not zero_m((n*n - I2).subs(on)):
        raise ValueError('cut')
    if not zero_m((chi*chi - Exp_cut(-2*eta, n)).subs(on)):
        raise ValueError('cycle is not the boost')
    if not zero_m((chi**4 - Exp_cut(-4*eta, n)).subs(on)):
        raise ValueError('two cycles')
    res['cycle'] = '(iota kappa)^2 = Exp(-2 eta n) with sinh(eta) = w : cycles add rapidity 2 eta each'

    # exact numbers: w = 3/4
    wv = sp.Rational(3, 4)
    ch, sh = sp.sqrt(1 + wv*wv), wv
    x1 = sh/ch
    x2 = tower(x1)
    clock = 1/(1 + 2*wv*wv)
    if not (x1 == sp.Rational(3, 5) and x2 == sp.Rational(15, 17) and clock == sp.Rational(8, 17)
            and sp.simplify(1 - x2*x2 - clock*clock) == 0):
        raise ValueError('numbers')
    kw = sp.Rational(5, 4)*K + wv*R
    cyc = (R*kw)**2
    a4, u4, v4, w4 = parts(cyc)
    if not (a4 == sp.Rational(17, 8) and u4*u4 + v4*v4 == sp.Rational(225, 64) and w4 == 0):
        raise ValueError('witness cycle')
    res['numbers_w_3_4'] = {'tanh_eta': '3/5', 'speed_of_one_cycle': '15/17 = tower(3/5)', 'clock_factor': '8/17 = 1/(1 + 2 w^2)'}

    # T3 for a response element: tanh(eta) = (B2 - B1) / sqrt((A - C)^2 + (B1 + B2)^2) = t / rho ; three cases
    A, B1, B2, C = sp.symbols('A B_1 B_2 C', real=True)
    L = sp.Matrix([[A, B1], [B2, C]])
    m_, uL, vL, wL = parts(L)
    D = sp.simplify(uL**2 + vL**2 - wL**2)
    if sp.simplify(D - (((A - C)/2)**2 + B1*B2)) != 0:
        raise ValueError('D')
    if sp.simplify(wL**2/(uL**2 + vL**2) - (B2 - B1)**2/((A - C)**2 + (B1 + B2)**2)) != 0:
        raise ValueError('tanh eta')
    cases = {}
    for name, (Av, B1v, B2v, Cv) in {'cut': (3, 1, 2, 1), 'null': (2, 0, 1, 2), 'turn': (1, -2, 2, 1)}.items():
        Lm = sp.Matrix([[Av, B1v], [B2v, Cv]])
        Xm = Lm - Lm.trace()/2*I2
        Dv = (Xm*Xm)[0, 0]
        cases[name] = str(Dv)
    if not (sp.Rational(cases['cut']) > 0 and sp.Rational(cases['null']) == 0 and sp.Rational(cases['turn']) < 0):
        raise ValueError('cases')
    res['response'] = {'tanh_eta': '(B_2 - B_1)/sqrt((A - C)^2 + (B_1 + B_2)^2) = t/rho',
                       'D>0': 'the response has a cut; its cycle is a boost', 'D=0': 'the cut degenerates (light cone)',
                       'D<0': 'no cut: the traceless part is a turn; the cycle is a rotation'}

    # T4 the tower does not change it (RMG10-T6)
    lam = sp.Symbol('lambda', real=True)
    L2 = L + lam*L*L
    m2, u2L, v2L, w2L = parts(sp.expand(L2))
    if sp.simplify(w2L**2*(uL**2 + vL**2) - wL**2*(u2L**2 + v2L**2)) != 0:
        raise ValueError('tower')
    res['tower'] = 'tanh(eta) is the same at every level of the tower'

    # T5 the cross-corner pair of UP6 on a gas with constant capacities
    a_, b_ = sp.symbols('a b', positive=True)
    cc, Rg = sp.symbols('c R_g', positive=True)
    U = b_**(-Rg/cc)*sp.exp(a_/cc)
    br = lambda f, g: sp.diff(f, a_)*sp.diff(g, b_) - sp.diff(f, b_)*sp.diff(g, a_)
    Tg = sp.diff(U, a_)
    DS = lambda f: a_*sp.diff(f, a_)
    DV = lambda f: b_*br(f, Tg)/br(b_, Tg)
    Ag, B1g, B2g, Cg = [sp.simplify(q/U) for q in (DS(DS(U)), DV(DS(U)), DS(DV(U)), DV(DV(U)))]
    x_gas = sp.simplify((B2g - B1g)/sp.sqrt((Ag - Cg)**2 + (B1g + B2g)**2))
    if sp.simplify(B1g - Rg/cc) != 0 or sp.simplify(B2g) != 0 or sp.simplify(Cg) != 0 or sp.simplify(Ag - a_*(a_ + cc)/cc**2) != 0:
        raise ValueError('gas pair: ' + str((Ag, B1g, B2g, Cg)))
    if sp.simplify(x_gas**2 - (Rg/cc)**2/(Ag**2 + (Rg/cc)**2)) != 0:
        raise ValueError('gas rapidity')
    x2 = sp.simplify(x_gas**2)
    if sp.simplify(x2 - (Rg*cc)**2/(a_**2*(a_ + cc)**2 + (Rg*cc)**2)) != 0:
        raise ValueError('gas closed form')
    if sp.limit(x2, a_, sp.oo) != 0 or sp.limit(x2, a_, 0) != 1:
        raise ValueError('gas limits')
    res['gas'] = {'closed_form': 'tanh^2(eta) = (R c)^2 / (S^2 (S + c)^2 + (R c)^2)', 'S_large': 'eta -> 0', 'S_to_0': 'tanh(eta) -> 1',
                  'A/U': str(Ag), 'B1/U': str(B1g), 'B2/U': str(B2g), 'C/U': str(Cg),
                  'tanh_eta_squared': str(sp.simplify(x_gas**2)), 'D': str(sp.simplify((Ag - Cg)**2/4 + B1g*B2g))}
    return res


if __name__ == '__main__':
    out = run()
    with open('TT1_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
