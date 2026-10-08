"""CI1: one measure for both sides - the information shared between the two readings of an element.

The owner's statement: measure space-time and thermo against something else that is present in both - not
against time; an invariant, something like entropy and the arrow of time.

For an element M = a + bK + cS + dR (EMK-1: det = (a^2 - b^2) + (d^2 - c^2) = seen channel + lost channel):
    I_cut  = (1/2) Log( seen / det )            - what the cut along K does not hold by itself
    I_0    = (1/2) Log( a^2 / det )             - the same, without choosing a cut
Sources read (unchanged):
  Publications papers/emk-ugd-algebra EMK-1 (determinant channels), coherence-first-thermodynamics CF-3 (7)
  response-geometry RMG1 (rapidity l, tanh(l/2) = rho), RMG2 T2-T3 (tower), RMG6 (memory weight r^2), RMG10 T1, T4 (det = m^2 (1 - delta); light cone)
  extra-ideas physics gb1 (unit block Exp(psi n), clock factor), ln1 (seen and lost), in1, gr1, ms1 (speed^2 + memory = 1);
  extra-ideas uncut up3 (chi), up4 (w), tt1 (cycle = boost)
No measurement.  sympy, exact; the logarithm is checked through its argument.  Python 3.12.
"""
import json
import math

import sympy as sp

I2 = sp.eye(2)
R = sp.Matrix([[0, -1], [1, 0]])
K = sp.Matrix([[1, 0], [0, -1]])
S = R*K


def parts(M):
    return ((M[0, 0] + M[1, 1])/2, (M[0, 0] - M[1, 1])/2, (M[0, 1] + M[1, 0])/2, (M[1, 0] - M[0, 1])/2)


def e2I_cut(M):
    """exp(2 I_cut) = seen / det, for the cut along K"""
    a, b, c, d = parts(M)
    return (a*a - b*b)/M.det()


def e2I_0(M):
    a = parts(M)[0]
    return a*a/M.det()


def zero(e):
    return sp.simplify(sp.expand(sp.sympify(e).rewrite(sp.exp))) == 0


def run():
    res = {}
    A, B, C = sp.symbols('A B C', positive=True)
    H = sp.Matrix([[A, B], [B, C]])

    # T1 thermo side: exp(2 I_cut) = AC/det = C_P/C_V = 1/(1 - r^2) ; r^2 is also the squared correlation of the two fluctuations
    if not zero(e2I_cut(H) - A*C/(A*C - B*B)):
        raise ValueError('thermo cut value')
    T_, V_ = sp.symbols('T V', positive=True)
    det = A*C - B*B
    CV, CP = T_/A, T_*C/det                                   # CF-3 (7)
    if not zero(e2I_cut(H) - CP/CV):
        raise ValueError('capacity ratio')
    Sig = H.inv()
    if not zero(Sig[0, 1]**2/(Sig[0, 0]*Sig[1, 1]) - B*B/(A*C)) or not zero(1/e2I_cut(H) - (1 - B*B/(A*C))):
        raise ValueError('correlation')
    res['thermo'] = 'exp(2 I_cut) = C_P/C_V = 1/(1 - r^2), r^2 = B^2/(AC) = squared correlation of the two fluctuations'

    # T2 the cut can be turned: I_cut runs from 0 (cut on the principal axes) to I_0 (cut on the diagonal)
    m, rho, th = sp.symbols('m rho theta', positive=True)
    Hp = m*(I2 + rho*K)
    W = sp.cos(th)*I2 + sp.sin(th)*R
    Ht = sp.simplify(W*Hp*W.inv())
    val = sp.simplify(e2I_cut(Ht))
    if not zero(val - (1 - rho**2*sp.cos(2*th)**2)/(1 - rho**2)):
        raise ValueError('turned cut')
    if not (zero(val.subs(th, 0) - 1) and zero(val.subs(th, sp.pi/4) - 1/(1 - rho**2)) and zero(e2I_0(Ht) - 1/(1 - rho**2))):
        raise ValueError('ends of the range')
    lp, lm = m*(1 + rho), m*(1 - rho)
    if not zero(e2I_0(Hp) - ((lp + lm)/2)**2/(lp*lm)):
        raise ValueError('mean form')
    res['range'] = 'I_cut from 0 (principal cut) to I_0 (diagonal cut); exp(I_0) = (mean of the two principal responses)/(their geometric mean)'

    # T3 space-time side: unit block Exp(psi n): exp(I_0) = cosh(psi) = 1/(clock factor)
    psi = sp.Symbol('psi', real=True)
    for n in (K, S, (3*K + 4*S)/5):
        G = sp.cosh(psi)*I2 + sp.sinh(psi)*n
        if not zero(G.det() - 1) or not zero(e2I_0(G) - sp.cosh(psi)**2):
            raise ValueError('unit block')
    Gs = sp.cosh(psi)*I2 + sp.sinh(psi)*S
    Gk = sp.cosh(psi)*I2 + sp.sinh(psi)*K
    if not zero(e2I_cut(Gs) - sp.cosh(psi)**2) or not zero(e2I_cut(Gk) - 1):
        raise ValueError('cut across / along the block')
    x = sp.Symbol('x', positive=True)                      # x = tanh(psi) = speed; x^2 = r_s/r for the fall
    if not zero((sp.cosh(psi)**2).subs(psi, sp.atanh(x)) - 1/(1 - x*x)):
        raise ValueError('speed form')
    ser = sp.series(-sp.log(1 - x)/2, x, 0, 3).removeO()   # with x -> r_s/r
    if sp.expand(ser - (x/2 + x*x/4)) != 0:
        raise ValueError('weak field')
    res['space_time'] = 'exp(I_0) = cosh(psi) = 1/N ; exp(-2 I_0) = 1 - speed^2 = 1 - r_s/r for the fall ; I_0 = r_s/2r + ... far away'

    # T4 the arrow: same-sense composition never loses it; a tower generation and a diagram cycle at least double it
    p1, p2 = sp.symbols('psi_1 psi_2', positive=True)
    gain = sp.cosh(p1 + p2) - sp.cosh(p1)*sp.cosh(p2)
    if not zero(gain - sp.sinh(p1)*sp.sinh(p2)):
        raise ValueError('composition')
    if not zero(sp.cosh(2*p1) - sp.cosh(p1)**2 - sp.sinh(p1)**2):
        raise ValueError('doubling')
    lam = sp.Symbol('lambda', real=True)
    kappa = (1 + 2*lam*m)/(1 + lam*m*(1 + rho**2))          # RMG2-T2
    if not zero(kappa - 1 - lam*m*(1 - rho**2)/(1 + lam*m*(1 + rho**2))):
        raise ValueError('tower factor')
    Hn = sp.Matrix([[3, 1], [1, 2]])
    up, down = Hn + sp.Rational(1, 5)*Hn*Hn, Hn - sp.Rational(1, 20)*Hn*Hn
    if not (e2I_0(up) > e2I_0(Hn) > e2I_0(down) > 1):
        raise ValueError('tower direction')
    res['arrow'] = {'composition': 'cosh(a + b) - cosh(a) cosh(b) = sinh(a) sinh(b) >= 0 for the same sense',
                    'generation_or_cycle': 'I_0(2 psi) >= 2 I_0(psi)',
                    'tower': 'lambda > 0 raises I_0, lambda < 0 lowers it, lambda = 0 keeps it',
                    'witness_exp_2I0': [str(e2I_0(down)), str(e2I_0(Hn)), str(e2I_0(up))]}

    # T5 with a turn-part: exp(-2 I_0) = 1 - rho^2 + t^2 ; the sign of I_0 is the causal sector
    u, v, t = sp.symbols('u v t', real=True)
    L = m*(I2 + u*K + v*S) + m*t*R
    if not zero(1/e2I_0(L) - (1 - u*u - v*v + t*t)):
        raise ValueError('with a turn-part')
    signs = {}
    for name, (uu, tt) in {'cut': (sp.Rational(3, 5), sp.Rational(1, 5)), 'cone': (sp.Rational(3, 5), sp.Rational(3, 5)),
                           'turn': (sp.Rational(1, 5), sp.Rational(3, 5))}.items():
        signs[name] = str(e2I_0(L.subs({u: uu, v: 0, t: tt, m: 1})))
    if not (sp.Rational(signs['cut']) > 1 and sp.Rational(signs['cone']) == 1 and sp.Rational(signs['turn']) < 1):
        raise ValueError('sectors')
    res['sectors'] = {'exp_2I0': signs, 'reading': 'I_0 > 0 where the response has a cut, 0 on the light cone, < 0 where it is a turn'}

    # numbers (illustration, floats)
    res['numbers'] = {'gas_5_3_I_cut': round(0.5*math.log(5/3), 4),
                      'same_I_as_speed': round(math.sqrt(1 - 3/5), 4),
                      'same_I_as_radius_in_r_s': 2.5,
                      'earth_surface_I': 6.96e-10, 'sun_surface_I': 2.12e-6}
    return res


if __name__ == '__main__':
    out = run()
    with open('CI1_RESULT.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
