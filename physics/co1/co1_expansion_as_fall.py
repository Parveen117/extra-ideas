"""CO1: an expanding space read as a fall frame; the line's law reproduces the coupled solution of the research branch.

Sources read (unchanged):
  Publications research branch, papers/ugd-kahler-propagation/COUPLED_COSMOLOGICAL_SOLUTION.md
      CS.7 (n = n0/v), CS.9 (rho, p), CS.10 (lapse form), CS.11, CS.13, CS.15
  this line: TP1 (law Q = I1/4 + I2/2 - I3 on the order defect), MO1 (fall speed, beta^2 = r_s/r), GB1, EN1.
Symbolic algebra with sympy.  Python 3.12.
"""
import json
import pathlib
import sys

import sympy as sp

here = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(here.parent / 'cv1'))
sys.path.insert(0, str(here.parent / 'tp1'))
import tp1_frame_defect_law as tp  # noqa: E402

t, x, y, z = sp.symbols('t x y z', real=True)
X = (t, x, y, z)
a = sp.Function('a', positive=True)(t)
N = sp.Function('N', positive=True)(t)
h = sp.diff(a, t)/a


def law(E):
    I1, I2, I3 = tp.invariants(tp.structure_in(X, E))
    return sp.simplify(sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3)


def comoving_frame(lapse=1):
    return sp.Matrix([[1/lapse, 0, 0, 0], [0, 1/a, 0, 0], [0, 0, 1/a, 0], [0, 0, 0, 1/a]])


def fall_frame():
    """proper distances X = a x; flat slices, one time, e0 = d_t + h X.d_X"""
    return sp.Matrix([[1, h*x, h*y, h*z], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def run():
    res = {}
    # T1 the comoving frame, written in proper distances, is a fall frame (one frame, two coordinate systems)
    Qc, Qf = law(comoving_frame()), law(fall_frame())
    if sp.simplify(Qc - Qf) != 0:
        raise ValueError('the two forms are not the same frame')
    if sp.simplify(Qc + 6*h**2) != 0:
        raise ValueError('law value is not -6 h^2')
    res['law_value'] = '-6 h^2 in both frames'      # same frame field, so this agreement is coordinate invariance

    # T2 lapse and scale variations of  N a^3 [ Q/(2 kappa) - Lambda/kappa - rho(a) ]  give CS.11 with CS.9
    kappa, Lam, m, n0 = sp.symbols('kappa Lambda m n_0', positive=True)
    QN = law(comoving_frame(N))
    if sp.simplify(QN - Qc.subs(sp.diff(a, t), sp.diff(a, t)/N)) != 0:
        raise ValueError('lapse dependence')
    A, Ad, Nn = sp.symbols('A Adot Nn', positive=True)
    rho_of = lambda s: m*n0/s**3 + sp.Rational(3, 16)*kappa*n0**2/s**6
    Lag = Nn*A**3*(-6*(Ad/(A*Nn))**2/(2*kappa) - Lam/kappa - rho_of(A))
    C = sp.simplify(-sp.diff(Lag, Nn).subs(Nn, 1)*kappa/A**3)
    hh = Ad/A
    if sp.simplify(C - (Lam + kappa*rho_of(A) - 3*hh**2)) != 0:
        raise ValueError('lapse equation is not CS.11')
    # scale equation: d/dt dL/dAdot - dL/dA = 0 at N = 1, with Addot free
    Add = sp.Symbol('Addot')
    L1 = Lag.subs(Nn, 1)
    EL = sp.diff(sp.diff(L1, Ad), A)*Ad + sp.diff(sp.diff(L1, Ad), Ad)*Add - sp.diff(L1, A)
    hdot = Add/A - hh**2
    target = 2*hdot + 3*hh**2 - Lam + kappa*sp.Rational(3, 16)*kappa*n0**2/A**6
    if sp.simplify(EL*kappa/(-3*A**2) - target) != 0:
        raise ValueError('scale equation is not CS.11')
    p_from_action = sp.simplify(-rho_of(A) - A*sp.diff(rho_of(A), A)/3)
    if sp.simplify(p_from_action - sp.Rational(3, 16)*kappa*n0**2/A**6) != 0:
        raise ValueError('pressure is not CS.9')
    res['pressure'] = str(p_from_action)

    # T3 the fall law shell by shell: beta = h X, beta^2 = r_s(X)/X with r_s = (Lambda + kappa rho) X^3 / 3
    Xr, rho = sp.symbols('X rho', positive=True)
    h2 = (Lam + kappa*rho)/3
    rs = (Lam + kappa*rho)*Xr**3/3
    if sp.simplify(h2*Xr**2 - rs/Xr) != 0:
        raise ValueError('fall law')
    G, M = sp.symbols('G M', positive=True)
    if sp.simplify((kappa*rho*Xr**3/3).subs({kappa: 8*sp.pi*G, rho: M/(sp.Rational(4, 3)*sp.pi*Xr**3)}) - 2*G*M) != 0:
        raise ValueError('enclosed mass')

    # T4 CS.13 and CS.15 from the lapse equation; the perfect-square case
    om = sp.sqrt(3*Lam)
    d = kappa*m*n0/(2*Lam)
    b = 3*kappa*n0/(4*om)
    v = d*(sp.cosh(om*t) - 1) + b*sp.sinh(om*t)
    lhs = sp.diff(v, t)**2
    rhs = 3*Lam*v**2 + 3*kappa*m*n0*v + sp.Rational(9, 16)*kappa**2*n0**2
    if sp.simplify((lhs - rhs).rewrite(sp.exp)) != 0:
        raise ValueError('CS.15 does not solve CS.13')
    hv = sp.diff(v, t)/(3*v)
    rho_v = m*n0/v + sp.Rational(3, 16)*kappa*n0**2/v**2
    if sp.simplify((3*hv**2 - Lam - kappa*rho_v).rewrite(sp.exp)) != 0:
        raise ValueError('CS.15 does not satisfy the lapse equation')
    disc = sp.simplify((3*kappa*m*n0)**2 - 4*3*Lam*sp.Rational(9, 16)*kappa**2*n0**2)
    if sp.simplify(disc.subs(m, om/2)) != 0 or sp.simplify((d - b).subs(m, om/2)) != 0:
        raise ValueError('perfect-square case')
    res['perfect_square_case'] = 'm = sqrt(3 Lambda)/2  ->  v = b (exp(omega t) - 1)'

    # T5 constant rate: redshift = (Doppler of the fall) x (1 / clock factor)
    be = sp.symbols('beta', positive=True)
    H0, Xe = sp.symbols('H X_e', positive=True)
    T = sp.integrate(1/(1 - H0*Xr), (Xr, 0, Xe), conds='none')
    one_plus_z = sp.simplify(sp.exp(H0*T).subs(Xe, be/H0))
    if sp.simplify(one_plus_z - 1/(1 - be)) != 0:
        raise ValueError('redshift')
    if sp.simplify((1/(1 - be))**2 - ((1 + be)/(1 - be))*(1/(1 - be**2))) != 0:
        raise ValueError('factorisation')
    res['redshift_constant_rate'] = '1 + z = 1/(1 - beta) = sqrt((1+beta)/(1-beta)) / sqrt(1 - beta^2)'
    return res


if __name__ == '__main__':
    out = run()
    with open('CO1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
