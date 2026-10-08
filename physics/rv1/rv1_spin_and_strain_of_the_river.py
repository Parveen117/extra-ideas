"""RV1: the frame's river has a strain and a spin; the law reads only the strain, a loop and a gyroscope read the spin.

The owner's statement: Coriolis is a coupling of the radial and the tangential; a new dimension coupled with
the frame of measurement produces curvature - native or local.

Sources read (unchanged):
  this line: TP1 (law Q = I1/4 + I2/2 - I3 on the order defect), MO1 (flat slices, one time, count form),
             CO1 (lapse equation), CV1 (connection of the frame), RMG9 (R-sector: quadratic form blind to w, cycle 2w*Area),
             DM1 (the third cut), GR2 (gyroscope).
Frame class: e_0 = d_t + v.d_X, e_i = d_i  (count dtau^2 = dt^2 - |dX - v dt|^2).  Symbolic algebra with sympy.  Python 3.12.
"""
import json
import pathlib
import sys

import sympy as sp

here = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(here.parent / 'cv1'))
sys.path.insert(0, str(here.parent / 'tp1'))
import cv1_frame_curvature as cv  # noqa: E402
import tp1_frame_defect_law as tp  # noqa: E402

t, x, y, z = sp.symbols('t x y z', real=True)
XS = (x, y, z)
COORDS = (t, x, y, z)


def river_frame(v, lapse=1):
    return sp.Matrix([[sp.Integer(1)/lapse, v[0]/lapse, v[1]/lapse, v[2]/lapse],
                      [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def grad(v):
    """A[i][j] = d_i v_j"""
    return sp.Matrix(3, 3, lambda i, j: sp.diff(v[j], XS[i]))


def strain_spin(v):
    A = grad(v)
    return (A + A.T)/2, (A - A.T)/2


def curl(v):
    return [sp.diff(v[2], y) - sp.diff(v[1], z), sp.diff(v[0], z) - sp.diff(v[2], x), sp.diff(v[1], x) - sp.diff(v[0], y)]


def div(v):
    return sum(sp.diff(v[i], XS[i]) for i in range(3))


def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]


def law_from_frame(v, lapse=1):
    I1, I2, I3 = tp.invariants(tp.structure_in(COORDS, river_frame(v, lapse)))
    return sp.simplify(sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3)


def law_from_strain(v):
    S, _ = strain_spin(v)
    return sum(S[i, j]**2 for i in range(3) for j in range(3)) - S.trace()**2


def zero(e):
    return sp.simplify(e) == 0


def run():
    res = {}
    V = [sp.Function(f'v{i}')(t, x, y, z) for i in (1, 2, 3)]
    rr = sp.sqrt(x*x + y*y + z*z)

    # T1 order defect = minus the gradient of the river; the law reads the strain only
    c = tp.structure_in(COORDS, river_frame(V))
    A = grad(V)
    for i in range(3):
        for j in range(3):
            if not zero(c[0][i + 1][j + 1] + A[i, j]):
                raise ValueError('order defect is not the gradient of the river')
    if not zero(law_from_frame(V) - law_from_strain(V)):
        raise ValueError('law is not strain^2 - expansion^2')
    res['law'] = 'Q = S:S - (tr S)^2, no spin term'

    # T2 a rigid turn of the river: order defect not zero, each invariant not zero, the law zero
    Om = sp.Symbol('Omega', positive=True)
    rigid = [-Om*y, Om*x, 0]
    I1, I2, I3 = tp.invariants(tp.structure_in(COORDS, river_frame(rigid)))
    if not (zero(I1 - 4*Om**2) and zero(I2 + 2*Om**2) and zero(I3) and zero(law_from_frame(rigid))):
        raise ValueError('rigid turn')
    if not zero(law_from_strain([V[i] + rigid[i] for i in range(3)]) - law_from_strain(V)):
        raise ValueError('law changes under a rigid turn')
    res['rigid_turn'] = {'I1': '4 Omega^2', 'I2': '-2 Omega^2', 'I3': '0', 'law': '0'}

    # T3 path of a slow free reading: pull = grad(v^2/2), sideways push = curl v x velocity
    U = sp.symbols('u1 u2 u3', real=True)
    B = curl(V)
    BxU = cross(B, U)
    v2 = sum(vi*vi for vi in V)/2
    for i in range(3):
        from_count = sp.diff(V[i], t) + sum(U[j]*sp.diff(V[i], XS[j]) for j in range(3)) \
            - sum((U[j] - V[j])*sp.diff(V[j], XS[i]) for j in range(3))
        if not zero(from_count - (sp.diff(V[i], t) + sp.diff(v2, XS[i]) + BxU[i])):
            raise ValueError('path law')
    rs = sp.Symbol('r_s', positive=True)
    fall = [-sp.sqrt(rs/rr)*q/rr for q in XS]
    if not all(zero(sp.diff(sum(f*f for f in fall)/2, XS[i]) + rs*XS[i]/(2*rr**3)) for i in range(3)) or not all(zero(b) for b in curl(fall)):
        raise ValueError('fall: pull is not r_s/2r^2 inward, or it has spin')
    if not (zero(sp.diff(sum(f*f for f in rigid)/2, x) - Om**2*x) and curl(rigid) == [0, 0, 2*Om]):
        raise ValueError('rigid: pull is not Omega^2 rho outward, or spin is not 2 Omega')
    res['path'] = 'acceleration = d_t v + grad(v^2/2) + (curl v) x velocity'

    # T4 gyroscope: the frame's connection along e_0 turns the triad at the spin of the river
    w = cv.connection(c)
    _, W = strain_spin(V)
    sign = None
    for i in range(3):
        for j in range(3):
            if i != j:
                if zero(w[i + 1][j + 1][0] - W[i, j]):
                    s = 1
                elif zero(w[i + 1][j + 1][0] + W[i, j]):
                    s = -1
                else:
                    raise ValueError('triad does not turn at the spin')
                if sign not in (None, s):
                    raise ValueError('inconsistent sign')
                sign = s
    res['gyroscope'] = f'turn of the triad along e_0 = {sign:+d} x spin W'

    # T5 loop: exact two-way time difference around a circle in a rigidly turning river
    rho = sp.Symbol('rho', positive=True)
    lag = sp.simplify(2*sp.pi/(1/rho - Om) - 2*sp.pi/(1/rho + Om))
    flux = 2*Om*sp.pi*rho**2
    if not zero(lag - 2*flux/(1 - (Om*rho)**2)):
        raise ValueError('loop lag')
    res['loop'] = 'two-way lag = 2 x flux(curl v) / (1 - beta^2)'

    # T6 the law's equation for the river itself (variation of v): curl curl v = 0
    Vs = [sp.Function(f'v{i}')(x, y, z) for i in (1, 2, 3)]
    Qs = law_from_strain(Vs)
    cc = curl(curl(Vs))
    for i in range(3):
        el = -sum(sp.diff(sp.diff(Qs, sp.Derivative(Vs[i], XS[j])), XS[j]) for j in range(3))
        if not zero(el - cc[i]):
            raise ValueError('river equation is not curl curl v')
    psi = sp.Function('psi')(x, y, z)
    if not all(zero(q) for q in curl(curl([sp.diff(psi, q) for q in XS]))):
        raise ValueError('irrotational rivers should be free')
    n = sp.Symbol('n')
    axial = lambda om: [-om*y, om*x, 0]
    cc_ax = curl(curl(axial(rr**n)))
    cond = sp.simplify(cc_ax[1]/(x*rr**(n - 2)))
    if not zero(sp.expand(cond) + n*(n + 3)) or not zero(cc_ax[2]):
        raise ValueError('axial exponents')
    res['axial_river_exponents'] = 'omega ~ r^n with n(n+3) = 0: n = 0 (rigid) or n = -3'

    # T7 the n = -3 river: its spin field is a dipole, with no source and no curl
    Jm = sp.Symbol('J', positive=True)
    drag = axial(2*Jm/rr**3)
    Bd = curl(drag)
    dip = [2*Jm*(3*z*q - (rr**2 if k == 2 else 0))/rr**5 for k, q in enumerate(XS)]
    if not all(zero(Bd[i] - dip[i]) for i in range(3)) or not zero(div(Bd)) or not all(zero(q) for q in curl(Bd)):
        raise ValueError('dipole')
    Sd, _ = strain_spin(drag)
    Qd = sp.simplify(sum(Sd[i, j]**2 for i in range(3) for j in range(3)) - Sd.trace()**2)
    if not zero(Qd - 18*Jm**2*(x*x + y*y)/rr**8):
        raise ValueError('drag strain')
    res['dragging_river'] = {'spin_field': 'dipole 2J(3 z X - r^2 z^)/r^5', 'strain_squared': '18 J^2 sin^2(theta) / r^6'}

    # T8 second order: with flat slices and one time, radial fall + this river cannot keep Q = 0
    bfun = sp.Function('b')(rr)
    both = [-bfun*q/rr + d for q, d in zip(XS, drag)]
    Qb = sp.simplify(law_from_strain(both))
    radial_part = sp.simplify(law_from_strain([-bfun*q/rr for q in XS]))
    if not zero(Qb - radial_part - Qd):
        raise ValueError('cross term')
    r_ = sp.Symbol('r', positive=True)
    bf = sp.Function('b')(r_)
    target = -(2/r_**2)*sp.diff(r_*bf**2, r_)
    if not zero(radial_part.subs({y: 0, z: 0, x: r_}).doit() - target.doit()):
        raise ValueError('radial part')
    res['second_order'] = 'Q = -(2/r^2)(r b^2)\' + 18 J^2 sin^2(theta)/r^6 : no b(r) makes it vanish at every angle'

    # T9 in a plane the river's gradient is RMG9's response element L = m(1 + uK + vS) + wR
    mm, uu, vv, ww = sp.symbols('m u v w', real=True)
    Rm, Km = sp.Matrix([[0, -1], [1, 0]]), sp.diag(1, -1)
    Sm = Rm*Km
    Lm = mm*(sp.eye(2) + uu*Km + vv*Sm) + ww*Rm
    lin = list(Lm*sp.Matrix([x, y])) + [0]                 # the linear river v = L X
    if not zero(curl(lin)[2] - 2*ww) or not zero(div(lin) - 2*mm):
        raise ValueError('plane element: spin or expansion')
    Lsym = (Lm + Lm.T)/2
    if not zero(law_from_strain(lin) + 2*Lsym.det()) or not zero(law_from_frame(lin) + 2*mm**2*(1 - uu**2 - vv**2)):
        raise ValueError('plane law')
    if not zero(Lm.det() - Lsym.det() - ww**2):
        raise ValueError('determinant split')
    res['plane'] = 'Q = -2 det(symmetric part) = -2 m^2 (1 - u^2 - v^2); w absent; curl = 2w; expansion = 2m'

    # lapse scaling (as CO1): Q ~ 1/N^2
    Nn = sp.Symbol('N', positive=True)
    if not zero(law_from_frame(rigid, Nn)) or not zero(law_from_frame([Om*x, Om*y, Om*z], Nn) + 6*Om**2/Nn**2):
        raise ValueError('lapse scaling')

    # numbers (illustration, floats)
    G, cl, JE, r0, yr = 6.674e-11, 2.998e8, 5.86e33, 7.02e6, 3.156e7
    rate = G*JE/(2*cl**2*r0**3)
    res['numbers'] = {'gyroscope_polar_orbit_7020km_mas_per_year': round(rate*yr*206264806.2, 1),
                      'loop_1m2_at_pole_two_way_lag_s': 4*7.292e-5*1.0/cl**2}
    return res


if __name__ == '__main__':
    out = run()
    with open('RV1_RESULT.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
