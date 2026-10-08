"""SW1: a centre that turns, to first order: the swirl of the fall frame from the law of CV1/TP1, and what it does to a carried direction.

Frame of the thesis with one more entry (the fall history also goes round the axis at the rate W(r)):
    e0 = d_t - beta d_r + W(r) d_phi ,  e1 = d_r ,  e2 = (1/r) d_theta ,  e3 = (1/(r sin theta)) d_phi ,  beta^2 = r_s/r .
Slices flat, one time (MO1's assumption, derived for the centre at rest in MA1; here it is put in and then tested).
Symbolic (sympy), with CV1's own functions.  Python 3.12."""
import json, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1'))
_cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'tp1'))
import cv1_frame_curvature as cv
import tp1_frame_defect_law as tp
os.chdir(_cwd)

t, r, th, ph = cv.X
rs, eps, a = sp.symbols('r_s epsilon a', positive=True)
W = sp.Function('W')(r)
z = lambda e: sp.simplify(e) == 0


def ricci(bb, swirl):
    E = sp.Matrix([[1, -bb, 0, swirl], [0, 1, 0, 0], [0, 0, 1/r, 0], [0, 0, 0, 1/(r*sp.sin(th))]])
    c = cv.structure(E); w = cv.connection(c); R = cv.curvature(E, c, w)
    return R, cv.contracted(R)


def run():
    out = {}
    b = sp.sqrt(rs/r)
    R, Ric = ricci(b, eps*W)
    co = lambda i, j, k: sp.simplify(sp.expand(Ric[i, j]).coeff(eps, k))
    law = (sp.diff(r**4*sp.diff(W, r), r))/r**3                    # = r W'' + 4 W'
    # T1: first order in the swirl: the only components are (time, round) and (radial, round), both the same operator
    first = {(i, j): co(i, j, 1) for i in range(4) for j in range(4)}
    out['T1_first_order_law'] = (z(first[(0, 3)] - law*sp.sin(th)/2) and z(first[(1, 3)] + b*law*sp.sin(th)/2)
                                 and all(z(v) for k, v in first.items() if set(k) not in ({0, 3}, {1, 3})))
    out['T1_order_zero_empty'] = all(z(co(i, j, 0)) for i in range(4) for j in range(4))
    sol = sp.dsolve(sp.Eq(law, 0), W).rhs
    C1, C2 = sp.symbols('C1 C2')
    out['T1_solution'] = z(sp.diff(sol*r**3, r, 2)*r - 2*sp.diff(sol*r**3, r)) or z(sp.diff(r**4*sp.diff(sol, r), r))
    out['T1_profile'] = z(law.subs(W, a/r**3).doit()) and z(law.subs(W, a).doit())
    # T2: a constant swirl is no field: with W constant every contracted component vanishes to all orders,
    #     and with r_s = 0 as well the whole curvature vanishes (flat space read in a turning frame)
    w0 = sp.Symbol('w0', positive=True)
    R0, Ric0 = ricci(b, w0)
    Rf, Ricf = ricci(sp.Integer(0), w0)
    out['T2_constant_swirl_no_field'] = all(z(v) for v in Ric0) and all(z(v) for v in Rf.values())
    # T3: second order: with flat slices the law is NOT met: the remainder is exact and not zero
    second = sp.Matrix(4, 4, lambda i, j: co(i, j, 2).subs(W, a/r**3).doit())
    q = 9*a**2*sp.sin(th)**2/(2*r**6)
    out['T3_second_order_remainder'] = z(second[0, 0] + q) and z(second[1, 1] - q) and z(second[3, 3] + q) and z(second[2, 2]) \
        and all(z(second[i, j]) for i in range(4) for j in range(4) if i != j)
    out['T3_no_higher_orders'] = all(z(sp.expand(Ric[i, j]).coeff(eps, k)) for i in range(4) for j in range(4) for k in (3, 4))
    # T4: flat slices, one time, ANY fall velocity u(x): the connection of the frame e0 = d_t + u.grad, e_i = d_i
    T, x, y, zz = sp.symbols('T x y z', real=True)
    co4 = (T, x, y, zz)
    u = [sp.Function('u%d' % i)(x, y, zz) for i in (1, 2, 3)]
    E = sp.Matrix([[1, u[0], u[1], u[2]], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    c = tp.structure_in(co4, E)
    w = cv.connection(c)
    d = lambda i, j: sp.diff(u[i-1], co4[j])                       # d_j u_i
    # a direction S carried with the frame: dS^i = -w^i_{j0} S^j = (Omega x S)^i with Omega = (1/2) curl u
    cu = [d(3, 2) - d(2, 3), d(1, 3) - d(3, 1), d(2, 1) - d(1, 2)]
    S = sp.symbols('S1:4')
    carried = [-sum(w[i][j][0]*S[j-1] for j in (1, 2, 3)) for i in (1, 2, 3)]
    cross = sp.Matrix(cu).cross(sp.Matrix(S))/2
    out['T4_turn_is_half_curl'] = all(z(carried[i] - cross[i]) for i in range(3))
    out['T4_frame_in_free_fall'] = all(z(w[i][0][0]) for i in range(4))
    out['T4_no_turn_along_the_slice'] = all(z(w[i][j][k]) for i in (1, 2, 3) for j in (1, 2, 3) for k in (1, 2, 3))
    out['T4_stretch_is_symmetric_part'] = all(z(w[i][0][j] - (d(i, j) + d(j, i))/2) for i in (1, 2, 3) for j in (1, 2, 3))
    # T5: the fall velocity of the turning centre, u = -beta rhat + (a/r^3) zhat x r : half its curl
    rr = sp.sqrt(x**2 + y**2 + zz**2)
    pos = sp.Matrix([x, y, zz]); zh = sp.Matrix([0, 0, 1])
    uu = -sp.sqrt(rs/rr)*pos/rr + (a/rr**3)*zh.cross(pos)
    curl = sp.Matrix([sp.diff(uu[2], y) - sp.diff(uu[1], zz), sp.diff(uu[0], zz) - sp.diff(uu[2], x), sp.diff(uu[1], x) - sp.diff(uu[0], y)])
    target = a*(3*(zz/rr)*pos/rr - zh)/rr**3
    out['T5_curl'] = all(z(v) for v in (curl - target))
    # T6: mean over a circle through the poles and over the circle round the axis
    psi, rho = sp.symbols('psi rho', positive=True)
    polar = (curl/2).subs({x: rho*sp.sin(psi), y: 0, zz: rho*sp.cos(psi)})
    mean_polar = polar.applyfunc(lambda e: sp.simplify(sp.integrate(sp.simplify(e), (psi, 0, 2*sp.pi))/(2*sp.pi)))
    equat = (curl/2).subs({x: rho*sp.cos(psi), y: rho*sp.sin(psi), zz: 0}).applyfunc(sp.simplify)
    out['T6_mean_over_the_poles'] = all(z(v) for v in (mean_polar - sp.Matrix([0, 0, a/(4*rho**3)])))
    out['T6_round_the_axis'] = all(z(v) for v in (equat - sp.Matrix([0, 0, -a/(2*rho**3)])))
    # numbers: the Earth, circuit over the poles at 7020 km; a = 2 G J / c^2 is put in
    GM, c_, R_e, w_e, k_e, rad = 3.986004418e14, 299792458.0, 6.3781e6, 7.292115e-5, 0.3307, 7.020e6
    a_e = 2*k_e*GM*R_e**2*w_e/c_**2                                 # m^3 / s
    rate = a_e/(4*rad**3)                                           # rad / s
    mas = rate*31557600.0*206264806.247
    num = dict(a_m3_per_s=a_e, mean_rate_mas_per_year=round(mas, 2), projected_cos_16p84deg=round(mas*0.95712, 2))
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, numbers=num), open(os.path.join(HERE, 'SW1_RESULT.json'), 'w'), indent=1)
    return out, num


if __name__ == '__main__':
    print(run())
