"""CS1: why the memory of the stretch equals the energy of the source: the source enters through its count.

OR1: a history with clock g contributes its count g * integral of e^0 (CL2-T1).  The frame law is TP1's Q.
Let the time step of the frame vary: e_0 = (1/N)(d_t + u.grad), e_i = d_i (N = 1 is SE1's family).
Stationarity of  (1/2k) * integral of N Q  -  g nu * integral of N   with respect to N, at N = 1.
Symbolic (sympy), with the functions of TP1.  Python 3.12."""
import json, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cv1'))
sys.path.insert(0, os.path.join(HERE, '..', 'tp1'))
_cwd = os.getcwd(); os.chdir(os.path.join(HERE, '..', 'tp1'))
import cv1_frame_curvature as cv
import tp1_frame_defect_law as tp
os.chdir(_cwd)
z = lambda e: sp.simplify(e) == 0


def run():
    out = {}
    T, x, y, zz = sp.symbols('T x y z', real=True); co = (T, x, y, zz)
    N = sp.Function('N')(x, y, zz)
    u = [sp.Function('u%d' % i)(x, y, zz) for i in (1, 2, 3)]
    E = sp.Matrix([[1/N, u[0]/N, u[1]/N, u[2]/N], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    c = tp.structure_in(co, E)
    I1, I2, I3 = tp.invariants(c)
    Q = sp.Rational(1, 4)*I1 + sp.Rational(1, 2)*I2 - I3
    L = sp.expand(N*Q)                                            # N is the volume factor of the frame
    d_ = lambda i, j: sp.diff(u[i], co[j + 1])
    K = sp.Matrix(3, 3, lambda i, j: (d_(i, j) + d_(j, i))/2)
    e2 = sum(K[i, i]*K[j, j] - K[i, j]**2 for i in range(3) for j in range(i + 1, 3))
    # C1: with the time step in, the law is  -2 e2(K)/N  plus a pure boundary term
    grad = [sp.diff(N, v) for v in (x, y, zz)]
    rest = sp.simplify(L + 2*e2/N)
    # Euler expression of the rest with respect to every field vanishes  <=>  boundary term
    def euler(expr, f):
        vars_ = (x, y, zz)
        val = sp.diff(expr, f)
        for v in vars_:
            val -= sp.diff(sp.diff(expr, sp.diff(f, v)), v)
        for v in vars_:
            for w_ in vars_:
                val += sp.diff(sp.diff(expr, sp.diff(f, v, w_)), v, w_)*(sp.Rational(1, 2) if v != w_ else 1)
        return sp.simplify(val)
    out['C1_law_with_time_step'] = all(z(euler(rest, f)) for f in [N] + u)
    # C2: response of the law to the time step, at N = 1:  d(N Q)/dN = + 2 e2(K)
    at1 = lambda e: sp.simplify(e.subs({sp.diff(N, v): 0 for v in (x, y, zz)}).subs(N, 1).doit())
    resp = at1(euler(L, N).subs({sp.diff(N, v, w_): 0 for v in (x, y, zz) for w_ in (x, y, zz)}))
    out['C2_response_to_the_time_step'] = z(resp - 2*e2)
    # C3: a source that falls with the frame, nu histories per volume, clock g each: count per frame volume g nu N.
    #     Stationarity of (1/2k) N Q - g nu N in N at N = 1:   e2(K) = k g nu
    k, g, nu = sp.symbols('k g nu', positive=True)
    total = L/(2*k) - g*nu*N
    stat = at1(euler(total, N).subs({sp.diff(N, v, w_): 0 for v in (x, y, zz) for w_ in (x, y, zz)}))
    out['C3_memory_equals_count_density'] = z(stat - (e2/k - g*nu))
    # C4: response to the fall velocity at N = 1: the (time, cut) law of SE1, curl curl u = 0 ; the source at rest in the frame adds nothing
    curl = [d_(2, 1) - d_(1, 2), d_(0, 2) - d_(2, 0), d_(1, 0) - d_(0, 1)]
    cc = [sp.diff(curl[2], y) - sp.diff(curl[1], zz), sp.diff(curl[0], zz) - sp.diff(curl[2], x), sp.diff(curl[1], x) - sp.diff(curl[0], y)]
    L1 = L.subs({sp.diff(N, v): 0 for v in (x, y, zz)}).subs(N, 1).doit()
    ru = [euler(L1, f) for f in u]
    out['C4_response_to_the_fall_velocity'] = all(z(ru[i] - cc[i]) for i in range(3)) or all(z(ru[i] + cc[i]) for i in range(3))
    # C5: a centre: e2 = (r m)'/r^2 for radial fall (SE1-S4), so the total count inside a sphere fixes r_s:  r m = k g * (number inside) / (4 pi)
    r = sp.Symbol('r', positive=True); m = sp.Function('m')(r)
    beta = sp.sqrt(m)
    e2r = (beta/r)**2*(2*r*sp.diff(beta, r)/beta + 1)
    out['C5_radial_memory'] = z(e2r - sp.diff(r*m, r)/r**2)
    # a ball of histories, n0 per volume, radius a:  inside m = k g n0 r^2/3, outside m = r_s/r with r_s = k g n0 a^3/3 = k g (number)/(4 pi)
    a_, n0 = sp.symbols('a n0', positive=True)
    m_in = k*g*n0*r**2/3
    rs = k*g*n0*a_**3/3
    e2_of = lambda mm: sp.simplify(sp.diff(r*mm, r)/r**2)
    out['C5_inside'] = z(e2_of(m_in) - k*g*n0)
    out['C5_outside_and_joining'] = z(e2_of(rs/r)) and z(m_in.subs(r, a_) - rs/a_) and z(rs - k*g*(n0*sp.Rational(4, 3)*sp.pi*a_**3)/(4*sp.pi))
    # inside, the fall speed is beta = h r with h^2 = k g n0/3, every rate equal: e2 = 3 h^2 and Q = -6 h^2 (CO1's value)
    h = sp.Symbol('h', positive=True)
    Kh = sp.eye(3)*h
    e2h = sum(Kh[i, i]*Kh[j, j] for i in range(3) for j in range(i + 1, 3))
    out['C5_uniform_is_CO1'] = z(e2h - 3*h**2) and z(sp.sqrt(m_in)/r - sp.sqrt(k*g*n0/3))
    out = {k_: bool(v) for k_, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out), open(os.path.join(HERE, 'CS1_RESULT.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(run())
