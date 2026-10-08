"""CZ1: the centre record of a closed surface: dual couplings add at fixed cosets; after the cosets, a factor d^chi;
the weight of a twist is a power for one face and an exponential of the turn coupling for a cube.

CT1-C5: on faces the centre marks form a two-valued gauge record with local coupling kappa c_p.
theorum/41 (CG1): the cut here is the centre flip; even channel cosh, odd channel sinh.
Exact quaternion algebra with the second moments of the sphere, an exact gluing engine, exact enumeration over centre
marks, mpmath for the Bessel sums.  Python 3.12."""
import itertools, json, math, os, random
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------- quaternions (the turn block as the unit sphere, YM-F1)
def qmul(a, b):
    return (a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3], a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
            a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1], a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0])


def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])


def sphere_mean_degree2(expr, u):
    """mean over the unit sphere of a polynomial of degree two in u: <u_a u_b> = delta_ab / 4 (YM-37: rational moments)."""
    expr = sp.expand(expr)
    out = 0
    for a_ in range(4):
        for b_ in range(4):
            co = expr.coeff(u[a_], 2) if a_ == b_ else sp.expand(expr).coeff(u[a_], 1).coeff(u[b_], 1)
            if a_ == b_:
                out += co*sp.Rational(1, 4)
    rest = expr.subs({x: 0 for x in u})
    return sp.expand(out + rest)


# ---------- gluing engine: faces as cyclic words of (link, +-1); returns the exponent e with <prod chi(U_p)> = d^e
def glue(faces):
    words = [list(f) for f in faces]
    links = sorted({l for f in faces for (l, s) in f})
    power = 0
    for l in links:
        hits = [(i, j) for i, w in enumerate(words) for j, (ll, s) in enumerate(w) if ll == l]
        assert len(hits) == 2
        (i1, j1), (i2, j2) = hits
        assert words[i1][j1][1] == -words[i2][j2][1]
        if i1 != i2:                                   # chi(A U) chi(U^-1 B) -> chi(A B)/d
            w1 = words[i1][j1 + 1:] + words[i1][:j1]
            w2 = words[i2][j2 + 1:] + words[i2][:j2]
            new = w1 + w2
            words = [w for k_, w in enumerate(words) if k_ not in (i1, i2)] + [new]
        else:                                          # chi(A U B U^-1) -> chi(A) chi(B)/d
            w = words[i1]
            lo, hi = sorted((j1, j2))
            inner, outer = w[lo + 1:hi], w[hi + 1:] + w[:lo]
            words = [x for k_, x in enumerate(words) if k_ != i1] + [inner, outer]
        power -= 1
    assert all(len(w) == 0 for w in words)
    return power + len(words)                          # each empty word is chi(1) = d


def box_surface(nx, ny, nz):
    """outward-oriented unit faces of the surface of an nx x ny x nz box; links named by their end points."""
    def link(p, q):
        return ((p, q), 1) if p < q else ((q, p), -1)
    faces = []
    def quad(a, b, c, d):
        faces.append([link(a, b), link(b, c), link(c, d), link(d, a)])
    for x in range(nx):
        for y in range(ny):
            quad((x, y, 0), (x, y + 1, 0), (x + 1, y + 1, 0), (x + 1, y, 0))
            quad((x, y, nz), (x + 1, y, nz), (x + 1, y + 1, nz), (x, y + 1, nz))
    for x in range(nx):
        for z in range(nz):
            quad((x, 0, z), (x + 1, 0, z), (x + 1, 0, z + 1), (x, 0, z + 1))
            quad((x, ny, z), (x, ny, z + 1), (x + 1, ny, z + 1), (x + 1, ny, z))
    for y in range(ny):
        for z in range(nz):
            quad((0, y, z), (0, y, z + 1), (0, y + 1, z + 1), (0, y + 1, z))
            quad((nx, y, z), (nx, y + 1, z), (nx, y + 1, z + 1), (nx, y, z + 1))
    return faces


def torus_surface(L):
    def link(p, q, name):
        return (name, 1)
    faces = []
    for x in range(L):
        for y in range(L):
            h = lambda xx, yy: ('h', xx % L, yy % L)     # link from (x,y) to (x+1,y)
            v = lambda xx, yy: ('v', xx % L, yy % L)     # link from (x,y) to (x,y+1)
            faces.append([(h(x, y), 1), (v(x + 1, y), 1), (h(x, y + 1), -1), (v(x, y), -1)])
    return faces


def run():
    out, num = {}, {}
    zz = lambda e: sp.simplify(e) == 0

    # Z1: the two gluing rules for the first content, from the second moments of the sphere
    A = sp.symbols('a0:4'); B = sp.symbols('b0:4'); u = sp.symbols('u0:4')
    c = lambda q: q[0]
    r1 = sphere_mean_degree2(c(qmul(A, u))*c(qmul(qconj(u), B)), u)
    out['Z1_join_two_faces_along_a_link'] = zz(r1 - c(qmul(A, B))/4)
    r2 = sphere_mean_degree2(c(qmul(qmul(qmul(A, u), B), qconj(u))), u)
    out['Z1_link_twice_in_one_face'] = zz(r2 - c(A)*c(B)*sum(x*x for x in u).subs({x: sp.Rational(1, 2) for x in u}))   # |u|^2 = 1
    # (in characters chi = 2c:  chi chi -> chi/2  and  chi -> chi chi/2 ; for content j the same with d_j, YM-3/4)

    # Z2: closed surfaces: <prod chi(U_p)> = d^(V - E) ; with the face weights the surface carries d^chi (f_j/f_0)^F
    cases = {'cube': (box_surface(1, 1, 1), 2, 6), 'box 2x1x1': (box_surface(2, 1, 1), 2, 10), 'box 2x2x1': (box_surface(2, 2, 1), 2, 16),
             'torus 2x2': (torus_surface(2), 0, 4), 'torus 3x3': (torus_surface(3), 0, 9)}
    ok = True
    for name, (faces, chi, F) in cases.items():
        ok &= glue(faces) == chi - F
    out['Z2_closed_surface_carries_d_to_the_euler_number'] = ok

    # Z3: at fixed cosets the centre sum over the twelve links of a cube is 2^12 (prod cosh + prod sinh): dual couplings add
    random.seed(3)
    edges = sorted({l for f in box_surface(1, 1, 1) for (l, s) in f})
    faces = [[edges.index(l) for (l, s) in f] for f in box_surface(1, 1, 1)]
    ev = [Fr(random.randint(2, 9), random.randint(1, 5)) for _ in range(6)]          # e^(k_p) > 0, exact
    tot = Fr(0)
    for signs in itertools.product((1, -1), repeat=12):
        w = Fr(1)
        for f, e in zip(faces, ev):
            sgn = math.prod(signs[i] for i in f)
            w *= e if sgn == 1 else 1/e
        tot += w
    ch = [(e + 1/e)/2 for e in ev]; sh = [(e - 1/e)/2 for e in ev]
    out['Z3_centre_sum_at_fixed_cosets'] = tot == 2**12*(math.prod(ch) + math.prod(sh))
    k = sp.symbols('k1:7', positive=True)
    ksum = sum(-sp.log(sp.tanh(x))/2 for x in k)
    out['Z3_dual_couplings_add_over_the_surface'] = all(abs(float((sp.prod([sp.tanh(x) for x in k]) - sp.exp(-2*ksum)).subs(dict(zip(k, vals))))) < 1e-14 for vals in ((0.3, 0.5, 0.7, 1.1, 1.3, 2.0), (1, 1, 1, 1, 1, 1)))

    # Z4: after the cosets.  Z_cube = sum_j d_j^2 f_j^6, f_j = 2 I_(2j+1)(kappa)/kappa: even contents n odd, odd contents n even (n = 2j+1)
    mp.mp.dps = 40
    def parts(v):
        N = int(4*v) + 60
        b = [mp.besseli(n, v)**6*n*n for n in range(1, 2*N + 2)]
        return mp.fsum(b[0::2]), mp.fsum(b[1::2])
    r = lambda v: mp.besseli(2, v)/mp.besseli(1, v)
    # strong end: Z_o/Z_e -> 4 r^6, not r^6
    oks = True
    for v in (mp.mpf(1)/10, mp.mpf(1)/2, mp.mpf(1)):
        ze, zo = parts(v)
        oks &= abs((zo/ze)/(4*r(v)**6) - 1) < 10*r(v)**6 + 20*(mp.besseli(3, v)/mp.besseli(1, v))**6
    out['Z4_strong_end_factor_four'] = bool(oks)
    # weight of a twist F = (Z_e - Z_o)/(Z_e + Z_o): one face flipped
    twist = lambda v: (lambda p: (p[0] - p[1])/(p[0] + p[1]))(parts(v))
    face_twist = lambda v: (1 - r(v))/(1 + r(v))
    out['Z4_face_twist_is_a_power'] = all(abs(face_twist(v)*4*v/3 - 1) < 2/v for v in (mp.mpf(50), mp.mpf(500)))
    # closed surface of F faces (sphere): Z = sum_j d_j^2 f_j^F.  The twist weight falls as exp(-E kappa) with
    #     E = F (1 - cos(pi/F)): every face turned by one F-th of a full turn (classical value; exponent checked numerically)
    def twistF(v, F):
        E0 = F*(1 - math.cos(math.pi/F))
        mp.mp.dps = int(E0*float(v)/2.3) + 40
        v = mp.mpf(v)
        N = int(math.sqrt(2*float(v)*mp.mp.dps*2.31/F)) + 12            # beyond this the terms are below the working precision
        b = [mp.besseli(n, v)**F*n*n for n in range(1, 2*N + 2)]
        ze, zo = mp.fsum(b[0::2]), mp.fsum(b[1::2])
        return (ze - zo)/(ze + zo)
    expo, oke = {}, True
    for F in (2, 4, 6, 10, 16):
        E0 = F*(1 - math.cos(math.pi/F))
        sl = lambda v: -(mp.log(twistF(v + 1, F)) - mp.log(twistF(v, F)))
        s1, s2 = sl(40), sl(80)
        est = float(2*s2 - s1)                              # removes the 1/kappa term of the slope
        expo[F] = dict(classical=E0, from_the_sums=est)
        oke &= abs(est - E0) < 3e-3
    mp.mp.dps = 40
    num['twist_exponent_by_number_of_faces'] = expo
    out['Z4_twist_weight_is_exponential_with_the_classical_exponent'] = oke
    out['Z4_cube_exponent_is_6_minus_3_root_3'] = abs(6*(1 - math.cos(math.pi/6)) - (6 - 3*math.sqrt(3))) < 1e-14 and abs(4*(1 - math.cos(math.pi/4)) - (4 - 2*math.sqrt(2))) < 1e-14
    # the classical configuration: F equal turns by 2 pi / F about one axis multiply to -1
    okc = True
    for F in (2, 3, 4, 6, 12):                              # turns about one axis: the pair (cos, sin) of the half angle pi/F
        turn = sp.exp(sp.I*sp.pi/F)
        okc &= sp.simplify(turn**F + 1) == 0 and sp.re(turn) == sp.cos(sp.pi/F)
    out['Z4_classical_twist'] = okc
    FF = sp.symbols('F', positive=True)
    out['Z4_large_surfaces_cost_less'] = sp.limit(FF*FF*(1 - sp.cos(sp.pi/FF)), FF, sp.oo) == sp.pi**2/2
    num['cube_twist'] = {v: float(twistF(v, 6)) for v in (5, 10, 20, 40, 80)}
    mp.mp.dps = 40
    num['cube_twist_times_exp_over_kappa'] = {v: float(twistF(v, 6)*mp.e**((6 - 3*mp.sqrt(3))*v)/v) for v in (40, 80, 160)}
    mp.mp.dps = 40
    twist = lambda v: twistF(v, 6)
    Kc = lambda v: -mp.log(twist(v))/2
    num['cube_centre_coupling_over_kappa'] = {v: float(Kc(v)/v) for v in (20, 40, 80)}
    num['half_the_exponent'] = (6 - 3*math.sqrt(3))/2
    mp.mp.dps = 40
    # where the cube is self-dual, and the size of number for a twist weight of 7.7e-20
    kd = mp.findroot(lambda v: (lambda p: p[1]/p[0])(parts(v)) - (mp.sqrt(2) - 1), 4)
    num['turn_coupling_where_the_cube_is_self_dual'] = float(kd)
    k20 = mp.findroot(lambda v: mp.log(twist(v)) - mp.log(mp.mpf('7.685e-20')), 60)
    num['turn_coupling_for_a_twist_weight_of_7.7e-20'] = float(k20)
    out['Z4_numbers'] = 3 < kd < 6 and 40 < k20 < 80

    # Z5: Hamiltonian form (theorum/41 iv): cut = centre flip of one link; generator sum Delta + theta sum W.
    #     Delta of a face function: sum over its four links of sum_a D_a^2, D_a f(U) = d/dt f(U exp(t e_a/2)):  (e_a/2)^2 = -1/4
    e_units = [(0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    lap = sum(c(qmul(qmul(A, qmul(e, e)), B)) for e in e_units)/4
    out['Z5_laplacian_of_a_face_per_link'] = zz(lap + sp.Rational(3, 4)*c(qmul(A, B)))
    num['seam_curvature_on_the_free_vacuum'] = '-3 theta x (faces through the link)'

    out['pass'] = all(out.values())
    json.dump({'checks': out, 'numbers': num}, open(os.path.join(HERE, 'CZ1_RESULT.json'), 'w'), indent=1, default=str)
    return out, num


if __name__ == '__main__':
    o, n = run()
    for k_, v in o.items(): print(k_, v)
    for k_, v in n.items(): print(k_, v)
